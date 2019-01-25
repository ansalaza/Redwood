package sketch

import java.io.File

import atk.ProgressBar.progress
import bloomfilter.CanGenerateHashFrom
import bloomfilter.mutable.BloomFilter
import boopickle.Default._
import utilities.FileHandling.{timeStamp, verifyDirectory, verifyFile, writeSerialized}
import utilities.SequenceUtils._
import utilities.SequenceFormatUtils.loadSequenceFile
import utilities.SketchUtils.RedwoodSketch
import utilities.NumericalUtils.{mean, stdDev}

import scala.collection.mutable
import scala.util.hashing.MurmurHash3
import utilities.SequenceUtils.ByteEncoded

import scala.annotation.tailrec

/**
  * Author: Alex N. Salazar
  * Created on 3-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object BuildSketch {

  case class Config(
                     readFile: Seq[File] = null,
                     kmerSize: Int = 21,
                     sampleName: String = null,
                     sketchSize: Int = 100000,
                     log: Int = 1000000,
                     trackCov: Boolean = false,
                     isAssembly: Boolean = false,
                     minCov: Int = 2,
                     outputDir: File = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("sketch") {
      opt[Seq[File]]('r', "read-files") valueName ("<file1>,<file2>,...") required() action { (x, c) =>
        c.copy(readFile = x)
      } text ("Read or assembly file(s) in FASTA, FASTQ, or FASTQ.GZ format.")
      opt[String]("name") required() action { (x, c) =>
        c.copy(sampleName = x)
      } text ("Sample name.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      note("\nOPTIONAL: SKETCH CONSTRUCTION\n")
      opt[Unit]("assembly") action { (x, c) =>
        c.copy(isAssembly = true)
      } text ("Input file is a genome assembly. Consider all kmers (e.g. set 'min-cov' to 1).")
      opt[Int]("min-cov") action { (x, c) =>
        c.copy(minCov = x)
      } text ("Minimum coverage of kmer before adding to sketch (default is 2).")
      opt[Int]("kmer-size") action { (x, c) =>
        c.copy(kmerSize = x)
      } text ("Size of kmer (default 21).")
      opt[Int]("sketch-size") action { (x, c) =>
        c.copy(sketchSize = x)
      } text ("Size of sketch (default is 100000).")
      opt[Unit]("track-cov") action {(x,c) =>
        c.copy(trackCov = true)
      } text("Track coverage of kmers in sketch")
      opt[Int]("log") hidden() action { (x, c) =>
        c.copy(log = x)
      } text ("Count for progress logger.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      config.readFile.foreach(x => verifyFile(x))
      buildSketch(config)
    }
  }

  def buildSketch(config: Config): Unit = {
    println(timeStamp + "Building sketch of " + config.sketchSize + " kmers of size " + config.kmerSize + " for "
      + config.sampleName + " using " + config.readFile.size + " sequence files")
    //create sketch
    var sketch = new mutable.PriorityQueue[(Int, ByteEncoded)]()(Ordering.by(_._1))
    //set mutable sketch size
    var current_sketch_size = 0

    /**
      * Function to determine whether a given hash is less than the current sketch's max. Automatically true if
      * sketch is empty.
      *
      * @return Boolean
      */
    def isLessThanMax: Int => Boolean = hash => (current_sketch_size < config.sketchSize) || (hash <= sketch.head._1)

    //hashmap for keeping track of frequency of hashes in sketch
    var inSketch = mutable.HashMap[Int, Int]()
    //hashmap for keeping track of frequency of candidate hashes
    var candidateHashes = mutable.HashMap[Int, Int]()

    /**
      * Function to get the count of a given candidate hash
      *
      * @return Int
      */
    def getCandidateCount: Int => Int = hash => if (config.minCov == 2) 1 else candidateHashes.getOrElse(hash, 1)

    //parameters for bloom filter
    val expectedElements = 100000000
    val falsePositiveRate = 0.1

    /**
      * Override generateHash to avoid double hashing
      */
    implicit object CanGenerateFromLong extends CanGenerateHashFrom[Long] {
      override def generateHash(from: Long): Long = from
    }

    //create bloom filter
    val bf = BloomFilter[Long](expectedElements, falsePositiveRate)

    /**
      * Function to determine whether a given hash is in the bloom filter considering input may be an assembly
      * @return Boolean
      */
    def isInBF: Int => Boolean = hash => config.isAssembly || bf.mightContain(hash)

    //set min coverage
    val min_cov = if(config.isAssembly) 1 else config.minCov

    //iterate through each read file and create sketch (in parallel if multiple threads specified
    config.readFile.foreach(read_file => {
      //load iterator
      val iterator = loadSequenceFile(read_file)
      //go through reads and construct sketch
      while (iterator.hasNext()) {
        //byte-encode sequence
        val seq = iterator.next()
        //get size of sequence
        val seq_length = seq.size
        //check sequence at least as big as the specified kmer
        if (seq_length >= config.kmerSize) {
          //obtain reverse complement
          val seq_reverse = reverseComplement(seq)
          //define stopping index
          val size_stop = (seq_length - config.kmerSize) + 1
          //mutable index variable
          var index = 0
          //iterate until reach
          while (index < size_stop) {
            //get forward kmer
            val forward = seq.slice(index, index + config.kmerSize)
            //only if forward is valid
            if (isValid(forward)) {
              //get reverse kmer
              val reverse = seq_reverse.slice(seq_length - config.kmerSize - index, seq_length - index)
              //only if reverse is valid also
              if (isValid(reverse)) {
                //get smallest kmer
                val smallest = smallestKmer(forward, reverse, config.kmerSize)
                //get hash
                val hash = MurmurHash3.bytesHash(smallest)
                //convert to long (for utilizing bloomfilter)
                val hashL = hash.toLong
                //sketch is less than current max, worth considering
                if (isLessThanMax(hash)) {
                  //not an assembly file and kmer has not been observed before
                  if (!isInBF(hash)) bf.add(hashL)
                  //kmer has been observed before
                  else {
                    //get in-sketch coverage, if it's in the sketch
                    val in_sketch_cov = inSketch.get(hash)
                    //in sketch
                    if (in_sketch_cov.nonEmpty) {
                      //update coverage if specified
                      if(config.trackCov) inSketch.update(hash, in_sketch_cov.get + 1)
                    }
                    //not in current sketch
                    else {
                      //get coverage of current hash
                      val true_cov = getCandidateCount(hash) + 1
                      //count does not meet min threshold, update
                      if (true_cov < min_cov) candidateHashes.update(hash, true_cov)
                      //coverage meets min threshold
                      else {
                        //sketch is not full
                        if (current_sketch_size < config.sketchSize) {
                          //remove from candidates
                          candidateHashes.remove(hash)
                          //add to insketch
                          inSketch.update(hash, true_cov)
                          //add to sketch
                          sketch.enqueue((hash, smallest))
                          //increment sketch size
                          current_sketch_size += 1
                        }
                        //sketch is full
                        else {
                          //remove current max from insketch
                          inSketch.remove(sketch.head._1)
                          //remove current max from sketch
                          sketch.dequeue()
                          //remove new max from candidates
                          candidateHashes.remove(hash)
                          //add new max to insketch
                          inSketch.update(hash, true_cov)
                          //add new max to sketch
                          sketch.enqueue((hash,smallest))
                        }
                      }
                    }
                  }
                }
              }
            }
            //update slicing index
            index += 1
          }
        }
        //log progress
        progress(config.log)
      }
    })
    println(timeStamp + "Completed sketch of size " + sketch.size)
    //get kmers in sketch
    val sketch_kmers = {
      //set sketch map
      val sketch_map = {
        //convert to map
        val tmp = sketch.toMap
        //sanity check
        assert(tmp.keySet == inSketch.keySet)
        //set coverage and kmer sequence
        tmp.map(x => (x._1, ((if(!config.trackCov) 0 else inSketch(x._1)), x._2)))
      }
      if(config.trackCov){
        //get mean coverage
        val mean_cov = mean(sketch_map.map(_._2._1))
        //get standard deviation
        val std_cov = stdDev(sketch_map.map(_._2._1))
        println(timeStamp + "Mean kmer coverage of " + mean_cov + " with standard deviation of " + std_cov)
      }
      sketch_map.map(identity)
    }
    //clear out sketch, bloom filter, hash-tables
    bf.dispose()
    sketch.clear()
    inSketch.clear()
    candidateHashes.clear()
    //create redwood sketch object
    println(timeStamp + "Writing to disk")
    val redwood_sketch = new RedwoodSketch(config.sampleName, config.kmerSize, sketch_kmers)
    //serialize and write to disk
    writeSerialized(
      Pickle.intoBytes(redwood_sketch).array(),
      new File(config.outputDir + "/" + config.sampleName + ".rdws"))
    println(timeStamp + "Successfully completed!")
  }

  /**
    * Function compare a given byte-encoded DNA sequence(kmers) forward and revers orientation and return the
    * smallest lexicgoraphic version
    *
    * @return DNAseq
    */
  def smallestKmer(forward: ByteEncoded, reverse: ByteEncoded, kmer_size: Int): ByteEncoded = {
    /**
      * Tail-recursive method
      *
      * @param i
      * @return
      */
    @tailrec def _smallestKmer(i: Int): ByteEncoded = {
      if (i == kmer_size) forward
      else {
        if (forward(i) < reverse(i)) forward
        else if (forward(i) > reverse(i)) reverse
        else _smallestKmer(i + 1)
      }
    }

    _smallestKmer(0)
  }
}

