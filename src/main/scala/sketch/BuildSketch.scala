package sketch

import java.io.{File, PrintWriter}

import atk.ProgressBar.progress
import bloomfilter.CanGenerateHashFrom
import bloomfilter.mutable.BloomFilter
import utilities.FileHandling.{writeSerialized, timeStamp, verifyDirectory, verifyFile}
import utilities.SequenceUtils._
import utilities.SequenceFormatUtils.loadSequenceFile
import utilities.SketchUtils.RedwoodSketch

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
    val parser = new scopt.OptionParser[Config]("build-sketch") {
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
      opt[Int]("log") hidden() action { (x, c) =>
        c.copy(log = x)
      } text ("Count for progress logger.")
      note("\nOPTIONAL: HASH WEIGHTS (CNV)\n")
      opt[Unit]("track-cov") hidden() action { (x, c) =>
        c.copy(trackCov = true)
      } text ("Keep track of coverage of kmers that are inside the sketch.")
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
      + config.sampleName + " using " + config.readFile.size + " read files")
    //create sketch
    var sketch = new mutable.PriorityQueue[Int]()
    //set mutable sketch size
    var current_sketch_size = 0

    /**
      * Function to determine whether a given hash is less than the current sketch's max. Automatically true if
      * sketch is empty.
      *
      * @return Boolean
      */
    def isLessThanMax: Int => Boolean = hash => (current_sketch_size < config.sketchSize) || (hash <= sketch.head)

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
                    //in sketch, update coverage
                    if (in_sketch_cov.nonEmpty) {
                      if (config.trackCov) inSketch.update(hash, in_sketch_cov.get + 1)
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
                          sketch.enqueue(hash)
                          //increment sketch size
                          current_sketch_size += 1
                        }
                        //sketch is full
                        else {
                          //remove current max from insketch
                          inSketch.remove(sketch.head)
                          //remove current max from sketch
                          sketch.dequeue()
                          //remove new max from candidates
                          candidateHashes.remove(hash)
                          //add new max to insketch
                          inSketch.update(hash, true_cov)
                          //add new max to sketch
                          sketch.enqueue(hash)
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
      //load hashes from sketch and their coverage
      val tmp = inSketch.toMap
      //sanity check
      assert(tmp.keySet == sketch.toSet)
      if (!config.trackCov) tmp.mapValues(_ => 1.0).map(identity)
      else {
        println(timeStamp + "Estimating kmer copy number")
        println(timeStamp + "--Removing kmers with 0.01 lowest/highest coverage")
        val drop_size = (current_sketch_size * 0.01).toInt
        val kmers_per_cov = tmp.values.toList.sorted //.drop(drop_size).dropRight(drop_size)
          .groupBy(identity).mapValues(_.size).toList.sortBy(_._1)
        val pw2 = new PrintWriter(config.outputDir + "/" + config.sampleName + ".cnvs.txt")
        kmers_per_cov.foreach(x => pw2.println(x._1 + "\t" + x._2))
        pw2.close
        //normalize weights of original kmers
        tmp.mapValues(_ => 1.0).map(identity)
      }
    }
    //clear out sketch, bloom filter, hash-tables
    bf.dispose()
    sketch.clear()
    inSketch.clear()
    candidateHashes.clear()
    //create redwood sketch object
    println(timeStamp + "Writing to disk")
    val redwood_sketch = new RedwoodSketch(config.sampleName, config.kmerSize, sketch_kmers)
    writeSerialized(redwood_sketch, new File(config.outputDir + "/" + config.sampleName + ".rdws"))
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

