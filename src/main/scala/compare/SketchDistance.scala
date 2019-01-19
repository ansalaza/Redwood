package compare

import java.io.{File, PrintWriter}

import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory, verifyFile}
import utilities.SketchUtils.loadRedwoodSketch
import scala.math.log
import scala.annotation.tailrec

import utilities.NumericalUtils.min

/**
  * Author: Alex N. Salazar
  * Created on 10-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object SketchDistance {

  case class Config(
                     sketchFiles: Seq[File] = null,
                     pathsFile: File = null,
                     outputDir: File = null,
                     prefix: String = null,
                     threads: Int = 1,
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("sketch-converter") {
      opt[Seq[File]]('s', "sketch-files") valueName ("stetch1,sketch2,...") action { (x, c) =>
        c.copy(sketchFiles = x)
      } text ("Comma-separated sketch files to compare. Performs-pairwise mash distance for all provided sketches.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      opt[String]("prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output matrix file.")
      note("\nOPTIONAL\n")
      opt[File]("paths-file") action { (x, c) =>
        c.copy(pathsFile = x)
      } text ("Alternatively, provide a file with the paths of all sketches to compare, one per line.")
      opt[Int]('t', "threads") action { (x, c) =>
        c.copy(threads = x)
      } text ("Max number of threads..")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists
      verifyDirectory(config.outputDir)
      //get list of sketches from input parameter
      val all_sketches = {
        if (config.sketchFiles != null && config.pathsFile != null) {
          assert(false, "Provide sketch files through '--sketch-files' or '--paths-file'")
          List[File]()
        }
        else if (config.sketchFiles != null) {
          config.sketchFiles.foreach(verifyFile(_))
          config.sketchFiles.toList
        }
        else {
          val sketches = openFileWithIterator(config.pathsFile).toList.map(new File(_))
          sketches.foreach(verifyFile(_))
          sketches
        }
      }
      assert(config.threads > 0, "Number of threads specified is a non-positive integer")
      sketchDistance(all_sketches, config)
    }
  }

  def sketchDistance(sketches: List[File], config: Config): Unit = {
    println(timeStamp + "Found " + sketches.size + " sketches")
    //obtain sketch size to use
    val (sketch_size, kmer_size, names) = {
      //load sketches into 3-tuples: (name, sketch size, kmer size)
      val tmp = sketches.map(sketch => {
        val s = loadRedwoodSketch(sketch);
        (s.name, s.kmers.size, s.kmerSize)
      })
      //group by kmer size
      val kmer_sizes = tmp.groupBy(_._3)
      //sanity check
      assert(kmer_sizes.size == 1, "One more sketches have different kmer size:\n" + kmer_sizes.toList.mkString("\n"))
      //group by name, get counts for each name
      val all_names = tmp.groupBy(_._1).mapValues(_.size)
      val non_unique = all_names.filter(_._2 > 1)
      //sanity check
      assert(non_unique.isEmpty, "The following sketch names are not unique: " + all_names.keys.mkString(","))
      //get min sketch size
      (tmp.map(_._2).min, kmer_sizes.head._1, all_names.keys.toList)
    }
    println(timeStamp + "Coputing pairwise-distances using sketch size of " + sketch_size + " kmers of size " + kmer_size)
    //create a map of all pairwise distances
    val distance_map = pairwiseDist(kmer_size)(sketches, List()).toMap
    println(timeStamp + "Writing to disk")
    //create distance matrix
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".matrix")
    //output header
    pw.println("$\t" + names.mkString("\t"))
    //iterate through each pairwise interaction in defined order and create rows
    names.foreach(subj => {
      //output row name
      pw.print(subj)
      names.foreach(target => {
        //get mash distance
        val md = if(subj == target) 0.0 else distance_map((subj, target))
        pw.print("\t" + md)
      })
      pw.println
    })
    pw.close()
    println(timeStamp + "Successfully completed!")

  }

  /**
    * Tail-recursive method to compute all pairwise distances in a given list of sketch files
    *
    * @param remaining Remaining list of sketch files
    * @param distances Accumulating pairwise distances as 3-tuples
    * @return List of 3-tuples: (subj, target, jaccard distance)
    */
  @tailrec def pairwiseDist(kmersize: Int)(remaining: List[File],
                            distances: List[((String, String), Double)]
                           ): List[((String, String), Double)] = {
    remaining match {
      case Nil => distances
      case (head :: tail) => {
        //load subj sketch
        val subj = loadRedwoodSketch(head)
        //iterate through remaining, and compute jaccard index
        val updated_distances = tail.foldLeft(distances)((acc, _target) => {
          //load target sketch
          val target = loadRedwoodSketch(_target)
          //get pairwise interaction
          val interaction = (subj.name, target.name)
          //compute union
          val union = broderUnion(subj.kmers.keySet, target.kmers.keySet)
          //compute intersection
          val intersect = union.intersect(subj.kmers.keySet).intersect(target.kmers.keySet).size.toDouble
          //compute jaccard distance
          val ji = (intersect / union.size)
          //compute mash distance
          val md = (-1.0/kmersize)*log((2*ji)/(1+ji))
          //return 3-tuple with jaccard distance
          (interaction.swap, md) :: ((interaction, md) :: acc)
        })
        pairwiseDist(kmersize)(tail, updated_distances)
      }
    }
  }

  /**
    * Method to obtain's Broder's union of two sets. In short, the smallest X hashes from the union of two sets.
    * @param x
    * @param y
    * @return
    */
  def broderUnion(x: Set[Int], y: Set[Int]): Set[Int] = {
    //get min sketch size to use
    val min_sketch_size = {
      //get size of each sketch
      val (x_size, y_size) = (x.size, y.size)
      //same sketch size
      if(x_size == y_size) x_size
      else {
        println("--Warning: uneven sketch sizes " + (x_size, y_size) + ". Using sketch size of " + min(x_size, y_size))
        //get min sketch size
        min[Int](x_size, y_size)
      }
    }
    //return Broder's union
    x.union(y).toList.sorted.take(min_sketch_size).toSet
  }

}