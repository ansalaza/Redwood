package compare

import java.io.{File, PrintWriter}

import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory}
import utilities.SketchUtils.{loadSketches, loadRedwoodSketch}
import utilities.DistanceUtils._
import utilities.NumericalUtils.choose
import atk.ProgressBar.progress
import scala.annotation.tailrec
import scala.collection.parallel.ForkJoinTaskSupport

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
                     log: Int = 1000
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("distance") {
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
      } text ("Max number of threads.")
      opt[Int]("log") hidden() action { (x, c) =>
        c.copy(log = x)
      } text ("Log process value.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists
      verifyDirectory(config.outputDir)
      //get list of sketches from input parameter
      println(timeStamp + "Loading sketches")
      val (all_sketches, kmer_length, sketch_size) = {
        if (config.sketchFiles == null && config.pathsFile == null) {
          assert(false, "Provide sketch files through '--sketch-files' or '--paths-file'")
          (Map[String, File](), 0, 0)
        }
        else if (config.sketchFiles != null) loadSketches(config.sketchFiles.toList)
        else loadSketches(openFileWithIterator(config.pathsFile).toList.map(x => new File(x)))
      }
      assert(config.threads > 0, "Number of threads specified is a non-positive integer")
      sketchDistance(all_sketches.mapValues(loadRedwoodSketch(_).sketch.keySet), kmer_length, config)
    }
  }

  def sketchDistance(sketches: Map[String, Set[Int]], kmer_length: Int, config: Config): Unit = {
    println(timeStamp + "Loaded " + sketches.size + " sketches with kmer-length of " + kmer_length)
    //set id order
    val ids = sketches.keys.toList
    //obtain map of all distances above diagonal
    val distance_map = {
      println(timeStamp + "Generating pairwise interactions")
      //get all pairwise interactions above diagonal
      val pairwise = generatePairwise(ids, List()).par
      //set max threads
      pairwise.tasksupport = new ForkJoinTaskSupport(new scala.concurrent.forkjoin.ForkJoinPool(config.threads))
      println(timeStamp + "Computing " + pairwise.size +  " mash distances")
      //obtain mash distances in parallel
      pairwise.map { case (s, t) => {
        progress(config.log)
        ((s, t), computeMashDist(sketches(s), sketches(t), kmer_length)
        )}}.seq.toMap
    }
    /**
      * Set curryied function to obtain distance of given two IDs
      */
    val getDist = fetchDistance(distance_map) _
    println(timeStamp + "Writing to disk")
    //create distance matrix
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".matrix")
    //output header
    pw.println("$\t" + ids.mkString("\t"))
    //iterate through each pairwise interaction in defined order and create rows
    ids.foreach(subj => {
      //output row name
      pw.print(subj)
      ids.foreach(target => {
        //get mash distance
        val md = if (subj == target) 0.0 else getDist(subj, target)
        pw.print("\t" + md)
      })
      pw.println
    })
    pw.close()
    println(timeStamp + "Successfully completed!")

  }

  /**
    * Method to generate all pairwise interactions from a given list of IDs
    *
    * @param ids          List of IDs
    * @param interactions (Accumulating) pairwise interactions
    * @return List[(String,String)
    */
  @tailrec def generatePairwise(ids: List[String],
                                interactions: List[(String, String)],
                               ): List[(String, String)] = {
    ids match {
      //no more IDs to process
      case Nil => interactions
      //generate pairwise with current ID
      case (subj :: tail) => generatePairwise(tail, tail.foldLeft(interactions)((acc, target) => (subj, target) :: acc))
    }
  }

}
