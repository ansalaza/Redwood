package compare

import java.io.{File, PrintWriter}

import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory}
import utilities.SketchUtils.{loadSketches, loadRedwoodSketch}
import utilities.DistanceUtils.{computeMashDist, fetchDistance,loadMatrix}
import atk.ProgressBar.progress
import scala.annotation.tailrec
import scala.collection.parallel.ForkJoinTaskSupport

/**
  * Author: Alex N. Salazar
  * Created on 24-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object UpdateMatrix {

  case class Config(
                     matrix: File = null,
                     oldSketches: File = null,
                     newSketches: File = null,
                     outputDir: File = null,
                     prefix: String = null,
                     threads: Int = 1,
                     log: Int = 1000
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("update") {
      opt[File]('m', "matrix") required() action { (x, c) =>
        c.copy(matrix = x)
      } text ("Matrix to update.")
      opt[File]('e', "existing-sketches") required() action { (x, c) =>
        c.copy(oldSketches = x)
      } text ("File containing path to all sketches in given matrix, one per line.")
      opt[File]('n', "new-sketches") required() action { (x, c) =>
        c.copy(newSketches = x)
      } text ("File containing path to new sketches to add to matrix, one per line.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      opt[String]("prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for updated matrix.")
      note("\nOPTIONAL\n")
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
      //get old sketches
      val old = loadSketches(openFileWithIterator(config.oldSketches).toList.map(new File(_)))
      //get new sketches
      val news = loadSketches(openFileWithIterator(config.newSketches).toList.map(new File(_)))
      //sanity check
      assert(old._2 == news._2, "Kmer-length from existing sketches is not the same as kmer-length from new sketches")
      //get new ids already in the matrix, if any
      val repeat_ids = old._1.keySet.intersect(news._1.keySet)
      if (repeat_ids.nonEmpty) {
        println(timeStamp + repeat_ids.size + " sketches are already in the matrix: " + repeat_ids.mkString(","))
        println(timeStamp + "--These sketches will be removed")
      }
      assert(config.threads > 0, "Number of threads specified is a non-positive integer")
      updateMatrix(old._1.mapValues(loadRedwoodSketch(_).sketch.keySet),
        news._1.filterNot(x => repeat_ids(x._1)).mapValues(loadRedwoodSketch(_).sketch.keySet), old._2, config)
    }
  }

  def updateMatrix(old: Map[String, Set[Int]], news: Map[String, Set[Int]], kmer_length: Int, config: Config): Unit = {
    println(timeStamp + "Loaded " + old.size + " existing sketches and " + news.size + " new sketches")
    //load matrix
    val (existing_distances, all_ids) = loadMatrix(config.matrix)
    //sanity check
    assert(all_ids.size == old.size, "Unexpected number of samples in matrix (" + all_ids.size + ") and samples in" +
      " existing sketches (" + old.size + ")")
    println(timeStamp + "Loaded matrix with " + existing_distances.size + " distances (upper diagonal)")
    //update distance map
    val updated_distances = {
      println(timeStamp + "Generating pairwise interactions")
      //generate new pairwise instances
      val pairwise = generateNewPairwise(news.keys.toList, old.keys.toList).par
      //set max threads
      pairwise.tasksupport = new ForkJoinTaskSupport(new scala.concurrent.forkjoin.ForkJoinPool(config.threads))
      println(timeStamp + "Computing " + pairwise.size + " mash distances")
      //obtain mash distances in parallel
      pairwise.map{ case (s,t) => {
        progress(config.log)
        //compute mash distance
        val d = computeMashDist(if(old.contains(s)) old(s) else news(s), if(old.contains(t)) old(t) else news(t), kmer_length)
        //set as 2-tuple with distance
        ((s,t),d)
      }}.seq.toMap ++ existing_distances
    }
    val getDist = fetchDistance(updated_distances) _
    println(timeStamp + "Writing to disk")
    //create distance matrix
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".matrix")
    //set id order
    val ids = old.keySet.toList ::: news.keySet.toList
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
    * @param newids List of new IDs
    * @param oldids List of old IDs
    * @return List[(String,String)
    */
  def generateNewPairwise(newids: List[String],
                          oldids: List[String]
                         ): List[(String, String)] = {

    /**
      *
      * @param remaining    Remaining IDs from new IDs
      * @param interactions (Accumulating) pairwise interactions
      * @return List[(String, String)]
      */
    @tailrec def _generateNewPairwise(remaining: List[String],
                                      interactions: List[(String, String)]
                                     ): List[(String, String)] = {
      remaining match {
        //no more IDs to process
        case Nil => interactions
        //generate pairwise with current ID
        case (s :: tail) =>
          //first generate pairwise using new ids, then generate using old ids
          _generateNewPairwise(tail,
            oldids.foldLeft(tail.foldLeft(interactions)((nacc, t) => (s, t) :: nacc))((oacc, t) => (s, t) :: oacc)
          )
      }
    }

    _generateNewPairwise(newids, List())
  }

}
