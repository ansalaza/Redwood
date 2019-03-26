package compare

import java.io.{File, PrintWriter}

import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyFile, verifyDirectory}
import utilities.SketchUtils.{loadSketches, loadRedwoodSketch}
import utilities.DistanceUtils._
import utilities.NumericalUtils.{choose, multiplyList}
import atk.ProgressBar.progress
import scala.annotation.tailrec

/**
  * Author: Alex N. Salazar
  * Created on 10-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object SketchDistance {

  case class Config(
                     sketchesFile: File = null,
                     outputDir: File = null,
                     prefix: String = null,
                     fullJI: Boolean = false,
                     log: Int = 1000
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("distance") {
      opt[File]("sketches") required() action { (x, c) =>
        c.copy(sketchesFile = x)
      } text ("File with the paths of all sketches to compare, one per line.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      opt[String]("prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output matrix file.")
      note("\nOPTIONAL\n")
      opt[Unit]("full-ji") action { (x, c) =>
        c.copy(fullJI = true)
      } text ("Do not express Jaccard Index as average genome size.")
      opt[Int]("log") hidden() action { (x, c) =>
        c.copy(log = x)
      } text ("Log process value.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists
      verifyDirectory(config.outputDir)
      verifyFile(config.sketchesFile)
      sketchDistance(config)
    }
  }

  def sketchDistance(config: Config): Unit = {
    println(timeStamp + "Processing sketches")
    //open sketches file and create map of id -> sketch path, kmer length, and sketch size
    val (sketches, klength, sketch_size) = {
      //load all paths as file objects
      val tmp = openFileWithIterator(config.sketchesFile).toList.map(new File(_))
      //verify file
      tmp.foreach(verifyFile(_))
      //load sketches
      loadSketches(tmp)
    }
    //set sketch names with unique int id
    val name2IDs = sketches.keySet.zipWithIndex.toMap
    //reverse above
    val id2Names = name2IDs.map(_.swap)
    println(timeStamp + "Found " + sketches.size + " sketches with kmer-length of " + klength + " and sketch " +
      "size of " + sketch_size)
    //set temp file
    val temp_file = new File(config.outputDir + "/.distances.tmp")
    val pw_temp = new PrintWriter(temp_file)

    @tailrec def computePairwiseDist(remaining: List[Int]): String = {
      remaining match {
        case Nil => "done"
        case (subj_id :: targets) => {
          //get subj name
          val subj_name = id2Names(subj_id)
          //load subj ID and kmers
          val subj_kmers = loadRedwoodSketch(sketches(subj_name)).sketch.keySet
          targets.par.foreach(target_id => {
            progress(config.log)
            //set target name
            val target_name = id2Names(target_id)
            //load target kmers
            val target_kmers = loadRedwoodSketch(sketches(target_name)).sketch.keySet
            //map as tuple of ids with mash distance
            pw_temp.println(subj_id + "\t" + target_id + "\t" + computeMashDist(subj_kmers, target_kmers, klength))
          })
          //compute mash distances in parallel (uses all threads) and move on with next
          computePairwiseDist(targets)
        }
      }
    }
    val (num, denom) = choose(sketches.size, 2)
    println(timeStamp + "Computing " + (multiplyList(num) / multiplyList(denom)) + " mash distances")
    //obtain map of all distances above diagonal
    val distance_map = computePairwiseDist(id2Names.keySet.toList)
    //close temp file
    pw_temp.close
    println(timeStamp + "Constructing distances map ")
    //set curried function to obtain distance of given two IDs
    val getDist = {
      val tmp = openFileWithIterator(temp_file).foldLeft(List[((Int,Int), Double)]())((acc,line) => {
        val columns = line.split("\t")
        ((columns(0).toInt, columns(1).toInt), columns(2).toDouble) :: acc
      }).toMap
      temp_file.delete()
      fetchDistance(tmp) _
    }
    println(timeStamp + "Writing to disk")
    //create distance matrix
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".matrix")
    //set name order
    val names = name2IDs.keySet
    //output header
    pw.println("$\t" + names.mkString("\t"))
    //iterate through each pairwise interaction in defined order and create rows
    names.foreach(subj => {
      //output row name
      pw.print(subj)
      names.foreach(target => {
        //get mash distance
        val md = if (subj == target) 0.0 else getDist(name2IDs(subj), name2IDs(target))
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
