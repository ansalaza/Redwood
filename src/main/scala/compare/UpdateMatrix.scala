package compare

import java.io.{File, PrintWriter}

import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory, verifyFile}
import utilities.DistanceUtils.{computeDistance, fetchDistance, loadMatrix}
import atk.ProgressBar.progress
import utilities.NumericalUtils.{choose, multiplyList}
import utilities.SketchUtils.{loadRedwoodSketch, loadSketches}

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
      opt[Int]("log") hidden() action { (x, c) =>
        c.copy(log = x)
      } text ("Log process value.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists
      verifyDirectory(config.outputDir)
      verifyFile(config.matrix)
      verifyFile(config.newSketches)
      verifyFile(config.oldSketches)
      updateMatrix(config)
    }
  }

  def updateMatrix(config: Config): Unit = {
    //load columns from matrix
    val (old_distances, old_names) = loadMatrix(config.matrix)
    println(timeStamp + "Processing existing sketches")
    //open new sketches file
    val (joined_sketches, new_names, klength) = {
      //open sketches file and create map of id -> sketch path, kmer length, and sketch size
      val (old_sketches, _klength, sketch_size) = {
        //load all paths as file objects
        val tmp = openFileWithIterator(config.oldSketches).toList.map(new File(_))
        //verify file
        tmp.foreach(verifyFile(_))
        //load sketches
        loadSketches(tmp)
      }
      println(timeStamp + "Processed existing sketches with kmer length of " + _klength + " and sketch size of " + sketch_size)
      println(timeStamp + "Processing new sketches")
      //load all paths as file objects
      val tmp = openFileWithIterator(config.newSketches).toList.map(new File(_))
      //verify file
      tmp.foreach(verifyFile(_))
      //load sketches
      val t = loadSketches(tmp)
      //sanity checks
      assert(_klength == t._2, "Kmer length from existing sketches (" + _klength + ") does not match kmer length of new" +
        " sketches")
      assert(sketch_size == t._3, "Sketch size of existing sketches (" + sketch_size + ") does not match sketch size " +
        "of new sketches")
      t._1.keys.foreach(x => assert(!old_sketches.contains(x), "Sketch " + x + " already exists existing sketches"))
      //return map
      (old_sketches ++ t._1, t._1.keys.toList, _klength)
    }
    println(timeStamp + "Processed " + new_names.size + " new sketches to add")
    //set name -> assigned ID
    val name2IDs = new_names.foldLeft(old_names)((acc, a) => acc.:+(a)).zipWithIndex.toMap
    //reverse above
    val id2Names = name2IDs.map(_.swap)
    //set temp file
    val temp_file = new File(config.outputDir + "/.distances.tmp")
    val pw_temp = new PrintWriter(temp_file)
    //set old name set
    val oldNamesSet = old_names.toSet

    @tailrec def computePairwiseDist(remaining: List[Int]): String = {
      remaining match {
        case Nil => "done"
        case (subj_id :: targets) => {
          //get subj name
          val subj_name = id2Names(subj_id)
          //load subj ID and kmers
          lazy val subj_kmers = loadRedwoodSketch(joined_sketches(subj_name)).sketch.keySet
          targets.par.foreach(target_id => {
            progress(config.log)
            //set target name
            val target_name = id2Names(target_id)
            //not already existing distances
            if(!(oldNamesSet(subj_name) && oldNamesSet(target_name))) {
              //load target kmers
              val target_kmers = loadRedwoodSketch(joined_sketches(target_name)).sketch.keySet
              //map as tuple of ids with mash distance
              pw_temp.println(subj_id + "\t" + target_id + "\t" + computeDistance(subj_kmers, target_kmers, klength))
            }
          })
          //compute mash distances in parallel (uses all threads) and move on with next
          computePairwiseDist(targets)
        }
      }
    }
    val (oldnum, olddenum) = choose(old_names.size, 2)
    val (newnum, newdenom) = choose(new_names.size, 2)
    println(timeStamp + "Verifying " + (multiplyList(oldnum) / multiplyList(olddenum)) + " existing distances and " +
      "computing " + ((multiplyList(newnum) / multiplyList(newdenom)) + (new_names.size * old_names.size)) +
      " additional mash distances")
    //obtain map of all distances above diagonal
    val distance_map = computePairwiseDist(id2Names.keySet.toList)
    //close temp file
    pw_temp.close
    println(timeStamp + "Constructing distances map ")
    System.gc()
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
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".updated.matrix")
    //set name order
    val names = name2IDs.toList.sortBy(_._2).map(_._1)
    //output header
    pw.println("$\t" + names.mkString("\t"))
    //iterate through each pairwise interaction in defined order and create rows
    names.foreach(subj => {
      //output row name
      pw.print(subj)
      names.foreach(target => {
        //get mash distance
        val md = {
          if (subj == target) 0.0
          else {
            if(old_distances.contains(subj, target)) old_distances(subj,target)
            else if(old_distances.contains(target, subj)) old_distances(target,subj)
            else getDist(name2IDs(subj), name2IDs(target))
          }
        }
        pw.print("\t" + md)
      })
      pw.println
    })
    pw.close()
    println(timeStamp + "Successfully completed!")

  }

}