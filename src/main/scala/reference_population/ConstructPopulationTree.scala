/**
package reference_population
import java.io._

import boopickle.Default._
import utilities.FileHandling.{timeStamp, verifyDirectory, verifyFile, writeSerialized, openFileWithIterator}
import utilities.ClusteringUtils.{createDendogram, hierchicalClustering}
import utilities.SketchUtils.{isKmerLengthCompatible, isUniqueNameCompatible, constructSketchesMap}
import utilities.DistanceUtils.loadMatrixColumns

/**
  * Author: Alex N. Salazar
  * Created on 4-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object ConstructPopulationTree {

  case class Config(
                     matrix: File = null,
                     sketchesFile: File = null,
                     prefix: String = null,
                     verbose: Boolean = false,
                     outputDir: File = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("population-tree") {
      opt[File]('m', "matrix-file") required() action { (x, c) =>
        c.copy(matrix = x)
      } text ("Tab-delimited distance matrix file.")
      opt[String]("prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output file.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      note("\nOPTIONAL\n")
      opt[File]("sketches") action { (x, c) =>
        c.copy(sketchesFile = x)
      } text ("File containing path to all sketches, one per line.")
      opt[Unit]("verbose") action { (x, c) =>
        c.copy(verbose = true)
      } text ("Verbose.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.matrix)
      if (config.sketchesFile != null) verifyFile(config.sketchesFile)
      constructPopulationTree(config)
    }
  }

  def constructPopulationTree(config: Config): Unit = {
    //perfom hierchical clustering of given distance matrix
    val reduced_kmer_tree = {
      println(timeStamp + "Loading distance matrix")
      //load matrix
      val columns = loadMatrixColumns(config.matrix)
      println(timeStamp + "--Loaded matrix with " + columns.size + " samples")
      //load sketches, if any
      val (all_sketches, universal_kmer_length) = {
        if(config.sketchesFile == null) (List[File](), -1)
        else {
          println(timeStamp + "Verifying sketches")
          val tmp = openFileWithIterator(config.sketchesFile).toList.map(new File(_))
          //first verify files actually exist
          tmp.foreach(verifyFile(_))
          println(timeStamp + "--Checking name uniqueness")
          //then, verify unique name
          val names = isUniqueNameCompatible(tmp)
          //assert all leaf and sketch names are accounted for
          assert(names.size == columns.size, "Not all samples/sketches are accounted for: " +
            names.size + " sketches vs " + columns.size + " columns")
          //then, verify kmer length compatability
          val kmer_length = isKmerLengthCompatible(tmp)
          println(timeStamp + "--Using universal kmer-length of " + kmer_length)
          (tmp, kmer_length)
        }
      }
      println(timeStamp + "Initializing dendogram")
      //create initial dendogram
      val dendo = createDendogram(config.matrix, universal_kmer_length, constructSketchesMap(config.sketchesFile), config.verbose)
      println(timeStamp + "Clustering samples and constructing kmer-tree")
      hierchicalClustering(dendo, 0, config.verbose)
    }
    println(timeStamp + "Writing to disk")
    writeSerialized(Pickle.intoBytes(reduced_kmer_tree).array(),
      new File(config.outputDir + "/" + config.prefix + ".rdwt"))
    println(timeStamp + "Successfully completed!")
  }
}
  */
