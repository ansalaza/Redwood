package reference_population

import java.io._
import boopickle.Default._
import utilities.FileHandling.{timeStamp, verifyDirectory, verifyFile}
import utilities.KmerTreeClusteringUtils.{createDistanceMatrix, hierchicalClustering}

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
      if(config.sketchesFile != null) verifyFile(config.sketchesFile)
      constructPopulationTree(config)
    }
  }

  def constructPopulationTree(config: Config): Unit = {
    //perfom hierchical clustering of given distance matrix
    val dendogram = {
      println(timeStamp + "Loading distance matrix")
      //load distance matrix
      val tmp = createDistanceMatrix(config.matrix, config.sketchesFile, config.verbose)
      println(timeStamp + "--Loaded matrix with " + tmp.labels.size + " samples and " + tmp.pairwise_leaf_dist.size +
        " values")
      println(timeStamp + "Clustering samples and constructing kmer-tree")
      hierchicalClustering(tmp, tmp.sketch_map, config.verbose)
    }
    println(timeStamp + "Writing to disk")
    //serialize to byte array
    val buf = Pickle.intoBytes(dendogram).array()
    //create output file
    val pw = new BufferedOutputStream(new FileOutputStream(config.outputDir + "/" + config.prefix + ".rdwt"))
    //write bytes to output
    pw.write(buf)
    pw.close
    println(timeStamp + "Successfully completed!")

    //println(timeStamp + "--Deserializing")
    //println(Unpickle[Node].fromBytes(
     // ByteBuffer.wrap(Files.readAllBytes(Paths.get(config.outputDir + "/" + config.prefix + ".rdw")))
    //))
  }
}
