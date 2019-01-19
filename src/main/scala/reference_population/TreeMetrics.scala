package reference_population

import java.io.{File, PrintWriter}

import utilities.FileHandling.{timeStamp, verifyDirectory, verifyFile}
import utilities.KmerTreeUtils.loadKmerTree

/**
  * Author: Alex N. Salazar
  * Created on 7-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object TreeMetrics {
  case class Config(
                     kmerTree: File = null,
                     prefix: String = null,
                     outputDir: File = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("tree-metrics") {
      opt[File]('t', "kmer-tree") required() action { (x, c) =>
        c.copy(kmerTree = x)
      } text ("Kmer-tree as constructed by redwood.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      opt[String]("prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output file.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.kmerTree)
      treeMetrics(config)
    }
  }

  def treeMetrics(config: Config): Unit = {
    //load kmer tree
    println(timeStamp + "Loading kmer-tree")
    val ktree = loadKmerTree(config.kmerTree)
    //get all leafs
    val leafs = ktree.getLeafNames()
    //get all kmer set sizes
    val set_sizes = ktree.getKmerSetSizes()
    println(timeStamp + "Found " + leafs.size + " leafs with " + set_sizes.size + " clusters")
    println(timeStamp + "Writing to disk")
    //create output file
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".txt")
    //output header
    pw.println("Kmers\tLeafs\tIdentifier")
    //output metrics
    set_sizes.foreach(cluster => pw.println(cluster._2 + "\t" + (cluster._1.count(_ == ',') + 1) + "\t" + cluster._1))
    pw.close
    println(timeStamp + "Successfully completed!")
  }


}
