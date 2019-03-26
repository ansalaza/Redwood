package reference_population

import java.io.{File, PrintWriter}

import utilities.FileHandling.{timeStamp, verifyDirectory, verifyFile}
import utilities.ReducedKmerTreeUtils.loadReducedKmerTree

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
    val ktree = loadReducedKmerTree(config.kmerTree)
    //get leaf id -> name
    val leafID2Name = ktree.tree.getLeafId2Name()
    //get node -> leaf names counts
    val (node2Leafs, node2Chilren) = {
      val tmp = ktree.tree.node2Children()
      (tmp.mapValues(x => x.filter(leafID2Name.contains(_)).size), tmp.mapValues(_.size))
    }
    //get all nodes
    val nodes = ktree.tree.getNodeIDsPostOrder().size
    //get all leafs
    val leafs = ktree.tree.getLeafNamesPostOrder().size
    //get node -> kmers
    val node2kmers = ktree.node2totalkmers
    println(timeStamp + "Found " + leafs + " leafs with " + nodes + " nodes")
    println(timeStamp + "Writing to disk")
    //create output file
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".txt")
    //output header
    pw.println("ID\tName\tNodes\tLeafs\tKmers")
    //output metrics
    node2kmers.toList.sortBy(_._2).foreach(node => {
      val name = leafID2Name.get(node._1)
      pw.println(node._1 + "\t" + (if(name.isEmpty) "" else name.get) + "\t" + node2Chilren(node._1) + "\t" +
        node2Leafs(node._1) + "\t" + node._2.toInt)
    })
    pw.close
    println(timeStamp + "Successfully completed!")
  }


}
