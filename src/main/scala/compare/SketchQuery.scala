package compare

import java.io.{File, PrintWriter}

import utilities.FileHandling.{timeStamp, verifyDirectory, verifyFile}
import utilities.SketchUtils.loadRedwoodSketch
import utilities.KmerTreeUtils.loadKmerTree

/**
  * Author: Alex N. Salazar
  * Created on 7-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object SketchQuery {

  case class Config(
                     sketchFile: File = null,
                     kmerTree: File = null,
                     prefix: String = null,
                     outputDir: File = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("classify") {
      opt[File]('s', "sketch") required() action { (x, c) =>
        c.copy(sketchFile = x)
      } text ("Sketch file to classify in redwood format.")
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
      verifyFile(config.sketchFile)
      sketchClassifier(config)
    }
  }

  def sketchClassifier(config: Config): Unit = {
    //load sketch
    val sketch = loadRedwoodSketch(config.sketchFile)
    //compute total kmers
     val total_kmers = sketch.sketch.size
    println(timeStamp + "Loaded sketch for " + sketch.name + " with " + total_kmers + " kmers of size " + sketch.kmer_length)
    println(timeStamp + "Loading kmer-tree")
    //load tree
    val ktree = loadKmerTree(config.kmerTree)
    println(timeStamp + "--" + ktree.getLeafNames().size + " leafs and " + ktree.getKmerSetSizes().size + " clusters ")
    println(timeStamp + "Querying sketch for " + sketch.name)
    //compute weights for each strain using LCA of each kmer
    val scores = sketch.sketch.keySet.toList.foldLeft(Map[String, Int]())((counts,kmer) => {
      //get lowest common ancestor
      val node = ktree.queryLCA(kmer)
      //update counts
      if(node.isEmpty) counts else counts + (node.get -> (counts.getOrElse(node.get, 0)+1))
    }).toList.map(x => (x._1, x._2.toDouble / total_kmers)).sortBy(-_._2)
    //explained percentage
    val explained_percentage = scores.foldLeft(0.0)((b,a) => a._2 + b) / total_kmers
    println(timeStamp + "Fraction of kmers explained: " + explained_percentage)
    //create output file
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".txt")
    pw.println("Strain\tWeight")
    scores.foreach{case(strain,weight) => pw.println(strain + "\t" + weight)}
    pw.close()
    println(timeStamp + "Successfully completed!")
  }

}
