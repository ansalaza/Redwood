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
object SketchClassifier {

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
     val total_kmers = sketch.kmers.size
    println(timeStamp + "Loaded sketch for " + sketch.name + " with " + total_kmers + " kmers of size " + sketch.kmerSize)
    println(timeStamp + "Loading kmer-tree")
    //load tree
    val ktree = loadKmerTree(config.kmerTree)
    println(timeStamp + "--" + ktree.getLeafNames().size + " leafs and " + ktree.getKmerSetSizes().size + " clusters ")
    println(timeStamp + "Querying sketch for " + sketch.name)
    //compute weights for each strain using LCA of each kmer
    val weights = sketch.kmers.keySet.toList.foldLeft(Map[String, Double]())((strain_weights,kmer) => {
      //get lowest common ancestor
      val leafs = ktree.queryLCA(kmer)
      //assign weight to kmer
      val weight = if(leafs.isEmpty) 0.0 else (1 / leafs.size.toDouble)
      //update weights
      leafs.foldLeft(strain_weights)((acc, strain) => acc + (strain -> (acc.getOrElse(strain, 0.0) + weight)))
    }).toList.sortBy(-_._2)
    //explained percentage
    val explained_percentage = weights.foldLeft(0.0)((b,a) => a._2 + b) / total_kmers
    println(timeStamp + "Fraction of kmers explained: " + explained_percentage)
    //create output file
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".txt")
    pw.println("Strain\tWeight")
    weights.foreach{case(strain,weight) => pw.println(strain + "\t" + weight)}
    pw.close()
    println(timeStamp + "Successfully completed!")
  }

}
