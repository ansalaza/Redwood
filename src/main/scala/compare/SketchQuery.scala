package compare

import java.io.{File, PrintWriter}

import utilities.FileHandling.{timeStamp, verifyDirectory, verifyFile, openFileWithIterator}
import utilities.SketchUtils.loadRedwoodSketch
import utilities.KmerTreeUtils._

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
                     outputDir: File = null,
                     labelsFile: File = null,
                     sketchesFile: File = null
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("query") {
      opt[File]('s', "sketch") action { (x, c) =>
        c.copy(sketchFile = x)
      } text ("Sketch file to query.")
      opt[File]('t', "kmer-tree") required() action { (x, c) =>
        c.copy(kmerTree = x)
      } text ("Kmer-tree as constructed by redwood.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      opt[String]("prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output file.")
      note("\nOPTIONAL\t")
      opt[File]("sketches") action { (x, c) =>
        c.copy(sketchesFile = x)
      } text ("Alternatively, provide a file with paths to multiple sketches, one per line. Output will contain " +
        "summary information jointly for all sketches.")
      opt[File]("labels") action { (x, c) =>
        c.copy(labelsFile = x)
      } text ("Tab-delimited file containing: sample ID, label. Output will contain summary information per label")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.kmerTree)
      verifyFile(config.sketchFile)
      assert(config.sketchesFile != null || config.sketchFile != null, "Provide sketch file or file containing path " +
        "to multiple sketches")
      assert(!(config.sketchesFile != null && config.sketchFile != null), "Provide sketch file or file containing " +
        "path to multiple sketches")
      querySketch(config)
    }
  }

  def querySketch(config: Config, loaded_tree: Tree[Kmers] = null): Unit = {
    //load sketch
    val sketches = {
      //single sketch file provided
      if (config.sketchFile != null) List(config.sketchFile)
      //multiple sketch files provided
      else openFileWithIterator(config.sketchesFile).toList.map(x => {
        val file = new File(x); verifyFile(file); file
      })
    }
    println(timeStamp + "Found " + sketches.size + " sketches to query")
    println(timeStamp + "Loading kmer-tree")
    //load tree
    val ktree = if(loaded_tree != null) loaded_tree else loadKmerTree(config.kmerTree)
    if(loaded_tree == null) {
      println(timeStamp + "--" + ktree.getLeafNames().size + " leafs and " + ktree.getKmerSetSizes().size + " clusters ")
    }
    println(timeStamp + "Querying sketches:")
    //compute weights for each strain using LCA of each kmer
    val (scores, total_kmers) = sketches.foldLeft((Map[String, Int](),0)){ case ((counts, totalk), _sketch) => {
      //load sketch
      val sketch = loadRedwoodSketch(_sketch)
      println(timeStamp + "--" + sketch.name)
      //iterate through each kmer and query
      sketch.sketch.keySet.foldLeft((counts, totalk)){ case ((acc_counts, acc_total), kmer) => {
        //get lowest common ancestor
        val node = ktree.queryLCA(kmer)
        //update counts
        if (node.isEmpty) (acc_counts, acc_total + 1)
        else (acc_counts + (node.get -> (acc_counts.getOrElse(node.get, 0) + 1)), acc_total + 1)
      }}
    }}
    //compute branch proportions based on scores above
    val branch_proportions = calcualteBranchWeights(ktree, scores).mapValues(_.toDouble / total_kmers)
    //explained percentage
    val explained_percentage = scores.foldLeft(0.0)((b, a) => a._2 + b) / total_kmers
    println(timeStamp + "Fraction of kmers explained: " + explained_percentage)
    //create output file
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".txt")
    pw.println("Strain\tWeight")
    branch_proportions.toList.sortBy(-_._2) foreach { case (strain, weight) => pw.println(strain + "\t" + weight) }
    pw.close()

    //labels provided
    if (config.labelsFile != null) {
      //load labels as map leaf -> label
      val leaf2Labels = {
        //open labels file
        val tmp = openFileWithIterator(config.labelsFile).toList.map(x => {
          val columns = x.split("\t")
          (columns.head, columns(1))
        }).toMap
        //sanity check
        assert(tmp.size == ktree.getLeafNames().size, "Unexpected number of leafs found in labels file")
        tmp
      }
      //construct map as node -> labels -> weight
      val node2Labels = constructNode2Labels(config.labelsFile, ktree)
      //for each node, split it counts to it's corresponding labels proportionally
      val label2Scores = scores.foldLeft(Map[String, Double]()) { case (props, (node, count)) => {
        //get labels for current node
        val local_labels = node2Labels(node)
        //update props via the sum of current prop and the product of local label weight and current node prop
        local_labels.foldLeft(props) { case (acc, (label, weight)) =>
          acc + (label -> (acc.getOrElse(label, 0.0) + weight * count))
        }
      } //proportion overall
      }.mapValues(_ / total_kmers)
      //create output file
      val pw2 = new PrintWriter(config.outputDir + "/" + config.prefix + ".labels.txt")
      pw2.println("Label\tWeight")
      label2Scores.foreach(x => pw2.println(x._1 + "\t" + x._2))
      pw2.close
    }

    println(timeStamp + "Successfully completed!")
  }

  def calcualteBranchWeights(t: Tree[Kmers], scores: Map[String, Int]): Map[String, Int] = {
    def _calcualteBranchWeights(current: Tree[Kmers], acc: Map[String, Int]): Map[String, Int] = {
      current match {
        case Leaf(a, b) => acc + (b -> (scores.getOrElse(b, 0)))
        case Node(i, k, d, l, r) => {
          //get cumulative weights of children
          val updated_prop = _calcualteBranchWeights(l, acc) ++ _calcualteBranchWeights(r, acc)
          //calculate current node's weight
          val node_weight = scores.getOrElse(i.toString, 0) + updated_prop(l.getID()) + updated_prop(r.getID())
          //update acc
          updated_prop + (i.toString -> node_weight)
        }
      }
    }

    _calcualteBranchWeights(t, Map())
  }

  def constructNode2Labels(file: File, tree: Tree[Kmers]): Map[String, Map[String, Double]] = {
    //construct map as node -> labels -> weight
    val node2labels = {
      //open labels file
      val tmp = openFileWithIterator(file).toList.map(x => {
        val c = x.split("\t");
        (c.head, c(1))
      }).toMap
      //sanity check
      assert(tmp.size == tree.getLeafNames().size, "Unexpected number of leafs found in labels file")
      //construct node -> labels
      tree.id2Leafnames().toMap.mapValues(_.map(tmp(_)).toSet)
    }
    //assign weights to the labels of each node
    node2labels.mapValues(x => {
      val max = x.size.toDouble;
      x.groupBy(identity).mapValues(_.size / max)
    })
  }

}
