package compare

import java.io.{File, PrintWriter}

import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory, verifyFile}
import utilities.SketchUtils.loadRedwoodSketch
import utilities.ReducedKmerTreeUtils._
import utilities.NumericalUtils.choose

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
      } text ("Prefix for output file. Also used as query identifier.")
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

  def querySketch(config: Config, loaded_tree: ReducedKmerTree = null): Unit = {
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
    val ktree = if(loaded_tree != null) loaded_tree else loadReducedKmerTree(config.kmerTree)
    //log info
    if(loaded_tree == null) {
      println(timeStamp + "--" + ktree.tree.getLeafId2Name().size + " leafs and " +  " clusters ")
      println(timeStamp + "--Global sketch size of tree is " + ktree.min_sketch_size)
    }
    //set curried-probability function
    //val getProb = computeProbability(ktree.kmer_length, ktree.min_sketch_size)
    println(timeStamp + "Querying sketches:")
    //compute weights for each strain using LCA of each kmer
    val (scores, total_kmers) = sketches.foldLeft((Map[Int, Double](),0)){ case ((counts, totalk), _sketch) => {
      //get kmers and sketch name
      val (kmers, name) = {
        //load sketch
        val sketch = loadRedwoodSketch(_sketch)
        println(timeStamp + "--" + sketch.name)
        //sketch size is the same as global sketch size in tree, return as is
        if(sketch.sketch.size == ktree.min_sketch_size) (sketch.sketch.keySet, sketch.name)
        //need to adjust kmers in sketch
        else {
          println("----Adjusting sketch size to global sketch size in tree")
          //adjust by returning the smallest X kmers in the sketch (X = global sketch size)
          (sketch.sketch.keySet.toList.sorted.take(ktree.min_sketch_size), sketch.name)
        }
      }
      //iterate through each kmer and query
        kmers.foldLeft((counts, totalk)){ case ((acc_counts, acc_total), kmer) => {
        //get lowest common ancestor
        val node = ktree.lca_map.get(kmer)
        //update counts
        if (node.isEmpty) (acc_counts, acc_total + 1)
        else (acc_counts + (node.get -> (acc_counts.getOrElse(node.get, 0.0) + 1)), acc_total + 1)
      }}
    }}
    //compute branch proportions based on scores above
    val branch_proportions = ktree.genericCumulative(scores).mapValues(_ / total_kmers)
    //explained percentage
    val explained_percentage = scores.foldLeft(0.0)((b, a) => a._2 + b) / total_kmers
    //get leaf ID -> name
    val leafid2name = ktree.tree.getLeafId2Name()
    println(timeStamp + "Fraction of kmers explained: " + explained_percentage)
    //create output file
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".txt")
    pw.println("ID\tName\tWeight")
    branch_proportions.toList.sortBy(-_._2) foreach { case (id, weight) => {
      val name = leafid2name.get(id)
      pw.println(id + "\t" + (if(name.isEmpty) "" else name.get) + "\t" + weight)
    } }
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
        assert(tmp.size == ktree.tree.getLeafId2Name().size, "Unexpected number of leafs found in labels file")
        tmp
      }
      //construct map as node -> labels -> weight
      val node2Labels = constructNode2Labels(config.labelsFile, ktree.tree)
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
      pw2.println("Name\tLabel\tWeight")
      label2Scores.foreach(x => pw2.println(config.prefix + "\t" + x._1 + "\t" + x._2))
      pw2.println(config.prefix + "\tUnique\t" + (1.0 - explained_percentage))
      pw2.close
    }

    println(timeStamp + "Successfully completed!")
  }

  /**
    * Curried-method to compute the probability of sharing at least X kmers in a given node
    * @param kmer_length
    * @param sketch_size
    * @param size1
    * @param size2
    * @param shared
    * @return
    */
  def computeProbability(kmer_length: Int, sketch_size: Int)
                        (size1: Int, size2: Int, shared: Int): Double =  {
    /**
      * Function to compute probability of a random kmer in a random genome of given length
      * @return Double
      */
    def probRandomKmer: Int => Double = x => (x.toDouble / (x + Math.pow(4, kmer_length)))

    /**
      * Function to compute expected jaccard index between to random genomes given their corresponding probability of
      * observing a random kmer
      * @return Double
      */
    def expectedJI: (Double, Double) => Double =(x,y) => (x*y)/(x+y - (x*y))

    //set expected jaccard index
    val eji = expectedJI(probRandomKmer(size1), probRandomKmer(size2))

    //compute probability
    //(0 to shared).view.foldLeft(0.0)((acc, i) => {acc + choose(sketch_size, i)*Math.pow(eji, i) * Math.pow(1 - eji, sketch_size - i)})
    ???
  }

  /**
    *
    * @param file
    * @param tree
    * @return
    */
  def constructNode2Labels(file: File, tree: ReducedTree): Map[Int, Map[String, Double]] = {
    //construct map as node -> List(labels)
    val node2labels = {
      //open labels file
      val leafname2label = openFileWithIterator(file).toList.map(x => {
        val c = x.split("\t");
        (c.head, c(1))
      }).toMap
      //sanity check
      assert(leafname2label.size == tree.getLeafId2Name().size, "Unexpected number of leafs found in labels file")
      //construct node -> labels
      tree.getId2LeafNames().mapValues(_.map(leafname2label(_)))
    }
    //assign weights to the labels of each node
    node2labels.mapValues(x => {
      val max = x.size.toDouble;
      x.groupBy(identity).mapValues(_.size / max)
    })
  }

}
