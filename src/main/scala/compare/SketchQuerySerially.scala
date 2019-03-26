package compare

import java.io.File

import utilities.FileHandling.{verifyDirectory, verifyFile, timeStamp, openFileWithIterator}
import utilities.SketchUtils.loadRedwoodSketch
import utilities.ReducedKmerTreeUtils.loadReducedKmerTree
/**
  * Author: Alex N. Salazar
  * Created on 27-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object SketchQuerySerially {

  case class Config(
                     kmerTree: File = null,
                     prefix: String = null,
                     outputDir: File = null,
                     labelsFile: File = null,
                     sketchesFile: File = null,
                     normalize: Boolean = false,
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("query-serially") {
      opt[File]("sketches") required() action { (x, c) =>
        c.copy(sketchesFile = x)
      } text ("File with paths to multiple sketches, one per line.")
      opt[File]('t', "kmer-tree") required() action { (x, c) =>
        c.copy(kmerTree = x)
      } text ("Kmer-tree as constructed by redwood.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      note("\nOPTIONAL\t")
      opt[File]("labels") action { (x, c) =>
        c.copy(labelsFile = x)
      } text ("Tab-delimited file containing: sample ID, label. Output will contain summary information per label")
      opt[Unit]("normalize") action { (x,c) =>
        c.copy(normalize = true)
      } text ("Normalize queries to non-unique proportion.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.kmerTree)
      verifyFile(config.sketchesFile)
      querySerially(config)
    }
  }

  def querySerially(config: Config): Unit = {
    //load sketch files
    val sketches = openFileWithIterator(config.sketchesFile).toList.map(x => {
      val file = new File(x)
      //sanity check
      verifyFile(file)
      file
    })
    println(timeStamp + "Found " + sketches.size + " sketches to query")
    println(timeStamp + "Pre-loading tree")
    val ktree = loadReducedKmerTree(config.kmerTree)
    println(timeStamp + "--" + ktree.tree.getLeafId2Name().size + " leafs and " +  " " + "clusters ")
    //iterate through each sketch and query
    sketches.foreach(sketch => {
      //get sketch name
      val name = loadRedwoodSketch(sketch).name
      //query current sketch
      SketchQuery.querySketch(
        new SketchQuery.Config(sketch, config.kmerTree, name, config.outputDir, config.labelsFile, null, config.normalize),
        ktree)
    })
  }


}
