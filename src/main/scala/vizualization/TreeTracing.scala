package vizualization

/**
  * Author: Alex N. Salazar
  * Created on 20-2-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */

import java.io.File

import utilities.FileHandling.{timeStamp, verifyDirectory, verifyFile}
import utilities.KmerTreeUtils._
import utilities.TreeDrawingUtils._
import doodle.syntax._
import doodle.jvm.Java2DFrame._
import doodle.backend.StandardInterpreter._
import doodle.jvm.FileFrame.pdfSave
import doodle.core.font.Font
import doodle.core.font.FontFace.Bold
import doodle.core.font.FontFamily.Monospaced
import doodle.core.font.FontSize.Points


object TreeTracing {

  case class Config(
                     tree: File = null,
                     outputDir: File = null,
                     canvasHeight: Int = 1000,
                     canvasWidth: Int = 1000,
                     fontSize: Int = 14,
                     coloursFile: File = null,
                     proportions: File = null,
                     prefix: String = null,
                     lineWidth: Int = 1,
                     minFreq: Double = 0.005,
                     equalDist: Boolean = false,
                     features: File = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("tree-tracing") {
      opt[File]('t', "tree") required() action { (x, c) =>
        c.copy(tree = x)
      } text ("Kmer tree as constructed by Redwood.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      opt[String]("prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output file name.")
      note("\nOPTIONAL METADATA")
      opt[File]("colours") action { (x, c) =>
        c.copy(coloursFile = x)
      } text ("Tab-delimited file containing sample ID or label and colour as hue, saturation, and lightness values.")
      opt[File]("query") action { (x, c) =>
        c.copy(proportions = x)
      } text ("Use Redwood's 'query' output to reflect frequency of most-frequent traversed path. ")
      opt[Double]("ignore-frequency") action {(x,c) =>
        c.copy(minFreq = x)
      } text ("Minimum frequency to display from a query file (default is 0.005).")
      note("\nOPTIONAL AESTHETICS")
      opt[Unit]("equal-dist") action { (x, c) =>
        c.copy(equalDist = true)
      } text ("Force equal branch distances.")
      opt[Int]("font-size") action { (x, c) =>
        c.copy(fontSize = x)
      } text ("Font-size of labels (default is 14).")
      opt[Int]("line-width") action { (x, c) =>
        c.copy(lineWidth = x)
      } text ("Line width for tree (default is 1).")
      opt[Int]("image-height") action { (x, c) =>
        c.copy(canvasHeight = x)
      } text ("Height of final image (default 1000 units).")
      opt[Int]("image-width") action { (x, c) =>
        c.copy(canvasWidth = x)
      } text ("Width of final image (default 1000 units).")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.tree)
      drawTree(config)
    }
  }


  def drawTree(config: Config): Unit = {
    println(timeStamp + "Loading tree")
    //load tree
    val tree = loadKmerTree(config.tree)
    //get all leafs
    val all_leafs = tree.getLeafNames()
    //get total nodes
    val all_nodes = tree.id2Leafnames().map(_._1)
    //set node levels
    val node2levels = if(!config.equalDist) Map.empty[String,Int] else tree.node2Levels().toMap
    //set max node level
    val dist_step = {
      if(!config.equalDist) 0
      else {
        //set max node level
        val max_level = node2levels.maxBy(_._2)
        println(timeStamp + "Using max tree level of " + max_level._2 + " (" + max_level._1 + ")")
        //get distance level step
        config.canvasWidth.toDouble / max_level._2
      }
    }
    /**
      * Set font
      */
    val font = Font(Monospaced, Bold, Points(config.fontSize))
    /**
      * Function to normalize a given distance based on canvas width
      *
      * @return
      */
    def normalize_dist: (String,Double) => Double = (id,value) => {
      //use original branch distances
      if(!config.equalDist) config.canvasWidth - ((value / tree.sum_dist) * config.canvasWidth)
      //use equal distant branches
      else dist_step * node2levels(id)
    }

    //get colours, if provided
    val id2Colour = {
      //attempt to laod colours file
      val tmp = loadColoursFile(config.coloursFile)
      //log
      if(config.coloursFile != null) {
        println(timeStamp + "Found colours for " + tmp.size + " leafs")
        //get missing leafs, if any
        val missing_leafs = all_leafs.toSet.diff(tmp.keySet)
        //sanity check
        assert(missing_leafs.isEmpty, "Could not find colour mapping for the following IDs: " + missing_leafs
          .mkString(","))
      }
      tmp
    }
    //get proportions, if provided
    val node2Proportions = {
      if(config.proportions == null) all_nodes.map(x => (x, config.lineWidth.toDouble)).toMap
      else {
        //set min line width
        val min_width = config.minFreq * config.lineWidth
        println(timeStamp + "Setting min frequency to " + config.minFreq + " (" + min_width + ")")
        //load proportions from query file and ignore given min frequency
        val tmp = loadProportionsFile(config.proportions, config.lineWidth, tree)
          .map(x => if(x._2 >= min_width) x else (x._1, 0.0))
        println(timeStamp + "Loaded branch proportions for " + tmp.size + " nodes/leafs")
        tmp
      }
    }
    println(timeStamp + "Found tree with " + all_nodes.size + " nodes (" + all_leafs.size +
      " leafs) and a total distance of " + tree.sum_dist)
    println(timeStamp + "Computing coordinates")
    //obtain x,y coordinates for every x,y coordinate
    val node2Coords = getCoords(tree, (config.canvasWidth, config.canvasHeight), normalize_dist)
    println(timeStamp + "Drawing tree")
    //draw tree
    val tree_image = treeDrawer(tree, font, config.fontSize/2, node2Coords, node2Proportions, id2Colour)
    println(timeStamp + "Writing image to disk")
    tree_image.save(config.outputDir + "/" + config.prefix + ".pdf")
    println(timeStamp + "Successfully completed!")
  }
}
