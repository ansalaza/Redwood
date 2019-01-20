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
import doodle.jvm.FileFrame.svgSave
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
      note("\n\nOPTIONAL AESTHETICS")
      opt[Int]("image-height") action { (x, c) =>
        c.copy(canvasHeight = x)
      } text ("Height of final image (default 1000 units).")
      opt[Int]("image-width") action { (x, c) =>
        c.copy(canvasWidth = x)
      } text ("Width of final image (default 1000 units).")
      opt[Int]("font-size") action { (x, c) =>
        c.copy(fontSize = x)
      } text ("Font-size of labels (default is 25).")
      opt[Int]("line-width") action { (x, c) =>
        c.copy(lineWidth = x)
      } text ("Line width for feature boxes (default is \"font-size\"/4).")
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
    /**
      * Set font
      */
    val font = Font(Monospaced, Bold, Points(config.fontSize))
    /**
      * Function to normalize a given distance based on canvas width
      *
      * @return
      */
    def normalize_dist: Double => Double = value => config.canvasWidth - ((value / tree.sum_dist) * config.canvasWidth)
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
      val tmp = all_nodes.map(x => (x, config.lineWidth.toDouble)).toMap
      if(config.proportions == null) tmp else loadProportionsFile(config.proportions, 10, tree)
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
    tree_image.save(config.outputDir + "/" + config.prefix + ".svg")
    println(timeStamp + "Successfully completed!")
  }
}
