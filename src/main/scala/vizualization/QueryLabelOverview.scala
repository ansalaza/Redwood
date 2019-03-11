package vizualization

import java.io.File

import doodle.core.{Color, Image, allBeside}
import doodle.core.Image.{rectangle, text}
import doodle.jvm.Java2DFrame._
import doodle.backend.StandardInterpreter._
import doodle.core.Color
import doodle.jvm.FileFrame.pdfSave
import doodle.core.font.Font
import doodle.core.font.FontFace.Bold
import doodle.core.font.FontFamily.Monospaced
import doodle.core.font.FontSize.Points
import doodle.syntax._
import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory, verifyFile}
import utilities.MetaDataUtils.{loadColoursFile, loadQueryLabelFile}

import scala.annotation.tailrec

/**
  * Author: Alex N. Salazar
  * Created on 30-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object QueryLabelOverview {

  case class Config(
                     queryFiles: File = null,
                     outputDir: File = null,
                     canvasHeight: Int = 1000,
                     canvasWidth: Int = 1000,
                     fontSize: Int = 14,
                     coloursFile: File = null,
                     prefix: String = null,
                     lineWidth: Int = 3
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("tree-tracing") {
      opt[File]('q', "query-files") required() action { (x, c) =>
        c.copy(queryFiles = x)
      } text ("File containing path to query labels, one per line.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      opt[File]("colours") required() action { (x, c) =>
        c.copy(coloursFile = x)
      } text ("Tab-delimited file containing label and colour as hue, saturation, and lightness values. Order in " +
        "visualization is inferred from order of labels in this file.")
      opt[String]("prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output file name.")
      note("\nOPTIONAL METADATA")
      opt[Int]("font-size") action { (x, c) =>
        c.copy(fontSize = x)
      } text ("Font-size of labels (default is 14).")
      opt[Int]("line-width") action { (x, c) =>
        c.copy(lineWidth = x)
      } text ("Line width for tree (default is 3).")
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
      verifyFile(config.queryFiles)
      queryLabelOverview(config)
    }

    def queryLabelOverview(config: Config): Unit = {
      //load label -> colours
      val label2colour = loadColoursFile(config.coloursFile) + ("Unique" -> Color.white)
      //create label indeces
      val label2index = {
        openFileWithIterator(config.coloursFile).toList.map(_.split("\t").head).zipWithIndex.toMap +
          ("Unique" -> (label2colour.size - 1))
      }
      //create list of sorted labels
      val sorted_labels = label2index.toList.sortBy(_._2).map(_._1)

      /**
        * Function to sort two given labelled-query maps
        * @return Boolean
        */
      def sortQueries: ((String, Map[String, Double]), (String,Map[String, Double])) => Boolean = (x,y) => {
        @tailrec def _sortQueries(labels: List[String]): Boolean = {
          labels match {
            case Nil => true
            case (head::tail) => {
              val (xval, yval) =(x._2(head),y._2(head))
              if(xval == yval) _sortQueries(tail) else xval > yval
            }
          }
        }
        _sortQueries(sorted_labels)
      }

      println(timeStamp + "Found " + label2colour.size + " total labels")
      //load query labels
      val queries = openFileWithIterator(config.queryFiles).toList.map(x => {
        //set as file object
        val file = new File(x)
        //verify
        verifyFile(file)
        (loadQueryLabelFile(file))
      }).sortWith(sortQueries)
      //set box heights
      val height = config.canvasHeight.toDouble / queries.size
      //set font
      val font = Font(Monospaced, Bold, Points(config.fontSize))
      println(timeStamp + "Processed " + queries.size + " query files")
      println(timeStamp + "Drawing")
      //iterate through each query and draw proportions
      val viz = queries.foldLeft((Image.empty, queries.size)){ case ((figure,index), (name, queries)) => {
        //sanity check
        assert(queries.size == label2colour.size, "Expected " + label2colour.size + " labels but found " + queries.size)
        (allBeside(queries.toList.sortBy(x => label2index(x._1)).map { case (label, width_prop) => {
          rectangle((config.canvasWidth * width_prop), height).lineWidth(config.lineWidth).fillColor(label2colour(label))
        }}).at(0, index * height).on(figure), index - 1)
      }}._1

      queries.map(_._1).foldLeft((Image.empty, queries.size)){ case ((drawing,index), name) => {
        (text(name).font(font).at(0, index * height).on(drawing), index-1)
      }}._1.beside(viz).save(config.outputDir + "/" + config.prefix + ".pdf")
      println(timeStamp + "Successfully completed!")
    }

  }
}
