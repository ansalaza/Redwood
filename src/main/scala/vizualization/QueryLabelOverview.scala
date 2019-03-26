package vizualization

import java.io.{File, PrintWriter}

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
import utilities.NumericalUtils.round
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
                     queryOrder: File = null,
                     normUnique: Boolean = false,
                     minProp: Double = 0.0,
                     lineWidth: Int = 2
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
      note("\nOPTIONAL")
      opt[Unit]("normalize") action { (x, c) =>
        c.copy(normUnique = true)
      } text ("Normalize query proportions to their non-unique proportion.")
      opt[Double]("min-proportion") action { (x, c) =>
        c.copy(minProp = x)
      } text ("Minimum required value of a query's non-unique proportion. Query's vizualization is left blank if " +
        "below this value (default is 0.0).")
      opt[File]("query-order") action { (x, c) =>
        c.copy(queryOrder = x)
      } text ("File containing query ID (one per line) representing order in which each query is drawn.")
      opt[Int]("font-size") action { (x, c) =>
        c.copy(fontSize = x)
      } text ("Font-size of labels (default is 14).")
      opt[Int]("line-width") action { (x, c) =>
        c.copy(lineWidth = x)
      } text ("Line width for tree (default is 2).")
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
      if(config.queryOrder != null) verifyFile(config.queryOrder)
      verifyFile(config.coloursFile)
      queryLabelOverview(config)
    }

    def queryLabelOverview(config: Config): Unit = {
      //load label -> colours
      val label2colour = {
        val tmp = loadColoursFile(config.coloursFile)
        if (config.normUnique) tmp else tmp + ("Unique" -> Color.white)
      }
      //create label indeces
      val label2index = {
        val tmp = openFileWithIterator(config.coloursFile).toList.map(_.split("\t").head).zipWithIndex.toMap
        if (config.normUnique) tmp else tmp + ("Unique" -> (label2colour.size - 1))
      }
      //create list of sorted labels
      val sorted_labels = label2index.toList.sortBy(_._2).map(_._1)
      //get query order, if any
      val query2Index = {
        if (config.queryOrder == null) Map[String, Int]()
        else openFileWithIterator(config.queryOrder).toList.zipWithIndex.toMap
      }
      println(timeStamp + "Found " + label2colour.size + " total labels")
      //load query labels
      val queries = {
        //construct list of 2-tuples as: (name, list of proportions)
        val tmp = openFileWithIterator(config.queryFiles).toList.map(x => {
          //set as file object
          val file = new File(x)
          //verify
          verifyFile(file)
          (loadQueryLabelFile(file))
        }).map { case (name, props) => {
          //get non-unique proportion
          val non_unique = (props - ("Unique")).map(_._2).sum
          //below minimum non-unique threshold, leave as blank
          if (non_unique < config.minProp)
            (name, (if(!config.normUnique) props else props - ("Unique")).mapValues(x => 0.0))
          //use proportions as is
          else if (!config.normUnique) (name, props)
          //normalize proportions to
          else (name, (props - ("Unique")).mapValues(_ / non_unique))
        }
        }
        //order accordingly
        if (query2Index.isEmpty) tmp else tmp.sortBy(x => query2Index(x._1))
      }
      //set box heights
      val height = config.canvasHeight.toDouble / queries.size
      println(height)
      //set font
      val font = Font(Monospaced, Bold, Points(config.fontSize))

      /**
        * Function to draw a proportion given a proportion value and a label name
        *
        * @return Image
        */
      def drawProp: (Double, String) => Image = (prop, label) => {
        //set colour
        val colour = {
          if (label2colour.contains(label)) label2colour(label)
          else {
            if (!config.normUnique) println(timeStamp + "WARNING: could not find colour for label " + label)
            Color.white
          }
        }
        rectangle((config.canvasWidth*prop), height).lineWidth(config.lineWidth).fillColor(colour)
      }

      println(timeStamp + "Processed " + queries.size + " query files")
      println(timeStamp + "Drawing")
      //iterate through each query and draw proportions
      val viz = queries.foldLeft((Image.empty, queries.size)) { case ((figure, index), (name, queries)) => {
        println(timeStamp + "--" + name)
        if (queries.forall(_._2 == 0.0)) (drawProp(1.0, "unknown").at(0, index * height).on(figure), index - 1)
        else {
          //sanity check
          assert(queries.size == label2colour.size, "Expected " + label2colour.size + " labels but found " + queries.size)
          //compute total proportion
          val sum = round(queries.map(_._2).sum)
          assert(sum == 1.0, "Expected total proportion sum to be equal to 1.0: " + (name, sum))
          //iterate through each label and draw
          (allBeside(queries.toList.sortBy(x => label2index(x._1)).map { case (label, width_prop) => {
            //draw proportion with corresponding color
            drawProp(width_prop, label)
          }
          }).at(0, index * height).on(figure), index - 1)
        }
      }
      }._1
      //in same manner, add names
      queries.map(_._1).foldLeft((Image.empty, queries.size)) { case ((drawing, index), name) => {
        (text(name).font(font).at(0, index * height).on(drawing), index - 1)
      }
      }._1.beside(viz).save(config.outputDir + "/" + config.prefix + ".pdf")
      println(timeStamp + "Writing queries to single matrix")
      //set output file for matrix of proportions
      val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".query.matrix")
      //set column names
      pw.println(label2index.toList.sortBy(_._2).map(_._1).mkString("\t"))
      //output matrix of proportions
      queries.foreach { case (name, props) => {
        pw.print(name + "\t")
        if (props.forall(_._2 == 0.0)) pw.print(label2colour.toList.map(x => 0).mkString("\t"))
        else pw.print(props.toList.sortBy(x => label2index(x._1)).map(_._2).mkString("\t"))
        pw.println
      }
      }
      pw.close
      println(timeStamp + "Successfully completed!")

    }

  }
}
