package utilities

import java.io.File

import doodle.core.{Color, Image}
import doodle.core.Image.{openPath, text, rectangle}
import doodle.core.Color.hsl
import doodle.syntax._
import doodle.core.PathElement.{lineTo, moveTo}
import doodle.core.font.Font
import utilities.FileHandling.openFileWithIterator
import utilities.KmerTreeUtils.{Kmers, Leaf, Node, Tree}

/**
  * Author: Alex N. Salazar
  * Created on 20-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object TreeDrawingUtils {

  //set default color of black
  private val default_color = hsl(0.degrees, 0.normalized, 0.normalized)

  /**
    * Method to draw tree given a tree, font specification, and a map of all node/leafs x,y coords
    *
    * @param tree
    * @param font
    * @param node2Coords
    * @return Image
    */
  def treeDrawer(tree: Tree[Kmers],
                 font: Font,
                 text_offset: Int,
                 node2Coords: Map[String, (Double, Double)],
                 node2Proportions: Map[String, Double] = Map(),
                 leaf2Colour: Map[String, Color] = Map()): Image = {

    def drawSubBranch(from: (Double, Double), to: (Double, Double), width: Double): Image = {
      if (width == 0.0) Image.empty
      else openPath(List(moveTo(from._1, from._2), lineTo(to._1, to._2))).lineWidth(width)
    }

    /**
      * Recursive method to traverse tree and draw node, branches, and leafs
      *
      * @param current
      * @param acc
      * @return Image
      */
    def _treeDrawer(current: Tree[Kmers], acc: Image): Image = {
      current match {
        case Leaf(a, b) => acc
        case Node(i, k, d, l, r) => {
          //update left right subtree
          val updated_image = _treeDrawer(l, acc).on(_treeDrawer(r, acc))
          //get current nodes prop
          val node_prop = node2Proportions(i.toString)
          //current node is empty, move on
          if (node_prop == 0.0) updated_image
          //current node is not empty
          else {
            //get x,y coordinates
            val (x, y) = node2Coords(i.toString)
            //updated with node
            val updated_with_node = (text(i.toString).font(font)).at(x + text_offset, y).on(updated_image)
            //get left and right line-widths
            val (left_width, right_width) = (node2Proportions(l.getID()), node2Proportions(r.getID()))
            //both sub-branches are empty, move on
            if (left_width == 0.0 && right_width == 0.0) updated_with_node
            //draw sub-branches
            else {
              //get left and right sub-tree x and y-coord
              val ((left_x, left_y), (right_x, right_y)) = (node2Coords(l.getID()), node2Coords(r.getID()))
              //draw left-branch
              val left_branch = {
                //draw horizontal
                drawSubBranch((x, right_y), (right_x, right_y), right_width).on(
                  //draw vertical
                  drawSubBranch((x, y), (x, right_y), right_width))
              }
              //draw right-branch
              val right_branch = {
                drawSubBranch((x, left_y), (left_x, left_y), left_width).on(
                  drawSubBranch((x, y), (x, left_y), left_width))
              }
              //update image
              left_branch.on(right_branch).on(updated_with_node)
            }
          }
        }
      }
    }


    //get root coordinates
    val (rx, ry) = node2Coords(tree.getID())
    //draw rooted line
    val rooted = openPath(List(moveTo(rx, ry), lineTo(rx - 50, ry))).lineWidth(node2Proportions(tree.getID()))

    rooted.beside(
      _treeDrawer(tree, Image.empty).beside(
        //iterate through leafs and draw sample names
        tree.getLeafNames().reverse.foldLeft(Image.empty)((acc, leaf) => {
          val leaf_prop = node2Proportions(leaf)
          //get x, y coordinates
          val y_coord = node2Coords(leaf)._2
          //get colour
          val colour = leaf2Colour.getOrElse(leaf, default_color)
          //draw leaf on y coord
          (if (leaf_prop == 0.0) Image.empty else (text(leaf).font(font).fillColor(colour))).at(0, y_coord).on(acc)
        })
      )
    )
  }

  /**
    * Obtain (x,y) coordinates for every node/leaf in the tree
    *
    * @param tree
    * @return Map[String, (Double,Double)]
    */
  def getCoords(tree: Tree[Kmers],
                canvas: (Int, Int),
                normalize: (String,Double) => Double
               ): Map[String, (Double, Double)] = {
    //get total number of leafs
    val total_leafs = tree.getLeafNames().size
    //get the vertical y-coordinate step based on total leafs
    val ycoord_step = canvas._2.toDouble / total_leafs
    //get all leafs and their corresponding y-coord
    val leaf2ycoord = tree.getLeafNames().zipWithIndex.map(x => (x._1, x._2 * ycoord_step)).toMap

    /**
      * Recursive method to traverse tree and obtain x and y-coordinates for nodes and leafs
      *
      * @param current
      * @param acc
      * @return
      */
    def _getCoords(current: Tree[Kmers],
                   acc: Map[String, (Double, Double)]): Map[String, (Double, Double)] = {
      //check current (sub)tree
      current match {
        //leaf, add x,y coord
        case Leaf(a, b) => acc + (b -> (canvas._1, leaf2ycoord(b)))
        //node
        case Node(i, k, d, l, r) => {
          //Update left and right sub-trees
          val updated_acc = _getCoords(l, acc) ++ _getCoords(r, acc)
          //set y-coord as average of left and right
          val y_coord = (updated_acc(l.getID())._2 + updated_acc(r.getID())._2) / 2
          //add x, y coord
          updated_acc + (i.toString -> (normalize(i.toString, d), y_coord))
        }
      }
    }

    _getCoords(tree, Map())
  }

  /**
    * Function to open a provided colours file
    *
    * @return Map[String, HSL]
    */
  def loadColoursFile: File => Map[String, Color] = file => {
    //no colours file, return empty map
    if (file == null) Map.empty[String, Color]
    //open
    else openFileWithIterator(file).toList.map(x => {
      //get columns
      val columns = x.split("\t")
      //set colour
      (columns.head,
        hsl(columns(1).toDouble.degrees, columns(2).toDouble.normalized, columns(3).toDouble.normalized))
    }).toMap
  }

  def loadProportionsFile(file: File, max_width: Int, tree: Tree[Kmers]): Map[String, Double] = {
    openFileWithIterator(file).drop(1).toList.foldLeft(Map[String, Double]())((acc, line) => {
      val columns = line.split("\t")
      acc + (columns.head -> (columns(1).toDouble * max_width))
    })
  }


}
