package utilities

import java.io.File

import doodle.core.{Color, Image}
import doodle.core.Image.{openPath, text}
import doodle.core.Color.hsl
import doodle.syntax._
import doodle.core.PathElement.{lineTo, moveTo}
import doodle.core.font.Font
import utilities.FileHandling.openFileWithIterator
import utilities.ReducedKmerTreeUtils._

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
  def treeDrawer(tree: ReducedTree,
                 font: Font,
                 text_offset: Int,
                 node2Coords: Map[Int, (Double, Double)],
                 node2Proportions: Map[Int, Double] = Map(),
                 node2Colour: Map[Int, Color] = Map()): Image = {

    def drawSubBranch(from: (Double, Double),
                      to: (Double, Double),
                      width: Double,
                      color: Color): Image = {
      if (width == 0.0) Image.empty
      else {
        openPath(List(moveTo(from._1, from._2), lineTo(to._1, to._2)))
          .lineWidth(width).fillColor(color).lineColor(color)
      }
    }

    /**
      * Recursive method to traverse tree and draw node, branches, and leafs
      *
      * @param current
      * @param acc
      * @return Image
      */
    def _treeDrawer(current: ReducedTree, acc: Image): Image = {
      current match {
        case ReducedLeaf(a,b,c) => acc
        case ReducedNode(i, d, l, r) => {
          //update left right subtree
          val updated_image = _treeDrawer(l, acc).on(_treeDrawer(r, acc))
          //get current nodes prop
          val node_prop = node2Proportions(i)
          //current node is empty, move on
          if (node_prop == 0.0) updated_image
          //current node is not empty
          else {
            //get x,y coordinates
            val (x, y) = node2Coords(i)
            //get color
            val node_colour = node2Colour.getOrElse(i, default_color)
            //updated with node
            val updated_with_node =
              (text(i.toString).font(font)).fillColor(node_colour).at(x + text_offset, y).on(updated_image)
            //get left and right line-widths
            val (left_width, right_width) = (node2Proportions(l.getID()), node2Proportions(r.getID()))
            //both sub-branches are empty, move on
            if (left_width == 0.0 && right_width == 0.0) updated_with_node
            //draw sub-branches
            else {
              //get left and right sub-tree x and y-coord
              val ((left_x, left_y), (right_x, right_y)) = (node2Coords(l.getID()), node2Coords(r.getID()))
              //get left and right colours
              val (left_colour, right_colour) =
                (node2Colour.getOrElse(l.getID(), default_color), node2Colour.getOrElse(r.getID(), default_color))
              //draw right-branch
              val right_branch = {
                //draw horizontal
                drawSubBranch((x, right_y), (right_x, right_y), right_width, right_colour).on(
                  //draw vertical
                  drawSubBranch((x, y), (x, right_y), right_width, right_colour))
              }
              //draw left-branch
              val left_branch = {
                drawSubBranch((x, left_y), (left_x, left_y), left_width, left_colour).on(
                  drawSubBranch((x, y), (x, left_y), left_width, left_colour))
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
        tree.getLeafId2Name().foldLeft(Image.empty){case(acc, (id,name)) => {
          val leaf_prop = node2Proportions(id)
          //get x, y coordinates
          val y_coord = node2Coords(id)._2
          //get colour
          val colour = node2Colour.getOrElse(id, default_color)
          //draw leaf on y coord
          (if (leaf_prop == 0.0) Image.empty else (text(name).font(font).fillColor(colour))).at(0, y_coord).on(acc)
        }}
      )
    )
  }

  /**
    * Obtain (x,y) coordinates for every node/leaf in the tree
    *
    * @param tree
    * @return Map[String, (Double,Double)]
    */
  def getCoords(tree: ReducedTree,
                canvas: (Int, Int),
                normalize: (Int, Double) => Double
               ): Map[Int, (Double, Double)] = {
    //get total number of leafs
    val total_leafs = tree.getLeafId2Name().size
    //get the vertical y-coordinate step based on total leafs
    val ycoord_step = canvas._2.toDouble / total_leafs
    //get all leafs and their corresponding y-coord
    val leaf2ycoord = tree.getLeafNamesPostOrder().zipWithIndex.map(x => (x._1, x._2 * ycoord_step)).toMap

    /**
      * Recursive method to traverse tree and obtain x and y-coordinates for nodes and leafs
      *
      * @param current
      * @param acc
      * @return
      */
    def _getCoords(current: ReducedTree,
                   acc: Map[Int, (Double, Double)]): Map[Int, (Double, Double)] = {
      //check current (sub)tree
      current match {
        //leaf, add x,y coord
        case ReducedLeaf(a,b,c) => acc + (a -> (canvas._1, leaf2ycoord(b)))
        //node
        case ReducedNode(i, d, l, r) => {
          //Update left and right sub-trees
          val updated_acc = _getCoords(l, acc) ++ _getCoords(r, acc)
          //set y-coord as average of left and right
          val y_coord = (updated_acc(l.getID())._2 + updated_acc(r.getID())._2) / 2
          //add x, y coord
          updated_acc + (i -> (normalize(i, d), y_coord))
        }
      }
    }

    _getCoords(tree, Map())
  }
}
