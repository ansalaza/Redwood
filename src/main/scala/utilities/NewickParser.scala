package utilities

import scala.annotation.tailrec

/**
  * Author: Alex N. Salazar
  * Created on 7-3-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object NewickParser {

  /**
    * Case class for storing/parsing newick strings
    * @param str Newick string representing curreent (sub-)tree(s)
    * @param distance Distance of current node
    */
  case class NwkString(str: String, distance: Double) {
    /**
      * Function to determine whether current newick string represents a leaf
      * @return Boolean
      */
    def isLeaf(): Boolean = str.forall(x => x != '(' && x != ')' && x != ':')

    /**
      * Load subtree as new NwkString object
      * @return NwkString
      */
    def loadSubtree: String => NwkString = entry => {
      //split at distance char indicator
      val tmp = entry.splitAt(entry.lastIndexOf(":"))
      new NwkString(tmp._1, Math.abs(tmp._2.substring(1).toDouble))
    }

    /**
      * Method to load the next sub-trees
      * @return 2-tuple of: (NwkString, NwkString)
      */
    def getNextSubtrees(): (NwkString, NwkString) = {
      /**
        * Tail-recursive method to parse sub-trees in current newick string. Assumes binary sub-trees
        * @param remaining Remaining characters in newick string
        * @param count Count representing cumulative sum of '(' and ')' chars
        * @param acc Accumulating characters for sub-tree being parsed
        * @param parsed All parsed sub-trees
        * @return 2-tuple: (NwkString, NwkString)
        */
      @tailrec def _getNextSubtree(remaining: List[(Char, Int)],
                          count: Int,
                          acc: List[Char],
                          parsed: List[String]): (NwkString, NwkString) = {
        remaining match {
            //no more chars to process
          case Nil => {
            //no accumulating sub-tree
            if (acc.isEmpty) {
              //assert binary subtrees
              assert(parsed.size == 2, "Expected binary subtree but found " + parsed.size + " subtrees in tree:\n" +
              parsed.mkString("\n")
              )
              (loadSubtree(parsed(0)), loadSubtree(parsed(1)))
            }
            //accumulating sub-tree to process
            else {
              //concatenate subtree and add to parsed list
              val t = acc.reverse.mkString("") :: parsed
              assert(t.size == 2, "Expected binary subtree but found " + t.size + " subtrees in tree:\n" +
              t.mkString("\n"))
              (loadSubtree(t(0)), loadSubtree(t(1)))
            }
          }
          case ((char, index) :: tail) => {
            //finished parsing a sub-tree, add to parsed list
            if (count == 0 && char == ',') _getNextSubtree(tail, count, List(), acc.reverse.mkString("") :: parsed)
            //entering new sub-tree
            else if (char == '(') _getNextSubtree(tail, count + 1, char :: acc, parsed)
            //exiting sub-tree
            else if (char == ')') _getNextSubtree(tail, count - 1, char :: acc, parsed)
            //regular info inside subtree
            else _getNextSubtree(tail, count, char :: acc, parsed)
          }
        }
      }
      val start = if(str.startsWith("(") && str.endsWith(")")) str.drop(1).dropRight(1) else str
      _getNextSubtree(start.zipWithIndex.toList, 0, List(), List())
    }
  }

}
