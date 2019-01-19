package utilities

import java.io.File
import java.nio.ByteBuffer
import java.nio.file.{Files, Paths}

/**
  * Author: Alex N. Salazar
  * Created on 4-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object KmerTreeUtils {

  type Kmers = Set[Int]

  /**
    *
    * @tparam Kmers
    */
  sealed trait Tree[+Kmers] {
    /**
      * Method to obtain leaf-names in post-order traversal of a given tree
      *
      * @param tree Tree of type Kmers
      * @return post-order traversal as a list of each node's name
      */
    def getLeafNames[Kmers](tree: Tree[Kmers] = this): List[String] = {
      /**
        * Recursive method to iterate through tree in post-order traversal
        *
        * @param current Current Tree
        * @param acc     Accumulating list of post-traversal
        * @return List of (reversed) post-order traversal
        */
      def _getLeafNames[Kmers](current: Tree[Kmers], acc: List[String]): List[String] = {
        current match {
          //append value to list
          case Leaf(a, b) => b :: acc
          //traverse to left then right
          case Node(a, d, l, r) => (_getLeafNames(r, acc) ::: _getLeafNames(l, acc))
        }
      }

      _getLeafNames(tree, List()).reverse
    }

    def getKmerSetSizes[Kmers](tree: Tree[Kmers] = this): List[(String, Int)] = {
      /**
        * Recursive method to iterate through tree in post-order traversal
        *
        * @param current Current Tree
        * @param acc     Accumulating list of post-traversal
        * @return List of (reversed) post-order traversal
        */
      def _getKmerSetSizes[Kmers](current: Tree[Kmers], acc: List[(String, Int)]): List[(String, Int)] = {
        current match {
          //append value to list
          case Leaf(a, b) => (b, a.size) :: acc
          //traverse to left then right
          case Node(a, d, l, r) =>
            (getNodeName(current), a.size) :: (_getKmerSetSizes(r, acc) ::: _getKmerSetSizes(l, acc))
        }
      }

      _getKmerSetSizes(tree, List()).reverse
    }

    def queryLCA[Kmers](kmer: Int, tree: Tree[Kmers] = this): List[String] = {
      def _queryLCA[Kmers](current: Tree[Kmers], acc: List[String]): List[String] = {
        current match {
          //append value to list
          case Leaf(a, b) => if (!b.contains(kmer)) acc else List(b)
          //traverse to left then right
          case Node(a, d, l, r) => {
            //check membership in left and right subtree
            val (l_contains, r_contains) = (l.data().contains(kmer), r.data().contains(kmer))
            //both contain
            if (l_contains && r_contains) acc
            //left contains kmer
            else if (l_contains && !r_contains) _queryLCA(l, l.getLeafNames())
            //right contains kmer
            else if (!l_contains && r_contains) _queryLCA(r, r.getLeafNames())
            //neither contain kmer
            else acc
          }
        }
      }
      //check whether root contains kmer
      if (!tree.data().contains(kmer)) List() else _queryLCA(tree, tree.getLeafNames())
    }

    /**
      * Function to create name of a given node by concatanating all leaf names in lexicographic order
      *
      * @param tree
      * @tparam Kmers
      * @return
      */
    def getNodeName[Kmers](tree: Tree[Kmers] = this): String = tree.getLeafNames().sorted.mkString(",")

    def isLeaf[Kmers](tree: Tree[Kmers] = this): Boolean = tree match {
      case Leaf(a, b) => true;
      case _ => false
    }

    def loadAsNode[Kmers](tree: Tree[Kmers] = this): Node = tree.asInstanceOf[Node]

    def loadAsLeaf[Kmers](tree: Tree[Kmers] = this): Leaf = tree.asInstanceOf[Leaf]

    def data[Kmers](tree: Tree[Kmers] = this): KmerTreeUtils.Kmers =
      if (tree.isLeaf()) tree.loadAsLeaf().data else tree.loadAsNode().data
  }


  /**
    * Leaf node of a Tree
    *
    * @param data value of leaf as type A
    */
  case class Leaf(data: Kmers, id: String) extends Tree[Kmers]

  /**
    * Node of a tree
    *
    * @param data  Value of a node as type A
    * @param left  Left branch of the node as Tree[A]
    * @param right Right branch of the node as Tree[A]
    */
  case class Node(data: Kmers, sum_dist: Double, left: Tree[Kmers], right: Tree[Kmers]) extends Tree[Kmers]


  def loadKmerTree: File => Node = file => {
    import boopickle.Default._
    Unpickle[Node].fromBytes(ByteBuffer.wrap(Files.readAllBytes(Paths.get(file.getAbsolutePath))))
  }

}
