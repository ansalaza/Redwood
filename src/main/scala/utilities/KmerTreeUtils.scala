package utilities

import java.io.File
import boopickle.Default._

import utilities.FileHandling.deserialize

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
      * Method to obtain the tree level (depth) of all nodes/leafs
      * @param tree
      * @tparam Kmers
      * @return List[(String,Int)]
      */
    def node2Levels[Kmers](tree: Tree[Kmers] = this): List[(String, Int)] = {
      def _node2Levels[Kmers](current: Tree[Kmers], level: Int, acc: List[(String,Int)]): List[(String,Int)] = {
        current match {
          case Leaf(a,b) => (b,level) :: acc
          case Node(i,a,d,l,r) =>
            (i.toString, level) :: (_node2Levels(l, level+1, acc) ::: _node2Levels(r, level+1, acc))
        }
      }
      _node2Levels(tree, 0, List())
    }

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
          case Node(i, a, d, l, r) => (_getLeafNames(r, acc) ::: _getLeafNames(l, acc))
        }
      }

      _getLeafNames(tree, List()).reverse
    }

    /**
      * Method to obtain all leaf names under each node
      *
      * @param tree
      * @tparam Kmers
      * @return List of 2-tuples: (node/leaf ID, leaf names in lexicographic order)
      */
    def id2Leafnames[Kmers](tree: Tree[Kmers] = this): List[(String, List[String])] = {
      /**
        * Recursive method to iterate through tree in post-order traversal
        *
        * @param current Current Tree
        * @param acc     Accumulating list of post-traversal
        * @return List of (reversed) post-order traversal
        */
      def _id2Leafnames[Kmers](current: Tree[Kmers],
                               acc: List[(String, List[String])]): List[(String, List[String])] = {
        current match {
          //append value to list
          case Leaf(a, b) => (b, List(b)) :: acc
          //traverse to left then right
          case Node(i, a, d, l, r) =>
            (i.toString(), getLeafNames(current).sorted) :: (_id2Leafnames(r, acc) ::: _id2Leafnames(l, acc))
        }
      }

      _id2Leafnames(tree, List()).reverse
    }

    /**
      * Method to traverse to perform post-order traversal and get kmer-set size for each node/leaf
      *
      * @param tree Tree
      * @tparam Kmers
      * @return List of 2-tuples: (ID, Kmer-set size)
      */
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
          case Node(i, a, d, l, r) => (i.toString, a.size) :: (_getKmerSetSizes(r, acc) ::: _getKmerSetSizes(l, acc))
        }
      }

      _getKmerSetSizes(tree, List()).reverse
    }

    /**
      * Method to obtain the total distance of a given (sub)tree
      *
      * @param tree
      * @tparam Kmers
      * @return Double
      */
    def totalDistance[Kmers](tree: Tree[Kmers] = this): Double = {
      def _totalDistance[Kmers](current: Tree[Kmers], acc: List[Double]): List[Double] = {
        current match {
          //append distance to list
          case Leaf(a, b) => acc
          case Node(i, a, d, l, r) => d :: (_totalDistance(l, acc) ::: _totalDistance(r, acc))
        }
      }
      _totalDistance(tree, List()).sum
    }

    /**
      * Method to obtain the lowest-common ancestor of a given kmer based on set-membership of left/right sub-trees
      *
      * @param kmer Kmer hash
      * @param tree Tree
      * @tparam Kmers
      * @return Node/leaf ID of the lowest common ancestor as Option[String]
      */
    def queryLCA[Kmers](kmer: Int, tree: Tree[Kmers] = this): Option[String] = {
      /**
        * Recursive method to traverse tree based on set-membership of left/right subtrees
        *
        * @param current Current node/leaf in tree
        * @param acc     current lowest common ancestor
        * @tparam Kmers
        * @return Node/leaf ID of lowest common ancestor
        */
      def _queryLCA[Kmers](current: Tree[Kmers], acc: String): String = {
        current match {
          //append value to list
          case Leaf(a, b) => if (!a.contains(kmer)) acc else b
          //traverse to left then right
          case Node(i, a, d, l, r) => {
            //sanity check
            assert(a(kmer))
            //check membership in left and right subtree
            val (l_contains, r_contains) = (l.loadKmers().contains(kmer), r.loadKmers().contains(kmer))
            //both contain
            if (l_contains && r_contains) acc
            //left contains kmer
            else if (l_contains && !r_contains) _queryLCA(l, i.toString)
            //right contains kmer
            else if (!l_contains && r_contains) _queryLCA(r, i.toString)
            //neither contain kmer
            else acc
          }
        }
      }
      //check whether root contains kmer
      if (!tree.loadKmers().contains(kmer)) None else Option(_queryLCA(tree, tree.loadAsNode().id.toString))
    }

    /**
      * Method to determine whether given tree is node or leaf
      *
      * @param tree
      * @tparam Kmers
      * @return Boolean
      */
    def isLeaf[Kmers](tree: Tree[Kmers] = this): Boolean = tree match {
      case Leaf(a, b) => true;
      case _ => false
    }

    /**
      * Load tree as node
      *
      * @param tree
      * @tparam Kmers
      * @return
      */
    def loadAsNode[Kmers](tree: Tree[Kmers] = this): Node = tree.asInstanceOf[Node]

    /**
      * Load tree as leaf
      *
      * @param tree
      * @tparam Kmers
      * @return
      */
    def loadAsLeaf[Kmers](tree: Tree[Kmers] = this): Leaf = tree.asInstanceOf[Leaf]

    /**
      * Get kmers from tree
      *
      * @param tree
      * @tparam Kmers
      * @return
      */
    def loadKmers[Kmers](tree: Tree[Kmers] = this): KmerTreeUtils.Kmers =
      if (tree.isLeaf()) tree.loadAsLeaf().kmers else tree.loadAsNode().kmers

    /**
      * Method to get ID of given tree
      *
      * @param tree
      * @tparam Kmers
      * @return
      */
    def getID[Kmers](tree: Tree[Kmers] = this): String =
      if (tree.isLeaf()) tree.loadAsLeaf().id else tree.loadAsNode().id.toString

  }


  /**
    * Leaf node of a Tree
    *
    * @param kmers value of leaf as type A
    */
  case class Leaf(kmers: Kmers, id: String) extends Tree[Kmers]

  /**
    * Node of a tree
    *
    * @param id    Node ID
    * @param kmers Value of a node as type A
    * @param left  Left branch of the node as Tree[A]
    * @param right Right branch of the node as Tree[A]
    */
  case class Node(id: Int, kmers: Kmers, sum_dist: Double, left: Tree[Kmers], right: Tree[Kmers]) extends Tree[Kmers]


  def loadKmerTree: File => Node = file => Unpickle[Node].fromBytes(deserialize(file))
    //Unpickle[Node].fromBytes(ByteBuffer.wrap(Files.readAllBytes(Paths.get(file.getAbsolutePath))))
    //deserialize(Files.readAllBytes(Paths.get(file.getAbsolutePath))).asInstanceOf[Node]

}
