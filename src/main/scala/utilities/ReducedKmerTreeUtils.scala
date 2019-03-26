package utilities

import java.io.File

import utilities.FileHandling.deserialize
import boopickle.Default._
import utilities.KmerTreeUtils.Leaf

/**
  * Author: Alex N. Salazar
  * Created on 28-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object ReducedKmerTreeUtils {


  sealed trait ReducedTree {

    def getId2Distances(tree: ReducedTree = this): List[(Int, Double)] = {
      def _getId2Distances(current: ReducedTree, acc: List[(Int,Double)]): List[(Int,Double)] = {
        current match {
          case ReducedLeaf(i, n, d, g) => (i,d) :: acc
          case ReducedNode(i,d,l,r) => ((i,d) :: _getId2Distances(l, acc)) ::: _getId2Distances(r, acc)
        }
      }
      _getId2Distances(tree, List())
    }

    /**
      * Method to get lowest common ancestor for a given list of node/leaf IDs
      *
      * @param tree
      * @param query
      * @return Node ID as Int
      */
    def lowestCommonAncestor(tree: ReducedTree = this, query: List[Int]): Int = {
      def _lowestCommonAncestor(current: ReducedTree, acc: Int): Int = {
        if (current.isLeaf()) acc
        else {
          //get current node's id
          val current_id = current.getID()
          //check how many query IDs are less than current and also check if any of the ids are the same as current id
          val (lessthans, islowest) = query.foldLeft((Set[Boolean](), false)) { case ((acc, same), id) => {
            (acc + (id < current_id), if (same) same else id == current_id)
          }
          }
          //one of the ID is the local id or it is a mixture of greater/less than, this is the LCA
          if (islowest || lessthans.size != 1) acc
          //all are less than, go left
          else if (lessthans.head) {
            val left = current.loadAsNode().left
            _lowestCommonAncestor(left, left.getID())
          }
          //all are greater than, go right
          else {
            val right = current.loadAsNode().right
            _lowestCommonAncestor(right, right.getID())
          }
        }
      }

      _lowestCommonAncestor(tree, tree.loadAsNode().getID())
    }

    def getLeaf2Parent(tree: ReducedTree = this): Map[Int,Int] = {
      def _getLeaf2Parent(current: ReducedTree, acc: List[(Int,Int)]): List[(Int,Int)] = {
        current match {
          case ReducedLeaf(i, n, d, g) => acc
          case ReducedNode(i, d, l, r) => {
            _getLeaf2Parent(l, if(!l.isLeaf()) acc else (l.loadAsLeaf().id, i) :: acc) :::
            _getLeaf2Parent(r, if(!r.isLeaf()) acc else (r.loadAsLeaf().id, i) :: acc)
          }
        }
      }
      _getLeaf2Parent(tree, List()).toMap
    }

    /**
      * Method to obtain the tree level (depth) of all nodes/leafs
      *
      * @param tree
      * @return List[(String,Int)]
      */
    def getId2Levels(tree: ReducedTree = this): Map[Int, Int] = {
      def _getId2Levels(current: ReducedTree, level: Int, acc: List[(Int, Int)]): List[(Int, Int)] = {
        current match {
          case ReducedLeaf(i, n, d, g) => (i, level) :: acc
          case ReducedNode(i, d, l, r) =>
            (i, level) :: (_getId2Levels(l, level + 1, acc) ::: _getId2Levels(r, level + 1, acc))
        }
      }

      _getId2Levels(tree, 0, List()).toMap
    }

    def getLeafNamesPostOrder(tree: ReducedTree = this): List[String] = {
      /**
        * Recursive method to iterate through tree in post-order traversal
        *
        * @param current Current Tree
        * @param acc     Accumulating list of post-traversal
        * @return List of (reversed) post-order traversal
        */
      def _getLeafId2Name(current: ReducedTree, acc: List[String]): List[String] = {
        current match {
          //append value to list
          case ReducedLeaf(i, n, d, g) => n :: acc
          //traverse to left then right
          case ReducedNode(i, d, l, r) => (_getLeafId2Name(r, acc) ::: _getLeafId2Name(l, acc))
        }
      }

      _getLeafId2Name(tree, List()).reverse
    }

    def getNodeIDsPostOrder(tree: ReducedTree = this): List[Int] = {
      /**
        * Recursive method to iterate through tree in post-order traversal
        *
        * @param current Current Tree
        * @param acc     Accumulating list of post-traversal
        * @return List of (reversed) post-order traversal
        */
      def _getNodeIDsPostOrder(current: ReducedTree, acc: List[Int]): List[Int] = {
        current match {
          //append value to list
          case ReducedLeaf(i, n, d, g) => acc
          //traverse to left then right
          case ReducedNode(i, d, l, r) => i :: (_getNodeIDsPostOrder(r, acc) ::: _getNodeIDsPostOrder(l, acc))
        }
      }

      _getNodeIDsPostOrder(tree, List()).reverse
    }

    /**
      * Method to obtain leaf-names in post-order traversal of a given tree
      *
      * @param tree Tree of type Kmers
      * @return post-order traversal as a list of each node's name
      */
    def getLeafId2Name(tree: ReducedTree = this): Map[Int, String] = {
      /**
        * Recursive method to iterate through tree in post-order traversal
        *
        * @param current Current Tree
        * @param acc     Accumulating list of post-traversal
        * @return List of (reversed) post-order traversal
        */
      def _getLeafId2Name(current: ReducedTree, acc: List[(Int, String)]): List[(Int, String)] = {
        current match {
          //append value to list
          case ReducedLeaf(i, n, d, g) => (i,n) :: acc
          //traverse to left then right
          case ReducedNode(i, d, l, r) => (_getLeafId2Name(r, acc) ::: _getLeafId2Name(l, acc))
        }
      }

      _getLeafId2Name(tree, List()).toMap
    }
    /*
    /**
      * Method to obtain all leaf names under each node/leaf
      *
      * @param tree
      * @return List of 2-tuples: (node/leaf ID, leaf names in lexicographic order)
      */
    def getId2LeafNames(tree: ReducedTree = this): Map[Int, List[String]] = {
      /**
        * Recursive method to iterate through tree in post-order traversal
        *
        * @param current Current Tree
        * @param acc     Accumulating list of post-traversal
        * @return List of (reversed) post-order traversal
        */
      def _getId2LeafNames(current: ReducedTree, acc: List[(Int, List[String])]): List[(Int, List[String])] = {
        current match {
          //append value to list
          case ReducedLeaf(i, n, d, g) => (i, List(n)) :: acc
          //traverse to left then right
          case ReducedNode(i, d, l, r) =>
            (i, getLeafId2Name(current).map(_._2).toList.sorted) :: (_getId2LeafNames(r, acc) ::: _getId2LeafNames(l, acc))
        }
      }

      _getId2LeafNames(tree, List()).toMap
    }
    */

    def node2Children(tree: ReducedTree = this): Map[Int, List[Int]] = {
      def _node2Children(current: ReducedTree, acc: Map[Int, List[Int]]): Map[Int, List[Int]] = {
        current match {
          case ReducedLeaf(i, n, d, g) => acc + (i -> List(i))
          case ReducedNode(i, d, l, r) =>{
            val updated = _node2Children(l, acc) ++ _node2Children(r, acc)
            updated + (i -> updated.keys.toList)
          }
       }
      }
      _node2Children(tree, Map())
    }

    /**
      * Recursive method to compute a map of node ID -> (largest genome, smallest sketch).
      * The two values are based on all children under each node.
      *
      * @param tree
      * @return Map[Int,(Int,Int)]
      */
    def getId2LargestGenome(tree: ReducedTree = this): Map[Int, Int] = {
      def _getId2LargestGenome(current: ReducedTree, acc: List[(Int, Int)]): List[(Int, Int)] = {
        current match {
          case ReducedLeaf(i, n, d, g) => (i, g) :: acc
          case ReducedNode(i, d, l, r) => {
            //set largest left and right sub-trees
            val downstream = _getId2LargestGenome(l, acc) ::: _getId2LargestGenome(r, acc)
            //set largest genome
            val largest_genome = downstream.map(_._2).max
            //append
            (i, largest_genome) :: acc
          }
        }
      }

      _getId2LargestGenome(tree, List()).toMap
    }

    /**
      * Method to determine whether a given tree is a leaf
      *
      * @param tree
      * @return Boolean
      */
    def isLeaf(tree: ReducedTree = this): Boolean = tree match {
      case ReducedLeaf(i, n, d, g) => true;
      case _ => false
    }

    /**
      * Load tree as leaf
      *
      * @param tree
      * @return ReducedLeaf
      */
    def loadAsLeaf(tree: ReducedTree = this): ReducedLeaf = tree.asInstanceOf[ReducedLeaf]

    /**
      * Load tree as node
      *
      * @param tree
      * @return ReducedNode
      */
    def loadAsNode(tree: ReducedTree = this): ReducedNode = tree.asInstanceOf[ReducedNode]

    /**
      * Get ID of current tree node
      *
      * @param tree
      * @return
      */
    def getID(tree: ReducedTree = this): Int = if (tree.isLeaf()) tree.loadAsLeaf().id else tree.loadAsNode().id
  }

  /**
    * Leaf of reduced kmer tree
    * @param id
    * @param name
    * @param genome_size
    */
  case class ReducedLeaf(id: Int, name: String, dist: Double, genome_size: Int) extends ReducedTree

  /**
    * Node of reduced kmer tree
    *
    * @param id    Node ID
    * @param left  Left branch of the node as Tree[A]
    * @param right Right branch of the node as Tree[A]
    */
  case class ReducedNode(id: Int, dist: Double, left: ReducedTree, right: ReducedTree) extends ReducedTree

  /**
    * Reduced representation of a kmer-tree. The tree is a balanced-binary search tree storing the LCA for every kmer
    * in the tree as well as the universal sketch size  and kmer length from all the leafs. This representation is more
    * efficient in speed, memory, and storage as every kmer from the universal set of kmers is stored only once and
    * the LCA of every kmer is pre-computed.
    * @param tree
    * @param lca_map
    * @param min_sketch_size
    */
  case class ReducedKmerTree(tree: ReducedTree, lca_map: Map[Int, Int], min_sketch_size: Int, kmer_length: Int) {

    /**
      * Recursive method to obtain the cumulative counts for each node/leaf by summing the counts of all the children of
      * under each node.
      *
      * @param counts Map[Int,Int]
      * @return Map[Int,Int]
      */
    private def cumulativeCounts(counts: Map[Int, Double]): Map[Int, Double] = {
      def _cumulativeCounts(current: ReducedTree, acc: Map[Int, Double]): Map[Int, Double] = {
        current match {
          case ReducedLeaf(i, n, d, g) => acc + (i -> (counts.getOrElse(i, 0.0)))
          case ReducedNode(i, d, l, r) => {
            //get cumulative weights of children
            val updated_prop = _cumulativeCounts(l, acc) ++ _cumulativeCounts(r, acc)
            //calculate current node's weight
            val node_weight = counts.getOrElse(i, 0.0) + updated_prop(l.getID()) + updated_prop(r.getID())
            //update acc
            updated_prop + (i -> node_weight)
          }
        }
      }

      _cumulativeCounts(tree, Map())
    }

    /**
      * Function to compute the cumulative sum of each node/leaf for some generic count map as node/leaf ID -> count
      * by summing the counts of all the children under each node.
      *
      * @return Map[Int,Int]
      */
    def genericCumulative: Map[Int, Double] => Map[Int, Double] = counts => cumulativeCounts(counts)

    /**
      * Obtain cumulative number of kmers under each node/leaf
      */
    lazy val node2totalkmers = cumulativeCounts(lca_map.map(_._2).groupBy(identity).mapValues(_.size))

  }

  /**
    * Function to deserialize a reduced kmer tree from a given file
    *
    * @return ReducedKmerTree
    */
  def loadReducedKmerTree: File => ReducedKmerTree = file => Unpickle[ReducedKmerTree].fromBytes(deserialize(file))

}
