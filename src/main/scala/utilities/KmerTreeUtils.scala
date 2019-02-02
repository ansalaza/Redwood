package utilities

import java.io.File

import utilities.ReducedKmerTreeUtils._

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
      * Method to compute the in-order traversal for a given tree
      * @param tree
      * @tparam Kmers
      * @return List of IDs
      */
    def inOrderTraversal[Kmers](tree: Tree[Kmers] = this): List[String] = {
      def _inOrderTraversal[Kmers](current: Tree[Kmers], acc: List[String]): List[String] = {
        current match {
          case Leaf(a, b, c, d) => b :: acc
          case Node(i, a, d, l, r) => {
            val left = (_inOrderTraversal(l, acc))
            val right = (_inOrderTraversal(r, acc))
            right ::: ((i.toString) :: left)
          }
        }
      }

      _inOrderTraversal(tree, List()).reverse
    }

    /**
      * Method to create a ReducedTree along with the LCA-kmer map
      * @param tree
      * @param old2new
      * @tparam Kmers
      * @return 2-tuple: (ReducedTree, LCA-kmer map)
      */
    def copy2ReducedKmerTree[Kmers](tree: Tree[Kmers] = this,
                                    old2new: Map[String, Int]
                                   ): (ReducedTree, Map[Int, Int]) = {

      def _copy2ReducedKmerTree[Kmers](current: Tree[Kmers],
                                       acc_tree: ReducedTree,
                                       acc_map: Map[Int, List[Int]]
                                      ): (ReducedTree, Map[Int, List[Int]]) = {

        /**
          * Function to add kmers to a kmer map given a kmer-map, kmers, and an ID
          * @return Map[Int,Lis[Int]
          */
        def addKmers:(Map[Int,List[Int]], Set[Int], Int) => Map[Int, List[Int]] = (map, kmers, id) => {
          kmers.foldLeft(map)((acc,kmer) => {
            val updated = id :: acc.getOrElse(kmer, List[Int]())
            acc + (kmer -> updated)
          })
        }
        current match {
          //node is a leaf
          case Leaf(a,b,c,d) => {
            //leaf id
            val id = old2new(b)
            //set reduced leaf node
            val leaf_node: ReducedTree = ReducedLeaf(id, b, c)
            //update kmer map with current leaf kmers
            (leaf_node, addKmers(acc_map, a, id))
          }
          //node is a node
          case Node(i, a, d, l, r) => {
            //get node id
            val id = old2new(i.toString)
            //set left and right subtree with updated map
            val (left, right, updated_map) = {
              //first set left subtree
              val tmp_l = _copy2ReducedKmerTree(l, acc_tree, acc_map)
              //then set right subtree using updated map from left
              val tmp_r = _copy2ReducedKmerTree(r, acc_tree, tmp_l._2)
              //return left and right with updated map from right
              (tmp_l._1, tmp_r._1, addKmers(tmp_r._2, a, id))
            }
            (ReducedNode(id, d, left, right), updated_map)
          }
        }
      }
      val (reduced_tree, kmer_map) = _copy2ReducedKmerTree(tree, null, Map())
      (reduced_tree,
        kmer_map.mapValues(x => if(x.size == 1) x.head else reduced_tree.lowestCommonAncestor(query = x)))
    }


    /**
      * Method to determine whether given tree is node or leaf
      *
      * @param tree
      * @tparam Kmers
      * @return Boolean
      */
    def isLeaf[Kmers](tree: Tree[Kmers] = this): Boolean = tree match {
      case Leaf(a, b, c, d) => true;
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
  case class Leaf(kmers: Kmers, id: String, genome_size: Int, sketch_size: Int) extends Tree[Kmers]

  /**
    * Node of a tree
    *
    * @param id    Node ID
    * @param kmers Value of a node as type A
    * @param left  Left branch of the node as Tree[A]
    * @param right Right branch of the node as Tree[A]
    */
  case class Node(id: Int, kmers: Kmers, sum_dist: Double, left: Tree[Kmers], right: Tree[Kmers]) extends Tree[Kmers]

  //Unpickle[Node].fromBytes(ByteBuffer.wrap(Files.readAllBytes(Paths.get(file.getAbsolutePath))))
  //deserialize(Files.readAllBytes(Paths.get(file.getAbsolutePath))).asInstanceOf[Node]

}
