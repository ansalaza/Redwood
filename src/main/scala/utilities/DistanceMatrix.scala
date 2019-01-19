package utilities

import java.io.File

import utilities.KmerTreeClusteringUtils._
import utilities.KmerTreeUtils.{Kmers, Node, Tree}

/**
  * Author: Alex N. Salazar
  * Created on 7-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
/**
    *
    * @param matrix
    * @param labels
    * @param pairwise_leaf_dist
    * @param sketch_map
    */
  class DistanceMatrix(val matrix: Matrix,
                       val labels: Labels,
                       val pairwise_leaf_dist: Map[(String, String), Double],
                       val sketch_map: Map[String, File] = Map()) {

  private val infinity = 1000000.0

    /**
      * Override to output matrix as a tab-delimited string
      * @return String
      */
    override def toString: String = {
      //set names of rows and columns
      val string_labels = labels.map(x => x._1.getNodeName())
      //set matrix entries
      val entries = string_labels.zip(matrix.map(_._1.map(_._1))).map(x => x._1 + "\t" + x._2.mkString("\t"))
      //create tab-delimited string
      "$" + "\t" + string_labels.mkString("\t") + "\n" + entries.mkString("\n")
    }

    /**
      * Method to find the pair of nodes in the matrix that are that have the smallest global distance
      *
      * @return 2-tuple ((Node, Node), distance)
      */
    def getMinPair(): (Double, ((Tree[Kmers], Int), (Tree[Kmers], Int))) = {
      //set fake initial minimum distance
      val initial_min = (infinity, (-1,-1))
      /**
        * Iterate through matrix and find min distance
        * Returns 2-tuple: (distance, (ith, jth))
        */
      val (min_dist, (row_index, column_index)) = {
        //iterate through rows using initial min distance
        matrix.foldLeft(initial_min) { case (min, (distances, i)) => {
          //iterate through columns
          distances.foldLeft(min) { case (acc_min, (distance, j)) => {
            //update min distance
            if (i != j && distance < acc_min._1) (distance, (i, j)) else acc_min
          }}
        }}
      }
      //return tuple
      (min_dist, ((labels(row_index)._1, row_index), (labels(column_index)._1, column_index)))
    }

    /**
      * Method to merge two given clusters (trees) and update labels (trees) and matrix with re-calculated
      * single-linkage distances
      *
      * @param i Tree (cluster) to merge
      * @param j Tree (cluster) to merge
      * @return 2-tuple: (updated matrix, updated labels)
      */
    def mergeClusters(i: (Tree[Kmers], Int), j: (Tree[Kmers], Int), d: Double, kmers: Kmers): (Matrix, Labels) = {
      /**
        * FUnction to determine whether a given index belongs to the two nodes being merged
        *
        * @return Boolean
        */
      def isOld: Int => Boolean = index => index == i._2 || index == j._2

      //create new array of labels by removing old nodes and appending new cluster to front of array wit updated index
      val new_labels = {
        //remove old labels
        val filtered = labels.filterNot(x => isOld(x._2)).map(_._1)
        //create new node cluster
        val new_cluster: Tree[Kmers] = {
          //place node with smallest lexicographic name to the left
          if (i._1.getNodeName() < j._1.getNodeName()) new Node(kmers, d,  i._1, j._1)
          else Node(kmers, d, j._1, i._1)
        }
        //append new cluster to filtered labels
        (new_cluster :: filtered.toList).toArray.zipWithIndex
      }
      //create new matrix by re-calculating distances using single-linkage
      val new_matrix = {
        //iterate through labels
        new_labels.foldLeft(List[List[Double]]())((acc_matrix, row) => {
          //fill new distances in the current row
          val updated_row = new_labels.foldLeft(List[Double]())((acc_distances, column) => {
            //same clusters
            if (column._2 == row._2) 0.0 :: acc_distances
            else {
              //compute dist between trees (clusters) using single-linkage
              val dist = singleLinkageDist(row._1, column._1)
              //append to list
              dist :: acc_distances
            }
          }).reverse
          //append row to list
          updated_row :: acc_matrix
        })
          //reverse row-order and reformat to Matrix
          .map(_.zipWithIndex.toArray).reverse.zipWithIndex.toArray
      }
      //return new matrix and new labels
      (new_matrix, new_labels)
    }

    /**
      * Function to compute single-linkage distance metric between two given trees (clusters)
      *
      * @return Double
      */
    private def singleLinkageDist: (Tree[Kmers], Tree[Kmers]) => Double = (t1, t2) => {
      //get all leafs in this tree
      val t1_leafs = t1.getLeafNames()
      //get all leafs in this tree
      val t2_leafs = t2.getLeafNames()
      //iterate through t1 with initial distance
      t1_leafs.foldLeft(infinity)((min, x) => {
        //iterate through t2
        t2_leafs.foldLeft(min)((acc_min, y) => {
          //the same leaf cannot be in two different clusters/trees
          assert(x != y, "Leaf " + x + " appears in two different clusters: " + (t1, t2))
          //fetch distance between the two leafs
          val dist = pairwise_leaf_dist((x, y))
          //update min distance
          if (dist < acc_min) dist else acc_min
        })
      })
    }

  //private def

    /**
      * Function to create name of a cluster by concatenating all leaf names in each node in lexicographic order
      *
      * @return String
      */
    private def getClusterName: (Tree[Kmers], Tree[Kmers]) => String = (i, j) =>
      (i.getLeafNames() ::: j.getLeafNames()).sorted.mkString(",")

  }

