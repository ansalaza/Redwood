package utilities

import java.io.File

import utilities.ClusteringUtils._
import utilities.KmerTreeUtils.{Kmers, Leaf, Node, Tree}
import utilities.SketchUtils.loadRedwoodSketch

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
  * @param clusters
  */
class Dendogram(val matrix: Matrix, val clusters: Clusters, val leaf2sketches: Map[String, File]) {

  /**
    * Function to obtain distance of a given pair of IDs aware that the matrix only stores distances in upper-diagonal
    * @return Double
    */
  private def fetchDist: (String, String) => Double = (x,y) =>
    if(x == y) 0.0 else if(matrix.contains((x,y))) matrix((x,y)) else matrix((y,x))

  /**
    * Override to output matrix as a tab-delimited string
    *
    * @return String
    */
  override def toString(): String = {
    //get ids for rows and clusters
    val index2ID = clusters.map(_._1).zipWithIndex.map(_.swap).toMap
    //set total ids
    val total_ids = index2ID.keys.max
    val entries = {
      //iterate through each row
      (0 to total_ids).foldLeft("")((string, i) => {
        //set row id
        val row_id = index2ID(i)
        //iterate through each column and add distances
        (0 to total_ids).foldLeft(string + row_id)((acc, j) => acc + "\t" + fetchDist(row_id, index2ID(j))) + "\n"
      })
    }
    //create tab-delimited string
    "$" + "\t" + index2ID.toList.sortBy(_._1).map(_._2).mkString("\t") + "\n" + entries
  }

  /**
    * Method to find the pair of nodes in the matrix that are that have the smallest global distance
    *
    * @return 2-tuple ((Node, Node), distance)
    */
  def getMinPair(): ((String, String), Double) = matrix.minBy(_._2)

  /**
    * Method to merge two given clusters (trees) and update labels (trees) and matrix with re-calculated
    * single-linkage distances
    *
    * @param a ID of cluster to merge (order independent)
    * @param b ID of cluster to merge (order independent)
    * @param id ID given for the new merged cluster
    * @return 2-tuple: (Matrix, Clusters)
    */
  def mergeClusters(a: String, b: String, id: Int): (Matrix, Clusters) = {

    /**
      * Function to determine whether a given pair of IDs corresponds to the IDs of the clusters being merged
      *
      * @return Boolean
      */
    def isInvolved: String => Boolean = x => x == a || x == b

    /**
      * Function to obtain the min distance (single-linkage) for a given cluster ID based on the IDs of the clusters
      * being merged
      *
      * @return
      */
    def singleLinkDist: String => Double = node_id => {
      matrix.filter { case ((x, y), dist) =>
        (x == node_id || y == node_id) && (isInvolved(x) || isInvolved(y))
      }.minBy(_._2)._2
    }

    /**
      * Method to update the kmers of clusters that are leafs by ovewriting them with the set difference of leaf kmer
      * set and the parent intersection
      * @return Clusters
      */
    def updateLeafClusters(): Clusters = {
      //function to get kmers depending on whether it is a leaf and a sketch map was provided
      def getKmers: Tree[Kmers] => Kmers = t => {
        if(leaf2sketches.isEmpty || !t.isLeaf()) t.loadKmers()
        else loadRedwoodSketch(leaf2sketches(t.getID())).sketch.keySet
      }
      //get trees for each cluster
      val (at, bt) = (clusters(a), clusters(b))
      //get kmers
      val (ak, bk) = (getKmers(at), getKmers(bt))
      //get node intersection if either node is a leaf
      val parent_kmers = if(at.isLeaf() || bt.isLeaf()) ak.intersect(bk) else empty_kmers
      //update clusters with node a
      val tmp = if(!clusters(a).isLeaf()) clusters else clusters + (a -> Leaf(ak.diff(parent_kmers), a))
      //update with node b
      if(!clusters(b).isLeaf()) tmp else tmp + (b -> Leaf(bk.diff(parent_kmers),b))
    }

    //create new array of labels by removing old nodes and appending new cluster to front of array wit updated index
    val new_clusters = {
      //update clusters with leafs, if needed
      val updated = updateLeafClusters
      //remove old clusters
      val filtered = updated.filterNot(x => x._1 == a || x._1 == b)
      //create new node cluster
      val new_cluster: Tree[Kmers] = {
        //set left and right sub-trees
        val (l, r) = if (a < b) (updated(a), updated(b)) else (updated(b), updated(a))
        //set parent kmers
        val intersection = clusters(a).loadKmers().intersect(clusters(b).loadKmers())
        //create new node
        Node(id, intersection, fetchDist(a,b), l, r)
      }
      filtered + (id.toString -> new_cluster)
    }

    //create new matrix by re-calculating distances using single-linkage
    val new_matrix = {
      //remove old matrix entry of the clusters being merged
      matrix.filterNot(x => isInvolved(x._1._1) && isInvolved(x._1._2))
        //update single-linkage distance
        .map { case ((x, y), dist) => {
        //current distance does not involve clusters being merged, leave as is
        if (!(isInvolved(x) || isInvolved(y))) ((x, y), dist)
        //update with single linkage
        else {
          val nid = id.toString
          //set new ids for matrix entry get id of the cluster not part of the ones being merged
          val (new_ids, not_merged_id) = if (isInvolved(x)) ((y, nid), y) else ((x, nid), x)
          (new_ids, singleLinkDist(not_merged_id))
        }
      }
      }
    }
    //return new matrix and clusters
    (new_matrix, new_clusters)
  }

}

