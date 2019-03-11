/**
package utilities

import java.io.File

import utilities.FileHandling.timeStamp
import utilities.ClusteringUtils._
import utilities.KmerTreeUtils.{Kmers, Leaf, Node, Tree}
import utilities.SketchUtils.{loadRedwoodSketch}
import atk.ProgressBar.progress

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
class Dendogram(val matrix: Matrix,
                val clusters: Clusters,
                val leaf2sketches: Map[String, File],
                val min_sketch_size: Int,
                val kmer_length: Int,
                val verbose: Boolean) {

  /**
    * Function to obtain distance of a given pair of IDs aware that the matrix only stores distances in upper-diagonal
    *
    * @return Double
    */
  private def fetchDist: (String, String) => Double = (x, y) =>
    if (x == y) 0.0 else if (matrix.contains((x, y))) matrix((x, y)) else matrix((y, x))

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
    * @param a  ID of cluster to merge (order independent)
    * @param b  ID of cluster to merge (order independent)
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
      * Generate all possible tuple configurations in which a given ID can exist in the matrix map with the either of
      * the cluster IDs being merged. There are 4 different ways:
      * ----(node_id, a), (node_id, b), (a, node_id), (b, node_id)
      * @return List[(String,String)]
      */
    def allPossibleTuples: String => List[(String,String)] = id => List((id, a), (id,b), (a,id),(b,id))

    /**
      * Function to obtain the min distance (single-linkage) between a given cluster ID and the IDs of the clusters
      * being merged
      *
      * @return Double
      */
    def singleLinkDist: String => Double = node_id => {
      //set infinity
      val infinity = Double.MaxValue
      //get min distance
      val min_dist = allPossibleTuples(node_id).map(x => matrix.getOrElse(x, infinity)).min
      assert(min_dist != infinity)
      min_dist
    }

    /**
      * Method to update the kmers of clusters that are leafs by ovewriting them with the set difference of leaf kmer
      * set and the parent intersection
      *
      * @return Clusters
      */
    def updateLeafClusters(): Clusters = {
      /**
        * Function to return a 3-tuple for a node: (Kmers, genome size, sketch size). Note the second to last applies
        * to leafs, otherwise set to -1 value.
        */
      def loadInfo: Tree[Kmers] => (Kmers, Int, Int) = t => {
        if (leaf2sketches.isEmpty || !t.isLeaf()) (t.loadKmers(), -1, -1)
        else {
          //load sketch
          val sketch = loadRedwoodSketch(leaf2sketches(t.getID()))
          //get size of current sketch
          val sketch_size = sketch.sketch.size
          //adjust kmers based on minimum sketch size, if needed
          val kmers = {
            //min sketch size is same as current sketch size
            if (min_sketch_size == sketch_size) sketch.sketch.keySet
            //adjust sketch by getting X smallest kmers (X = min sketch size)
            else {
              println(timeStamp + "-----Adjusting sketch to size " + min_sketch_size + " for " + sketch.name)
              sketch.sketch.keySet.toList.sorted.take(min_sketch_size).toSet
            }
          }
          (kmers, sketch.genome_size, min_sketch_size)
        }
      }

      if (verbose) println(timeStamp + "--Loading clusters")
      //get trees for each cluster
      val (at, bt) = (clusters(a), clusters(b))
      if (verbose) println(timeStamp + "--Loading sketches, if any")
      //load sketches
      val (as, bs) = (loadInfo(at), loadInfo(bt))
      //get node intersection if either node is a leaf
      val parent_kmers = if (at.isLeaf() || bt.isLeaf()) as._1.intersect(bs._1) else empty_kmers
      //update clusters with node a
      val tmp = if (!clusters(a).isLeaf()) clusters else clusters + (a -> Leaf(as._1.diff(parent_kmers), a, as._2, as._3))
      //update with node b
      if (!clusters(b).isLeaf()) tmp else tmp + (b -> Leaf(bs._1.diff(parent_kmers), b, bs._2, bs._3))
    }

    if (verbose) println(timeStamp + "--Updating clusters")
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
        Node(id, intersection, fetchDist(a, b), l, r)
      }
      filtered + (id.toString -> new_cluster)
    }

    if (verbose) println(timeStamp + "--Updating matrix")
    //create new matrix by re-calculating distances using single-linkage
    val new_matrix = {
      //iterate through current clusters IDs using matrix with the distances of clusters being merged removed
      clusters.keys.foldLeft(matrix - ((a,b)))((updated_matrix, cluster_id) => {
        //current cluster is one of the clusters being merged, move on
        if(isInvolved(cluster_id)) updated_matrix
        else {
          //set new cluster ID
          val nid = id.toString
          //get tuple configuration that exists in the matrix map
          val existing_tuple_id = allPossibleTuples(cluster_id).filter(x => matrix.contains(x))
          assert(existing_tuple_id.size == 2, existing_tuple_id)
          //update with single linkage clustering
          existing_tuple_id.foldLeft(updated_matrix)((acc, old) => acc - old) +
            ((nid, cluster_id) -> singleLinkDist(cluster_id))
        }
      })
    }
    //return new matrix and clusters
    (new_matrix, new_clusters)
  }

}

*/