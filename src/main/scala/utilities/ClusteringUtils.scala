package utilities

import java.io.File

import utilities.FileHandling.{timeStamp}
import utilities.SketchUtils.{constructSketchesMap, loadRedwoodSketch}
import utilities.KmerTreeUtils.{Kmers, Leaf, Node, Tree}
import utilities.DistanceUtils.loadMatrix

/**
  * Author: Alex N. Salazar
  * Created on 4-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object ClusteringUtils {

  type Matrix = Map[(String, String), Double]
  type Clusters = Map[String, Tree[Kmers]]
  val empty_kmers: Kmers = Set()

  /**
    * Method to create an DistanceMatrix object given a matrix file object
    *
    * @param matrix_file
    * @param sketches_file
    * @param verbose
    * @return
    */
  def createDendogram(matrix_file: File,
                      sketches_file: File = null,
                      verbose: Boolean = false): Dendogram = {
    //load matrix
    val (matrix, all_leafs) = loadMatrix(matrix_file)
    //load sketches, if they exist
    val sketches_map = constructSketchesMap(sketches_file)
    //create initial cluster map where key is node/leaf ID -> Tree[Kmers] (initial just the leaf)
    val leafs: Map[String, Tree[Kmers]] = {
      //sketch file not provided
      if (sketches_file == null) all_leafs.map(x => (x, new Leaf(empty_kmers, x))).toMap
      else {
        /**
          * Function to fetch the sketch for a given ID
          * @return Kmers
          */
        def fetchSketch: String => Kmers = id => {
          //fetch sketch file
          val sketch = sketches_map.get(id)
          //sanity check
          assert(sketch.nonEmpty, "Could not find sketch file for sample: " + id)
          //load sketch
          loadRedwoodSketch(sketch.get).sketch.keySet
        }
        println(timeStamp + "--Loading kmers for " + all_leafs.size + " leafs")
        //iterate through each leaf and construct leaf object (with unique kmers if specified)
        all_leafs.map(x => (x, new Leaf(fetchSketch(x), x))).toMap
      }
    }
    if (verbose) leafs.foreach(x => println(timeStamp + "--Leafs: " + (x._2.loadAsLeaf().id, x._2.loadKmers().size)))
    //return DistanceMatrix object
    new Dendogram(matrix, leafs)
  }

  /**
    * Method to perform hierchical clustering given a DistanceMatrix object using single-linkage distance metric
    *
    * @param matrix DistanceMarix object
    * @return DistanceMatrix
    */
  def hierchicalClustering(dendogram: Dendogram, acc_node_id: Int, verbose: Boolean = false): Node = {
    //log
    if (verbose) println("Clusters remaining: " + dendogram.clusters.size)
    //final two clusters, merge manually
    if (dendogram.clusters.size == 2) {
      //get min dist as sum dist
      val sum_dist = dendogram.matrix(dendogram.getMinPair()._1)
      //get left and right trees based on lexicographic order
      val (l, r) = {
        val tmp = dendogram.clusters.toList.sortBy(_._1).map(_._2)
        (tmp.head, tmp(1))
      }
      //set kmer set
      val kmer_set = l.loadKmers().union(r.loadKmers())
      //create tree root as node
      new Node(acc_node_id, kmer_set, sum_dist, l, r)
    }
    //more clusters to merge
    else {
      if (verbose) println("--Finding closest pair")
      //obtain pair of trees (clusters) with smallest distance
      val closest_pair = dendogram.getMinPair()._1
      //log
      if (verbose) println("--Merging clusters")
      //merge clusters and update matrix and labels with new distances
      val (new_matrix, new_clusters) = dendogram.mergeClusters(closest_pair._1, closest_pair._2, acc_node_id)
      //continue clustering
      hierchicalClustering(new Dendogram(new_matrix, new_clusters), acc_node_id + 1, verbose)
    }
  }

}
