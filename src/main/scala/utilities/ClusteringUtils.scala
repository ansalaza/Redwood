package utilities

import java.io.File

import utilities.FileHandling.{timeStamp}
import utilities.SketchUtils.{constructSketchesMap}
import utilities.KmerTreeUtils.{Kmers, Leaf, Tree}
import utilities.ReducedKmerTreeUtils.ReducedKmerTree
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
    val leafs = all_leafs.map(id => (id, new Leaf(empty_kmers, id))).toMap
    if (verbose) leafs.foreach(x => println(timeStamp + "--Leafs: " + (x._2.loadAsLeaf().id, x._2.loadKmers().size)))
    //return DistanceMatrix object
    new Dendogram(matrix, leafs, sketches_map)
  }

  /**
    * Method to perform hierchical clustering given a DistanceMatrix object using single-linkage distance metric
    *
    * @param matrix DistanceMarix object
    * @return DistanceMatrix
    */
  def hierchicalClustering(dendogram: Dendogram,
                           acc_node_id: Int,
                           verbose: Boolean = false): ReducedKmerTree = {
    //log
    if (verbose) println("Clusters remaining: " + dendogram.clusters.size)
    //final two clusters, merge manually
      if(dendogram.clusters.size == 1) {
        println(timeStamp + "--Constructing reduced kmer-tree")
        //construct map of old node ID to new node ID
        val old2new = dendogram.clusters.head._2.inOrderTraversal().zipWithIndex.toMap
        //create reduced tree along with kmer map
        val (rt, km) = dendogram.clusters.head._2.copy2ReducedKmerTree[Kmers](old2new = old2new)
        ReducedKmerTree(rt, km)
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
      hierchicalClustering(new Dendogram(new_matrix, new_clusters, dendogram.leaf2sketches), acc_node_id + 1, verbose)
    }
  }

}
