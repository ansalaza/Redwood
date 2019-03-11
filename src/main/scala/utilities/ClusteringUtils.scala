/**
package utilities

import java.io.File

import utilities.FileHandling.{timeStamp}
import utilities.SketchUtils.loadMinimumSketchSize
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
    * Method to create an initial dendogram given a matrix and column identifiers. Optionally, universal kmer-length, a
    * map of sketch name -> sketch file, and verbosity level.
    *
    * Note: the 'l' char is added to leaf ID's to ensure uniqueness between assigned ID and leaf names when leaf
    * names are integers.
    * @param matrixFile
    * @param kmer_length
    * @param sketches_map
    * @param verbose
    * @return Dendogram
    */
  def createDendogram(matrixFile: File,
                      kmer_length: Int = -1,
                      sketches_map: Map[String, File] = Map(),
                      verbose: Boolean = false): Dendogram = {
    //load matrix
    val (matrix, columns) = loadMatrix(matrixFile)
    println(timeStamp + "--Loaded matrix with " + matrix.size + " values")
    //set minimum sketch size to use
    val min_sketch_size = if(sketches_map.isEmpty) -1 else loadMinimumSketchSize(sketches_map.values.toList)
    //create initial cluster map where key is node/leaf ID -> Tree[Kmers] (initial just the leaf)
    val leafs = columns.map(id => ("l" + id, new Leaf(empty_kmers, "l" + id, -1, min_sketch_size))).toMap
    if (verbose) leafs.foreach(x => println(timeStamp + "--Leafs: " + (x._2.loadAsLeaf().id, x._2.loadKmers().size)))
    //return DistanceMatrix object
    new Dendogram(
      matrix.map(x => (("l" + x._1._1, "l" + x._1._2),x._2)),
        leafs,
      sketches_map.map(x => ("l" + x._1, x._2)),
      min_sketch_size,
      kmer_length,
      verbose)
  }

  /**
    * Method to perform hierarchical clustering for given dendogram, initial node ID identifier. Optional, verbosity
    * level. Outputs a reduced-kmer tree
    * @param dendogram
    * @param acc_node_id
    * @param verbose
    * @return ReducedKmerTree
    */
  def hierchicalClustering(dendogram: Dendogram,
                           acc_node_id: Int,
                           verbose: Boolean = false): ReducedKmerTree = {
    //log
    if (verbose) println(timeStamp + "Clusters remaining: " + dendogram.clusters.size)
    //final two clusters, merge manually
      if(dendogram.clusters.size == 1) {
        println(timeStamp + "--Constructing reduced kmer-tree")
        //construct map of old node ID to new node ID
        val old2new = dendogram.clusters.head._2.inOrderTraversal().zipWithIndex.toMap
        //create reduced tree along with kmer map
        val (rt, km) = dendogram.clusters.head._2.copy2ReducedKmerTree[Kmers](old2new = old2new)
        ReducedKmerTree(rt, km, dendogram.min_sketch_size, dendogram.kmer_length)
      }
    //more clusters to merge
    else {
      if (verbose) println(timeStamp + "--Finding closest pair")
        if(verbose) dendogram.matrix.toList.sortBy(_._2).foreach(println)
      //obtain pair of trees (clusters) with smallest distance
      val closest_pair = dendogram.getMinPair()._1
      //log
      if (verbose) println(timeStamp + "--Merging clusters: " + closest_pair)
      //merge clusters and update matrix and labels with new distances
      val (new_matrix, new_clusters) = dendogram.mergeClusters(closest_pair._1, closest_pair._2, acc_node_id)
      //continue clustering
      hierchicalClustering(
        new Dendogram(new_matrix, new_clusters, dendogram.leaf2sketches,
          dendogram.min_sketch_size,dendogram.kmer_length, verbose), acc_node_id + 1, verbose)
    }
  }

}
  */
