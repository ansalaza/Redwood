package utilities

import java.io.File

import utilities.FileHandling.{openFileWithIterator, timeStamp}
import utilities.SketchUtils.{constructSketchesMap, loadRedwoodSketch}
import utilities.KmerTreeUtils.{Kmers, Leaf, Node, Tree}

/**
  * Author: Alex N. Salazar
  * Created on 4-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object ClusteringUtils {

  type Matrix = Array[(Array[(Double, Int)], Int)]
  type Clusters = Array[(Tree[Kmers], Int)]
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
    //load matrix as array
    val matrix = openFileWithIterator(matrix_file).toArray.map(_.split("\t"))
    //load sketches, if they exist
    val sketches_map = constructSketchesMap(sketches_file)
    //get labels and their indeces
    val leafs: Array[(Tree[Kmers], Int)] = {
      //get all columns (leafs)
      val all_leafs = matrix.head.tail.zipWithIndex
      //log
      if (sketches_file != null)
        println(timeStamp + "--Loading kmers for " + all_leafs.size + " leafs")
      //iterate through each leaf and construct leaf object (with unique kmers if specified)
      all_leafs.map(x => {
        //sketches provided, move on
        if (sketches_map.isEmpty) (new Leaf(empty_kmers, x._1), x._2)
        else {
          //fetch sketch file
          val sketch = sketches_map.get(x._1)
          //sanity check
          assert(sketch.nonEmpty, "Could not find sketch file for sample: " + x._1)
          //obtain kmers of current leaf
          (new Leaf(loadRedwoodSketch(sketch.get).sketch.keySet, x._1), x._2)
        }
      })
    }
    if (verbose) {
      leafs.foreach(x => println(timeStamp + "--Leafs: " + (x._1.loadAsLeaf().id, x._1.loadKmers().size)))
    }
    //get distance values as array of doubles
    val distances = matrix.drop(1).map(_.tail.map(_.toDouble).zipWithIndex).zipWithIndex
    //create map of all pairwise distances
    val pairwise_map = {
      //create mapt of index to id
      val idx2map = leafs.map(x => (x._2, x._1.asInstanceOf[Leaf].id)).toMap
      //iterate through each row and column and append pairwise distance
      distances.foldLeft(List[((String, String), Double)]()) { case (map, (row, i)) => {
        row.foldLeft(map) { case (acc_map, (dist, j)) => ((idx2map(i), idx2map(j)), dist) :: acc_map }
      }
      }.toMap
    }
    //return DistanceMatrix object
    new Dendogram(distances, leafs, pairwise_map, sketches_map)
  }

  /**
    * Method to perform hierchical clustering given a DistanceMatrix object using single-linkage distance metric
    *
    * @param matrix DistanceMarix object
    * @return DistanceMatrix
    */
  def hierchicalClustering(matrix: Dendogram,
                           acc_node_id: Int,
                           sketches_map: Map[String, File] = Map(),
                           verbose: Boolean = false): Node = {

    /**
      * Function to load kmer set from a given tree. For leaf nodes, this is all kmers unique to the leaf sample. For
      * tree nodes, its the union of left and right sub-tree. If a sub-tree is a node, simply loads it's kmer values
      * as it assumes it has already constructed kmers unique to all children.
      *
      * @return Kmers
      */
    def loadKmersFromTree: Tree[Kmers] => Kmers = tree => {

      //no sketches provided, return empty
      if (sketches_map.isEmpty) empty_kmers
      else {
        //tree is a leaf, obtain all unique kmer to current leaf sample
        if (tree.isLeaf()) tree.loadKmers()
        else {
          //load node
          val node = tree.loadAsNode()
          //obtain union of left and right
          node.left.loadKmers().union(node.right.loadKmers())
        }
      }
    }

    if (verbose) println("Clusters remaining: " + matrix.clusters.size)
    //final two clusters, merge manually
    if (matrix.clusters.size == 2) {
      //get min dist as sum dist
      val sum_dist = matrix.getMinPair()._1
      //get left and right trees based on lexicographic order
      val (l, r) = {
        val tmp = matrix.clusters.sortBy(x => x._1.getID())
        (tmp.head._1, tmp(1)._1)
      }
      //set kmer set
      val kmer_set = l.loadKmers().union(r.loadKmers())
      //create tree root as node
      new Node(acc_node_id, kmer_set, sum_dist, l, r)
    }
    else {
      if (verbose) println("--Finding closest pair")
      //obtain pair of trees (clusters) with smallest distance
      val (min_dist, closest_pair) = matrix.getMinPair()
      //log
      if (verbose) println("--Computing kmer set union")
      //get set union of kmers, if sketch map provided
      val kmer_set = loadKmersFromTree(closest_pair._1._1).union(loadKmersFromTree(closest_pair._2._1))
      if (verbose) println("--Set union of " + kmer_set.size)
      //log
      if (verbose) println("--Merging clusters")
      //merge clusters and update matrix and labels with new distances
      val (new_matrix, new_labels) = matrix.mergeClusters(closest_pair._1, closest_pair._2, min_dist, kmer_set, acc_node_id)
      if (verbose) {
        println("--New clusters: \n" + new_labels.map(x => "----" + x._1.getLeafNames().sorted.mkString("\n")))
      }
      //continue clustering
      hierchicalClustering(
        new Dendogram(new_matrix, new_labels, matrix.pairwise_leaf_dist),
        acc_node_id + 1,
        sketches_map,
        verbose)
    }
  }

}
