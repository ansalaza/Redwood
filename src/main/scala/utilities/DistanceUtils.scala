package utilities

import java.io.File

import utilities.ClusteringUtils.Matrix
import utilities.FileHandling.{openFileWithIterator, timeStamp}
import utilities.NumericalUtils.min

import scala.math.log
import scala.collection.mutable

/**
  * Author: Alex N. Salazar
  * Created on 24-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object DistanceUtils {

  /**
    * Method to compute MASH distance given two IDs
    *
    * @return Double
    */
  def computeMashDist: (Set[Int], Set[Int], Int) => Double = (x, y, kmer_length) => {
    //compute union
    val union = broderUnion(x, y)
    //compute intersection
    val intersect = union.filter(kmer => x(kmer) && y(kmer)).size.toDouble
    //compute jaccard distance
    val ji = (intersect / union.size)
    //compute mash distance
    (-1.0 / kmer_length) * log((2 * ji) / (1 + ji))
  }

  /**
    * Function to get the mash distance for a given subj and target ID accounting that the distance map only
    * includes those above the diagonal
    *
    * @return Double
    */
  def fetchDistance(distance_map: Map[(String, String), Double])(x: String, y: String): Double = {
    //check forward order, if not, swap order
    if (distance_map.contains((x, y))) distance_map((x, y)) else distance_map((y, x))
  }

  /**
    * Method to load a given matrix file into a map consisting of all distances in the upper diagonal.
    * @param file
    * @return 2-tuple: (Map where key is 2-tuple of sample ID -> distance, List of all IDs)
    */
  def loadMatrix(file: File): (Matrix, List[String]) = {
    //get identifiers
    val index2ID = openFileWithIterator(file).next().split("\t").drop(1).zipWithIndex.map(_.swap).toMap
    //total identfiers
    val total_ids = index2ID.size
    println(timeStamp + "Found matrix with " + total_ids + " samples")
    //construct matrix
    val matrix = openFileWithIterator(file).drop(1)
      //iterate through each row as i
      .foldLeft((List[((String, String), Double)](), 0)) { case ((distances, i), line) => {
      //get row id
      val row_id = index2ID(i)
      //iterate through each value as j
      val (updated, total_values) = line.split("\t").drop(1).foldLeft((distances, 0)) { case ((acc,j), dist) => {
        //sanity check
        //only add distances from upper diagonal
        if(j <= i) (acc, j+1) else (((row_id, index2ID(j)), dist.toDouble) :: acc, j + 1)
      }}
      //sanity check
      assert(total_values == total_ids, "Unexpected number of distances row: " + line)
      //return updated distances
      (updated, i + 1)
    }}._1
    (matrix.toMap, index2ID.values.toList)
  }

  /**
    * Method to obtain's Broder's union of two sets. In short, the smallest X hashes from the union of two sets.
    *
    * @param x
    * @param y
    * @return
    */
  private def broderUnion(x: Set[Int], y: Set[Int]): Set[Int] = {
    //get min sketch size to use
    val min_sketch_size = {
      //get size of each sketch
      val (x_size, y_size) = (x.size, y.size)
      //same sketch size
      if (x_size == y_size) x_size
      else {
        println("--Warning: uneven sketch sizes " + (x_size, y_size) + ". Using sketch size of " + min(x_size, y_size))
        //get min sketch size
        min[Int](x_size, y_size)
      }
    }
    //create broder union priority queue
    var union_pq = new mutable.PriorityQueue[Int]()
    //set union size
    var current_union_size = 0

    /**
      * Function to determine whether a given hash is less than current union's max
      *
      * @return
      */
    def lessThanMax: Int => Boolean = hash => current_union_size < min_sketch_size || hash < union_pq.head

    //create hashset for kmers in union
    var isInUnion = new mutable.HashSet[Int]()
    //iterate through each set and add kmers to priority queue
    List(x, y).foreach(set => {
      //iterate through kmer hashses in current set and add to priority queue
      set.foreach(hash => {
        //less than current max and not seen before
        if (lessThanMax(hash) && !isInUnion(hash)) {
          //priority queue not filled yet
          if (current_union_size < min_sketch_size) {
            //add hash
            union_pq.enqueue(hash)
            isInUnion.add(hash)
            current_union_size += 1
          }
          //priority queue filled, remove current max and add hash
          else {
            isInUnion.remove(union_pq.head)
            union_pq.dequeue()
            union_pq.enqueue(hash)
            isInUnion.add(hash)
          }
        }
      })
    })
    //set broder_union
    val broder_union = union_pq.toSet
    assert(broder_union.size == min_sketch_size)
    //return Broder's union
    broder_union
  }


}
