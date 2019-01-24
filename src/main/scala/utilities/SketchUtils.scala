package utilities

import java.io.File
import java.nio.file.{Files, Paths}

import utilities.FileHandling._
import utilities.SequenceUtils.ByteEncoded

/**
  * Author: Alex N. Salazar
  * Created on 7-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object SketchUtils {

  /**
    * Case class for storing a sketch as an object instance
    *
    * @param name        Name of sample
    * @param kmer_length Kmer size
    * @param sketch      Kmer set
    */
  case class RedwoodSketch(name: String, kmer_length: Int, sketch: Map[Int, (Int, ByteEncoded)])

  /**
    * Function to load a sketch file into a RedwoodSketch object
    *
    * @return RedwoodSketch
    */
  def loadRedwoodSketch: File => RedwoodSketch = file => {
    //deserialize as instance of RedwoodSketch
    deserialize(Files.readAllBytes(Paths.get(file.getAbsolutePath))).asInstanceOf[RedwoodSketch]
  }

  /**
    * Load and verify sketches (existence, uniqueness, kmer-length compatability) for given list of files
    *
    * @param _sketches
    * @return 2-tuple: (Map[String, Set[Int], Kmer length)
    */
  def loadVerifySketches(_sketches: List[File]): (Map[String, Set[Int]], Int) = {
    //verify file exists
    _sketches.foreach(verifyFile(_))
    //load sketches
    val sketches = _sketches.map(loadRedwoodSketch(_))
    //group by kmer length
    val kmer_lengths = sketches.map(x => (x.name, x.kmer_length)).groupBy(_._2).mapValues(_.map(_._1)).toList.sortBy(_._2.size)
    //sanity check
    assert(kmer_lengths.size == 1, "One more sketches have different kmer lengths:\n" +
      kmer_lengths.drop(1).map(x => "--kmer length of " + x._1 + " for samples: " + x._2.mkString(",")))
    //group by name, get counts for each name
    val all_names = sketches.map(_.name).groupBy(identity).mapValues(_.size)
    //sanity check
    assert(all_names.forall(_._2 == 1), "The following sketch names are not unique: " +
      all_names.filter(_._2 > 1).map(_._1).mkString(","))
    (sketches.map(x => (x.name, x.sketch.keySet)).toMap, kmer_lengths.head._1)
  }

  /**
    * Function to construct sketches map given a list of redwoodsketch files
    *
    * @return Map[String, File]
    */
  def constructSketchesMap: File => Map[String, File] = file => {
    //no sketch file provide, return empty map
    if (file == null) {
      println(timeStamp + "--Sketches file was not provided")
      Map[String, File]()
    }
    //sketch file provided
    else {
      println(timeStamp + "--Sketches file was provided")
      //obtain list of ids and sketch paths
      val tmp = openFileWithIterator(file).toList.foldLeft(List[(String, File)]())((list, line) => {
        //get file path
        val sketch_file = new File(line)
        //get id
        val id = getFileName(sketch_file)
        //verify sketch path exist
        verifyFile(sketch_file)
        //append
        (id, sketch_file) :: list
      })
      //verify all ids are unique
      tmp.foldLeft(Set[String]()) { case (set, (id, sketch)) => {
        assert(!set(id), "Non-unique ID found: " + id)
        set + (id)
      }
      }
      println(timeStamp + "--Found " + tmp.size + " sketches")
      //return map
      tmp.toMap
    }
  }

  /**
    * Method to obtain all kmers unique to a given list of samples. First constructs the set intersection between the
    * given list of samples. Then, iteratively subtracts kmers from the constructed intersection using kmers from all
    * other samples.
    *
    * @param names        Name of samples
    * @param sketches_map Map of sample names -> sketch file
    * @return Kmer-set intersection as Set[Long]
    */
  def uniqueKmers(names: List[String], sketches_map: Map[String, File]): Set[Int] = {
    /**
      * Function to corresponding kmer set of given sample/key
      *
      * @return Set[Long]
      */
    def fetch: String => Set[Int] = key => {
      //attempt to get sketch file
      val sketch = sketches_map.get(key)
      assert(sketch.nonEmpty, "Could find corresponding sketch file for sample: " + key)
      //load sketch and extract kmer hashes
      loadRedwoodSketch(sketch.get).sketch.keySet
    }

    //iterate through all given samples and construct intersecting set
    val intersect = {
      //obtain initial kmer set
      val initial_set = fetch(names.head)
      names.tail.foldLeft(initial_set)((acc_set, name) => acc_set.union(fetch(name)))
    }
    //set of all intersecting names
    val names_set = names.toSet
    //iteratively perform set difference between the sample intersect above and all other samples
    sketches_map.toList.filter(x => !names_set(x._1))
      .foldLeft(intersect) { case (acc_intersect, (name, sketch)) => acc_intersect.diff(fetch(name)) }
  }
}
