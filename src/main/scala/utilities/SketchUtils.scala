package utilities

import java.io.File

import boopickle.Default._
import utilities.FileHandling._
import utilities.SequenceUtils.ByteEncoded

import scala.annotation.tailrec

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
  case class RedwoodSketch(name: String, kmer_length: Int, genome_size: Int, sketch: Map[Int, (Int, ByteEncoded)])

  /**
    * Function to load a sketch file into a RedwoodSketch object
    *
    * @return RedwoodSketch
    */
  def loadRedwoodSketch: File => RedwoodSketch = file => Unpickle[RedwoodSketch].fromBytes(deserialize(file))

  /**
    * Method to obtain a 3-tuple: (map of sketch id -> sketch file path, kmer length, min sketch size)
    * @param sketches List of sketch files
    * @return 3-tuple as described above
    */
  def loadSketches(sketches: List[File]): (Map[String,File], Int, Int) = {
    @tailrec def _loadSketches(remaining: List[File],
                               smap: Map[String,File],
                               kmer_lengths: Set[Int],
                               sketch_sizes: Set[Int]
                              ): (Map[String, File], Int, Int) = {
      remaining match {
        case Nil => {
          //assert single kmer length
          assert(kmer_lengths.size == 1, "Found different kmer lengths in sketches: " + kmer_lengths.mkString(","))
          if(sketch_sizes.size > 1)
            println(timeStamp + "WARNING: different sketch sizes found: " + sketch_sizes.mkString(","))
          (smap, kmer_lengths.head, sketch_sizes.min)
        }
        case (head::tail) => {
          //get sketch identifier
          val sketch = loadRedwoodSketch(head)
          //assert uniqueness
          assert(!smap.contains(sketch.name), "Dual instance of sketch identifier: " + sketch.name)
          //add to set
          _loadSketches(tail, smap + (sketch.name -> head),
            kmer_lengths + sketch.kmer_length, sketch_sizes + sketch.sketch.size)
        }
      }
    }
    _loadSketches(sketches, Map(), Set(), Set())
  }



  /**
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
    //get universal kmer length
    val kmer_length = isKmerLengthCompatible(_sketches)
    //verify sketch name uniqueness
    isUniqueNameCompatible(_sketches)
    (sketches.map(x => (x.name, x.sketch.keySet)).toMap, kmer_length)
  }

  /**
    * Function to compute verify that all given sketches have the same kmer-length. If true, returns universal
    * kmer-length
    * @return Int
    */
  def isKmerLengthCompatible: List[File] => Int = sketches => {
    //construct list of 2-tuples: (kmer-length, list of sample names)
    val kmer_lengths = sketches.map(x => {val sketch = loadRedwoodSketch(x); (sketch.name, sketch.kmer_length)})
      .groupBy(_._2).mapValues(_.map(_._1)).toList.sortBy(_._2.size)
    //verify all sketches have the same kmer length
    assert(kmer_lengths.size == 1, "One or more sketches have different kmer lengths:\n" +
      kmer_lengths.drop(1).map(x => "--kmer length of " + x._1 + " for samples: " + x._2.mkString(",")))
    //return universal kmer length
    kmer_lengths.head._1
  }

  /**
    * Function to compute and verify that all given sketches have unique name. If true, return list of all sketch names
    * @return List[String]
    */
  def isUniqueNameCompatible: List[File] => List[String] = sketches => {
    val all_names =
      sketches.map(x => {val sketch = loadRedwoodSketch(x); sketch.name}).groupBy(identity).mapValues(_.size)
    //sanity check
    assert(all_names.forall(_._2 == 1), "The following sketch names are not unique: " +
      all_names.filter(_._2 > 1).map(_._1).mkString(","))
    all_names.keys.toList
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
    * Function to obtain minimum sketch size given a list of sketch files
    *
    * @return Int
    */
  def loadMinimumSketchSize: List[File] => Int = sketches => {
    //construct map of sketch size -> List[sketch names]
    val size2names = sketches.map(x => {
      val sketch = loadRedwoodSketch(x)
      (sketch.name, sketch.sketch.size)
    }).groupBy(_._2).mapValues(_.map(_._1))
    //all sketches have the same size
    if (size2names.size == 1) size2names.head._1
    //different sketches have different sizes
    else {
      //get min sketch size
      val min_sketch_size = size2names.minBy(_._1)
      println(timeStamp + "--Input sketches contain different sizes: " + size2names.map(_._1).mkString(","))
      println(timeStamp + "----Using sketch size of " + min_sketch_size._1 + " where appropriate")
      //return min sketch size
      min_sketch_size._1
    }
  }
    */
}
