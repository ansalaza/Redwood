package utilities

import java.io.File
import java.nio.file.{Files, Paths}

import utilities.FileHandling._

import scala.sys.process.Process

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
    * @param name     Name of sample
    * @param kmerSize Kmer size
    * @param kmers    Kmer set
    */
  case class RedwoodSketch(name: String, kmerSize: Int, kmers: Map[Int, Double])

  /**
    * Function to load a sketch file into a RedwoodSketch object
    *
    * @return RedwoodSketch
    */
  def loadRedwoodSketch: File => RedwoodSketch = file => {
    //deserialize as instance of RedwoodSketch
    deserialise(Files.readAllBytes(Paths.get(file.getAbsolutePath))).asInstanceOf[RedwoodSketch]
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
      loadRedwoodSketch(sketch.get).kmers.keySet
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
      .foldLeft(intersect) { case (acc_intersect, (name, sketch)) => acc_intersect.diff(fetch(name))}
  }

  /**
    * Function to load a mash sketch into a RedwoodSketch object
    *
    * @return RedwoodSketch
    *
  def mash2Rdw: File => RedwoodSketch = file => {
    //get name
    val name = getFileName(file)
    //dump sketch using mash
    val line = Process(Seq("mash", "info", "-d", file.getAbsolutePath)).!!
    //get kmer size
    val k_size = line.split("\n").find(_.replaceAll("\\s+", "").startsWith("\"kmer\":")).get.filter(_.isDigit).toInt
    //construct hash set
    val k_set = {
      //remove all ines up to start of hashes
      line.split("\n").dropWhile(x => !x.replaceAll("\\s+", "").startsWith("\"hashes\":")).drop(2)
        //iterate through each line and parse as hash (INT), append to list
        .foldLeft(List[Long]())((kmers, line) => {
        val filtered = line.filter(_.isDigit)
        if (filtered.isEmpty) kmers else filtered.toLong :: kmers
      }).toSet
    }
    new RedwoodSketch(name, k_size, k_set)
  }
  */
}
