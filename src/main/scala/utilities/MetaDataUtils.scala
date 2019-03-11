package utilities

import java.io.File

import doodle.core.Color
import doodle.core.Color.hsl
import utilities.FileHandling.openFileWithIterator
import utilities.ReducedKmerTreeUtils.ReducedTree
import doodle.syntax._

/**
  * Author: Alex N. Salazar
  * Created on 30-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object MetaDataUtils {

  /**
    * Function to load a colours file provided colours file
    *
    * @return Map[String, HSL]
    */
  def loadColoursFile: File => Map[String, Color] = file => {
    //no colours file, return empty map
    if (file == null) Map.empty[String, Color]
    //open
    else openFileWithIterator(file).toList.map(x => {
      //get columns
      val columns = x.split("\t")
      //set colour
      (columns.head,
        hsl(columns(1).toDouble.degrees, columns(2).toDouble.normalized, columns(3).toDouble.normalized))
    }).toMap
  }

  /**
    * Method to load a query file
    * @param file
    * @param max_width
    * @return Map[Int,Double]
    */
  def loadQueryFile(file: File, max_width: Double): Map[Int, Double] = {
    openFileWithIterator(file).drop(1).toList.foldLeft(Map[Int, Double]())((acc, line) => {
      val columns = line.split("\t")
      acc + (columns.head.toInt -> (columns(2).toDouble * max_width))
    })
  }

  /**
    * Method to load a query label file
    * @param file
    * @return
    */
  def loadQueryLabelFile(file: File): (String, Map[String, Double]) = {
    val m = openFileWithIterator(file).drop(1).toList.foldLeft((Map[String, Double](), Set[String]()))((acc, line) => {
      val columns = line.split("\t")
      (acc._1 + (columns(1) -> (columns(2).toDouble)), acc._2 + columns.head)
    })
    assert(m._2.size == 1, "Expected identical names but found different ones: " + m._2.mkString(","))
    (m._2.head, m._1)
  }

}
