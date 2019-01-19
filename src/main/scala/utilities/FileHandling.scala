package utilities

import java.io._
import java.util.Calendar
import java.util.zip.GZIPInputStream

import utilities.SketchUtils.RedwoodSketch

import scala.io.Source

/**
  * Author: Alex N. Salazar
  * Created on 3-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object FileHandling {

  /**Method to get current time stamp*/
  def timeStamp = "(" + Calendar.getInstance().getTime() + "): "

  /**Method to get name of file*/
  def getFileName: File => String = file => file.getName.substring(0, file.getName.lastIndexOf("."))

  /**Method to get extension of file*/
  def getFileExtension: File => String = file => file.getName.substring(file.getName.lastIndexOf(".")+1)

  /**Method to open file with Iterator*/
  def openFileWithIterator(x: File, remove_header_comment_lines: Boolean = false): Iterator[String] = {
    if(remove_header_comment_lines) Source.fromFile(x).getLines.dropWhile(_.startsWith("#"))
    else Source.fromFile(x).getLines
  }

  /**Method to open gziped file with Iterator*/
  def openGzipFileWithIterator(x: File): Iterator[String] = {
    Source.fromInputStream(new GZIPInputStream(new BufferedInputStream(new FileInputStream(x)))).getLines
  }

  /**Check whether directory exists and if it is valid. If it does not exist, create it.*/
  def checkOutputDirectory(inputDir: File) {
    if (inputDir.exists()) {
      assert(inputDir.isDirectory(), "The following directory does not exist or it is invalid: " + inputDir.getAbsolutePath())
    } else {
      inputDir.mkdir()
      assume(inputDir.exists(), "Could not create the following directory: " + inputDir.getAbsolutePath())
    }
  }

  /**Check whether directory exists and if it is valid. If it does not exist, create it.*/
  def verifyDirectory(inputDir: File, message: String = "The following directory does not exist or it is invalid") {
    assert(inputDir.exists() && inputDir.isDirectory(), message + ": " + inputDir.getAbsolutePath())
  }

  /**Check whether directory exists and if it is valid. If it does not exist, create it.*/
  def verifyFile(inputFile: File, message: String = "The following file does not exist or it is invalid") {
    assert(inputFile.exists() && inputFile.isFile(), message + ": " + inputFile.getAbsolutePath())
  }

  def outputSketch(sketch: RedwoodSketch, outputFile: File): Unit = {
    //create output file
    val pw = new BufferedOutputStream(new FileOutputStream(outputFile))
    pw.write(serialise(sketch))
    pw.flush()
    pw.close
  }

  def serialise: Any => Array[Byte] = some_object => {
    val stream: ByteArrayOutputStream = new ByteArrayOutputStream()
    val oos = new ObjectOutputStream(stream)
    oos.writeObject(some_object)
    oos.reset()
    oos.close()
    stream.toByteArray
  }


  def deserialise: Array[Byte] => Any = some_byte_array => {
    val ois = new ObjectInputStream(new ByteArrayInputStream(some_byte_array))
    val value = ois.readObject
    ois.close()
    value
  }

}
