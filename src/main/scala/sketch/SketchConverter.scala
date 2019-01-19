package sketch

import java.io.{BufferedOutputStream, File, FileOutputStream}

import utilities.FileHandling._

import scala.collection.parallel.ForkJoinTaskSupport

/**
  * Author: Alex N. Salazar
  * Created on 7-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object SketchConverter {

  case class Config(
                     sketch: File = null,
                     outputDir: File = null,
                     multiSketch: Boolean = false,
                     threads: Int = 1,
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("sketch-converter") {
      opt[File]('s', "sketch-file") required() action { (x, c) =>
        c.copy(sketch = x)
      } text ("Sketch file (MASH format.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      note("\nOPTIONAL\n")
      opt[Unit]("multi-sketch") action { (x, c) =>
        c.copy(multiSketch = true)
      } text ("Input for '--sketch-file' is a file containing paths to multiple sketches, one per line.")
      opt[Int]('t', "threads") action { (x, c) =>
        c.copy(threads = x)
      } text ("Max number of threads..")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.sketch)
      assert(config.threads > 0, "Number of threads specified is a non-positive integer")
      if (config.multiSketch) openFileWithIterator(config.sketch).foreach(file => verifyFile(new File(file)))
      sketchConverter(config)
    }
  }


  def sketchConverter(config: Config): Unit = {
    //set number of threads
    val forkJoinPool = new java.util.concurrent.ForkJoinPool(config.threads)
    //load single sketch or multiple sketches into a list
    val sketches = {
      //convert single sketch file
      if (!config.multiSketch) {
        println(timeStamp + "Loading single sketch")
        List(config.sketch)
      }
      //convert multiple sketches
      else {
        println(timeStamp + "Loading multiple sketches")
        openFileWithIterator(config.sketch).toList.map(new File(_))
      }
    }.par
    //specify threads in collection
    sketches.tasksupport = new ForkJoinTaskSupport(forkJoinPool)
    //log output
    if (config.multiSketch) println(timeStamp + "--Loaded " + sketches.size + " sketch files")
    println(timeStamp + "Converting:")
    //iterate through each sketch and convert
    /**
    sketches.foreach(sketch => {
      //convert sketch
      val rdw_sketch = mash2Rdw(sketch)
      println(timeStamp + "--Converted sketch for " + rdw_sketch.name + " containing " + rdw_sketch.kmers.size +
        " kmers of size " + rdw_sketch.kmerSize)
      //create output file
      val pw = new BufferedOutputStream(new FileOutputStream(config.outputDir + "/" + rdw_sketch.name + ".rdws"))
      //serialize and write to disk
      pw.write(serialise(rdw_sketch))
      pw.close
    })
      */

    println(timeStamp + "Successfully completed!")
  }

}
