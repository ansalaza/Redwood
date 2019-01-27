package sketch

import java.io.File

import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory, verifyFile}

import scala.collection.parallel.ForkJoinTaskSupport

/**
  * Author: Alex N. Salazar
  * Created on 26-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object BuildSketchSerially {

  case class Config(
                     libraries: File = null,
                     outputDir: File = null,
                     kmerSize: Int = 21,
                     sketchSize: Int = 100000,
                     minCov: Int = 2,
                     isAssembly: Boolean = false,
                     threads: Int = 1,
                     trackCov: Boolean = false)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("sketch-serially") {
      opt[File]('l', "libraries") required() action { (x, c) =>
        c.copy(libraries = x)
      } text ("Tab-delimited file containing sample name, assembly, forward reads, and/or reverse reads (if paired " +
        "end). One entry per line.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      note("\nOPTIONAL\n")
      opt[Int]("kmer-size") action { (x, c) =>
        c.copy(kmerSize = x)
      } text ("Kmer size (default is 21).")
      opt[Int]("sketch-size") action { (x, c) =>
        c.copy(sketchSize = x)
      } text ("Sketch size (default is 100000).")
      opt[Int]("min-cov") action { (x, c) =>
        c.copy(minCov = x)
      } text ("Minimum coverage of kmer before adding to sketch (default is 2).")
      opt[Unit]("assembly") action { (x, c) =>
        c.copy(isAssembly = true)
      } text ("Input files are genome assemblies.")
      opt[Unit]("track-cov") action {(x,c) =>
        c.copy(trackCov = true)
      } text("Track coverage of kmers in sketch")
      opt[Int]("threads") action { (x, c) =>
        c.copy(threads = x)
      } text ("Max threads to use (default is 1).")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.libraries)
      constructSketchesSerially(config)
    }

}

  def constructSketchesSerially(config: Config): Unit = {
    //open libraries
    val libraries = openFileWithIterator(config.libraries).toList.map(x => {
      val columns = x.split("\t")
      //create 2-tuple: (sample ID, sequence of sequence files)
      val (id, sequences) = (columns.head, columns.drop(1).map(new File(_)).toSeq)
      //verify sequence files
      sequences.foreach(file => verifyFile(file))
      //return tuple
      (id, sequences)
    }).par
    //set max threads
    libraries.tasksupport = new ForkJoinTaskSupport(new scala.concurrent.forkjoin.ForkJoinPool(config.threads))
    println(timeStamp + "Found " + libraries.size + " samples")
    println(timeStamp + "Sketching:")
    //sketch libraries in parallel
    libraries.foreach{ case (id, seq_files) => {
      //set local output dir
      val local_dir = new File(config.outputDir + "/" + id)
      //create local dir
      local_dir.mkdir()
      //set local config
      val local_config = {
        new BuildSketch.Config(seq_files, id, local_dir,
          config.kmerSize, config.minCov, config.sketchSize,
          config.trackCov, config.isAssembly)
      }
      //build sketch
      BuildSketch.buildSketch(local_config)
    }}
    println(timeStamp + "Successfully completed!")
  }
}
