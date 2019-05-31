package sketch

/**
  * Author: Alex N. Salazar
  * Created on 9-1-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */

import java.io.{File, PrintWriter}

import utilities.FileHandling.{openFileWithIterator, verifyDirectory, verifyFile}

object BuildSketchUsingCluster {

  case class Config(
                     libraries: File = null,
                     redwoodBinary: File = null,
                     outputDir: File = null,
                     kmerSize: Int = 21,
                     maxMemory: Int = 20000,
                     sketchSize: Int = 100000,
                     minCov: Int = 2,
                     isAssembly: Boolean = false,
                     trackCov: Boolean = false,
                     clusterConfig: File = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("sketch-cluster") {
      opt[File]('l', "libraries") required() action { (x, c) =>
        c.copy(libraries = x)
      } text ("Tab-delimited file containing sample name, assembly, forward reads, and/or reverse reads (if paired " +
        "end). One entry per line.")
      opt[File]("redwood-binary") required() action { (x, c) =>
        c.copy(redwoodBinary = x)
      } text ("Full path to redwood jar file.")
      opt[File]("cluster-config") required() action { (x, c) =>
        c.copy(clusterConfig = x)
      } text ("Scheduler configuration file. Sketeches can be computed in parallel via a cluster scheduler if a " +
        "configuration file is provided with the native scheduler parameters. See README.md for format specification.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      note("\nOPTIONAL\n")
      opt[Int]("kmer-size") action { (x, c) =>
        c.copy(kmerSize = x)
      } text ("Kmer size (default is 21).")
      opt[Int]("min-cov") action { (x, c) =>
        c.copy(minCov = x)
      } text ("Minimum coverage of kmer before adding to sketch (default is 2).")
      opt[Int]("sketch-size") action { (x, c) =>
        c.copy(sketchSize = x)
      } text ("Sketch size (default is 100000).")
      opt[Unit]("assembly") action { (x, c) =>
        c.copy(isAssembly = true)
      } text ("Input files are genome assemblies.")
      opt[Unit]("track-cov") action {(x,c) =>
        c.copy(trackCov = true)
      } text("Track coverage of kmers in sketch")
      opt[Int]("memory") action { (x, c) =>
        c.copy(maxMemory = x)
      } text ("Max memory for JVM in mb (default is 2000).")

    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.libraries)
      verifyFile(config.redwoodBinary)
      constructSketches(config)
    }
  }

  def constructSketches(config: Config): Unit = {
    /**
      * Function to create redwood sketch command given a sequence of read files and the sample name
      * assuming given samples are one or more FASTQ.GZ files
      * @return Redwood sketch command with user-specified parameters
      */
    def makeRedwoodSketchCommand: (Seq[File],String) => String = (files,name) => {
      Seq("java", "-Xmx" + config.maxMemory + "m", "-jar", config.redwoodBinary.getAbsolutePath, "sketch",
        "--seq-files", files.map(_.getAbsolutePath).mkString(","),
        "--kmer-size", config.kmerSize,
        "--sketch-size", config.sketchSize,
        "--min-cov", config.minCov,
        "--name", name,
        "-o", (config.outputDir.getAbsolutePath + "/" + name)
      ).mkString(" ") + (if(!config.isAssembly) "" else " --assembly") + (if(!config.trackCov) "" else " --track-cov")
    }

    //open target file
    val libraries =
      openFileWithIterator(config.libraries).toList.map(_.split("\t")).map(x => (x.head, x.tail.toSeq.map(new File(_))))
    //open configuration file, if provided. If not, return empty list
    val cluster_config = {
      if(config.clusterConfig == null) List("\\$COMMAND")
      else openFileWithIterator(config.clusterConfig).toList
    }
    println("Found " + libraries.size + " samples/libraries to sketch.")
    //iterate through each library and create mash sketch command
    libraries.foreach(library => {
      println("--Processing: " + library._1)
      //check if the files provided are valid, if not, throw warning and move on
      if(!library._2.forall(x => x.exists() && x.isFile))
        println("----WARNING: one of the libraries for sample " + library._1 + " is not a valid file.")
      else {
        //get local directory path
        val local_dir_path = config.outputDir.getAbsolutePath + "/" + library._1
        //create local directory
        val local_dir = new File(local_dir_path)
        local_dir.mkdir()
        //create local slurm script
        val local_pw = new PrintWriter(local_dir_path + "/runRedwoodSketch.sh")
        //create mash command
        val rdw_command =  makeRedwoodSketchCommand(library._2, library._1)
        //output mash command to local script
          cluster_config.foreach(line => {
            //for each line, perform the following
            local_pw.println(
              //add file path to stdout file
              line.replaceAll("\\$STDOUT", local_dir_path + "/" + library._1 + ".out")
                //add mash command
                .replaceAll("\\$COMMAND", rdw_command)
            )
          })
        //local local script
        local_pw.close
      }
    })
    println("Successfully completed!")
  }


}
