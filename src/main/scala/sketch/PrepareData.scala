package sketch

import java.io.File

import utilities.FileHandling.{verifyDirectory, verifyFile, openFileWithIterator}

/**
  * Author: Alex N. Salazar
  * Created on 23-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object PrepareData {

  case class Config(
                     inputFile: File = null,
                     outputDir: File = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("build-sketch") {
      opt[File]('i', "input-file") required() action { (x, c) =>
        c.copy(inputFile = x)
      } text ("List of all sequence files, one per line.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.inputFile)
      prepareData(config)
    }
  }

  def prepareData(config: Config): Unit = {
    openFileWithIterator(config.inputFile).foldLeft((Set[String](), 0))((set, line) => {
      val file = new File(line)
      verifyFile(file)
      val id = {
        val tmp = openFileWithIterator(file).foldLeft(List[String]())((acc, fasta_line) => {
          if(!fasta_line.startsWith(">")) acc else fasta_line.substring(1).split("\\s+").head :: acc
        })
        val _id = tmp.last + "+" + set._2
        assert(tmp.nonEmpty, "Empty FASTA headers: " + line)
        assert(!set._1(_id), "Non-unique ID: " + _id + " " + line)
        _id
      }
      println(id + "\t" + line)
      (set._1 + (id), set._2 + 1)
    })
  }


}
