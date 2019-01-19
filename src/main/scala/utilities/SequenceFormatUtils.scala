package utilities

import java.io.{File}

import utilities.FileHandling.{getFileExtension, openFileWithIterator, openGzipFileWithIterator}
import utilities.SequenceUtils.{ByteEncoded, Nucleotide, encodeSequence}

import scala.annotation.tailrec
import scala.collection.mutable

/**
  * Author: Alex N. Salazar
  * Created on 18-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object SequenceFormatUtils {

  /**
    * Function to load sequence iterator based on file type (extension)
    *
    * @return SequenceIterator
    */
  def loadSequenceFile: File => SequenceIterator = file => {
    //load sequence iterator based on file extension
    getFileType(file) match {
      //fasta type
      case `fasta_type` => new FastaIterator(file)
      //fastq type
      case `fastq_type` => new FastqIterator(file)
      //fastq gz type
      case `fastq_gz_type` => new FastqIterator(file, true)
      //unrecognized
      case _ => {assert(false, "Unrecognized file type for file " + file.getName); new FastaIterator(file)}
    }
  }

  /**
    * Class for FASTQ.GZ iterator
    *
    * @param file
    * @param gz File is GZ compressed
    */
  class FastqIterator(file: File, gz: Boolean = false) extends SequenceIterator {
    //load iterator (either regular or gz)
    private val iterator = (if (gz) openGzipFileWithIterator(file) else openFileWithIterator(file)).grouped(4)

    def hasNext(): Boolean = iterator.hasNext

    /**
      * Get next FASTQ entry as byte-encoded sequence
      *
      * @return ByteEncoded
      */
    def next(): ByteEncoded = {
      //load next entry
      val entry = iterator.next()
      //sanity echeck
      assert(entry.size == 4, "Unexpected number of lines in fastq entry:\n" + entry.mkString("\n"))
      //byte encode sequence
      encodeSequence(entry(1))
    }
  }


  /**
    * Class for FASTA iterator
    *
    * @param file
    */
  class FastaIterator(file: File) extends SequenceIterator {
    //load iterator
    private val iterator = openFileWithIterator(file)

    /**
      *
      * @return
      */
    def hasNext(): Boolean = iterator.hasNext

    /**
      *
      */
    private var last_entry: String = ""
    private var first: Boolean = true

    /**
      *
      * @return
      */
    def next(): ByteEncoded = {

      /**
        *
        * @param seq
        * @return
        */
      @tailrec def _next(seq: mutable.ArrayBuffer[Nucleotide]): ByteEncoded = {
        //very first line in the file
        if (first) {
          //get first line
          last_entry = iterator.next()
          //update boolean flag
          first = false
          //sanity check
          assert(last_entry.startsWith(">"))
          //go to next iteration
          _next(seq)
        }
        //end of file
        else if (!iterator.hasNext) seq.toArray
        //more lines to process
        else {
          //get current line
          val current = iterator.next()
          //still in the same fasta entry, append to seq
          if (!current.startsWith(">")) _next(seq ++= encodeSequence(current))
          //start of new fasta entry
          else {
            //update last entry header
            last_entry = current
            seq.toArray
          }
        }
      }

      //looad next fasta entry
      val next_seq = _next(mutable.ArrayBuffer.empty[Nucleotide])
      next_seq
    }

  }

  /**
    * Internal strings for classifying fasta or fastq iles
    */
  private val fasta_type = "fasta"
  private val fastq_type = "fastq"
  private val fastq_gz_type = "fastq.gz"

  /**
    * Abstract class for FASTA, FASTQ, FASTQ.GZ class iterators
    */
  abstract class SequenceIterator {
    def hasNext(): Boolean

    def next(): ByteEncoded
  }

  /**
    * Function to get file type for a given file
    *
    * @return String
    */
  private def getFileType: File => String = file => {
    //all fasta extensions
    val fasta_types = List("fa", "fasta", "fna")
    //all fastq extensions
    val fastq_types = List("fastq", "fq")
    //all fastq gz extensions
    val fastq_gz_types = List("fastq.gz", "fq.gz")
    //get file name
    val name = file.getName

    //check for fasta type
    if (fasta_types.exists(x => name.endsWith(x))) fasta_type
    //check for fastq type
    else if (fastq_types.exists(x => name.endsWith(x))) fastq_type
    //check for fastq gz type
    else if (fastq_gz_types.exists(x => name.endsWith(x))) fastq_gz_type
    //unknown type
    else {
      assert(false, "Uncrecognized file type for file " + file.getName)
      ""
    }
  }


}
