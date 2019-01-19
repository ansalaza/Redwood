package utilities


/**
  * Author: Alex N. Salazar
  * Created on 3-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object SequenceUtils {

  /**
    * Byte encodings for individual nucleotides
    */
  private val A_byte = 0.toByte
  private val C_byte = 1.toByte
  private val G_byte = 2.toByte
  private val T_byte = 3.toByte
  private val N_byte = 4.toByte

  /**
    * Type alias for encoding a DNA sequence into bytes
    */
  type ByteEncoded = Array[Nucleotide]
  type Nucleotide = Byte

  /**
    * Method to determine whether a given byte-encoded sequence only contains valid nucleotides (e.g. A,C,G,T)
    *
    * @return Boolean
    */
  def isValid: ByteEncoded => Boolean = seq => seq.forall(_ != N_byte)


  /**
    * Method to compute the reverse-complement of a given byte-encoded sequence
    *
    * @return Array[Byte]
    */
  def reverseComplement: ByteEncoded => ByteEncoded = seq => {

    /**
      * Function to get complemental byte of a given byte
      *
      * @return Byte
      */
    def complement: Byte => Byte = nt => nt match {
      case A_byte => T_byte
      case T_byte => A_byte
      case C_byte => G_byte
      case G_byte => C_byte
      case _ => N_byte
    }

    seq.reverse.map(complement(_))
  }

  /**
    * Encode a DNA sequence into bytes
    *
    * @return
    */
  def encodeSequence: String => ByteEncoded = seq => {
    seq.toCharArray.map(nt => {
      nt.toUpper match {
        case 'A' => A_byte
        case 'T' => T_byte
        case 'C' => C_byte
        case 'G' => G_byte
        case _ => N_byte
      }
    })
  }

}
