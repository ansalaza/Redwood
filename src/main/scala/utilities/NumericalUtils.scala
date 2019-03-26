package utilities

import scala.math.BigDecimal

/**
  * Author: Alex N. Salazar
  * Created on 16-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object NumericalUtils {

  import Numeric.Implicits._

  /**
    * Function to round double to 2-decimal points
    * @return Double
    */
  def round: Double => Double = value => BigDecimal(value).setScale(2, BigDecimal.RoundingMode.HALF_UP).toDouble

  def multiplyList: List[Int] => Int = values => values.foldLeft(1)((acc,x) => acc * x)

  /**
    * Recursive method to compute binomial coefficient
    *
    * @param n
    * @param k
    * @return Int
    */
  def choose(n: Int, k: Int): (List[Int], List[Int]) = {

    /**
      * Functioni to get a list of number representing the factorial expression
      * @param f Number
      * @return List[Int]
      */
    def factorialExpression(f: Int): List[Int] = (1 to f).toList

    //compute n choose k
    if (k == 0 || k == n) (List(1), List(1))
    else {
      //construct map integer -> integer count representing the entire denominator factorial expression
      val denom_map = (factorialExpression(n - k) ::: factorialExpression(k)).groupBy(identity).mapValues(_.size)
      //cancel out numerator and get final count of denominator
      val (remaining_denom, remaining_num) = {
        //iterate through numerator factorial expression
        factorialExpression(n).foldLeft((denom_map, List[Int]())){case ((denom_remaining, num_remaining), number) => {
          //check count in the denominator
          val count = denom_map.get(number)
          //does not appear in the denominator, stays in numerator
          if(count.isEmpty) (denom_remaining, number :: num_remaining)
          //appears in the denominator, subtract count by 1
          else (denom_remaining + (number -> (count.get - 1)), num_remaining)
        }}
      }
      //compute product of remaining numerator over remaining denominator
      (remaining_num, remaining_denom.filter(_._2 > 0).keys.toList)
    }
  }


  /**
    * Function to compute minimum of two numbers
    *
    * @param num Implicit numeric class type
    * @tparam T Numeric type
    * @return Minimum of two given numbers
    */
  def min[T](a: T, b: T)(implicit num: Numeric[T]): T = {
    import num._
    if (a < b) a else b
  }

  /**
    * Method to compute the median of a given collection of numbers
    *
    * @param xs Collection of numbers
    * @tparam T Numeric type
    * @return Double
    */
  def median[T: Numeric](xs: Iterable[T]): Double = {
    //sort values
    val sorted = xs.toList.sorted
    //size
    val size = sorted.size
    //odd number of values, get middle
    if (sorted.size % 2 == 1) sorted(size / 2).toDouble
    else {
      //get left and right values
      val (left, right) = (sorted((size / 2) - 1).toDouble, sorted(size / 2).toDouble)
      //get middle
      (left + right) / 2
    }
  }

  /**
    * Method to compute mean of a collection of numbers
    *
    * @param xs Collection of numbers
    * @tparam T Numeric type
    * @return Mean value of collection
    */
  def mean[T: Numeric](xs: Iterable[T]): Double = xs.sum.toDouble / xs.size

  /**
    * Method to compute variance in a collection of numbers (note: converts collection to parallel collection)
    *
    * @param xs Collection of numbers
    * @tparam T Numeric type
    * @return Variance of collection
    */
  def variance[T: Numeric](xs: Iterable[T]): Double = {
    val mu = mean(xs)
    xs.map(a => {
      math.pow(a.toDouble - mu, 2)
    }).sum / xs.size
  }

  /**
    * Method to compute standard deviation of a colletion of numbers
    *
    * @param xs Collection of numbers
    * @tparam T Numeric type
    * @return Standard deviation of collection
    */
  def stdDev[T: Numeric](xs: Iterable[T]): Double = math.sqrt(variance(xs))


}
