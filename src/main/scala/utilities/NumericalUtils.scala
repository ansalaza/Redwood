package utilities

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
    * Recursive method to compute binomial coefficient
    *
    * @param n
    * @param k
    * @return Int
    */
  def choose(n: Int, k: Int): Int = {
    /**
      * Method to compute factorial of an INT
      *
      * @param f Input number
      * @return Int
      */
    def fact(f: Int): Int = (1 to f).foldLeft(1)(_ * _)
    //compute n choose k
    if (k == 0 || k == n) 1 else fact(n) / ((fact(k) * fact(n - k)))
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
