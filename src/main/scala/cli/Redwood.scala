package cli

/**
  * Author: Alex N. Salazar
  * Created on 3-1-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object Redwood {

  def main(args: Array[String]): Unit = {

    val help = (
      "Usage: java -jar ptolemy2.jar [tool]\n\n" +
        "REDWOOD TOOLS\n" +
        "sketch            Create sketch of a given FASTQ file(s).\n" +
        "sketch-cluster    Creates automated scripts for 'sketch' command for a cluster.\n" +
        "distance          Compute the jaccard distance matrix of two or more sketches.\n"+
        "tree              Build kmer-based population tree from a given distance matrix.\n" +
        "tree-metrics      Obtain summary metrics for given kmer tree.\n" +
        "classify          Query kmers from a sketch file againts a given kmer-tree.\n\n"

      )

    if (args.length == 0) {
      println(help)
    } else {
      args(0) match {
        case "sketch" => sketch.BuildSketch.main(args.drop(1))
        case "sketch-cluster" => sketch.BuildSketchUsingCluster.main(args.drop(1))
        case "distance" => compare.SketchDistance.main(args.drop(1))
        case "tree" => reference_population.ConstructPopulationTree.main(args.drop(1))
        case "tree-metrics" => reference_population.TreeMetrics.main(args.drop(1))
        case "classify" => compare.SketchClassifier.main(args.drop(1))
        case _ => println(help)
      }
    }
  }
}
