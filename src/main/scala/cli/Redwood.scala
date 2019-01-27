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
        "SKETCH\n" +
        "sketch            Create a sketch for given sequence file.\n" +
        "sketch-cluster    Create automated scripts for 'sketch' command for a cluster.\n" +
        "sketch-serially   Create sketch for multiple libraries serially.\n\n" +
        "DISTANCE\n" +
        "distance          Compute the mash distance matrix of two or more sketches.\n"+
        "update            Update existing mash distance matrix with new sketches.\n\n" +
        "KMER TREE\n" +
        "tree              Build kmer-based population tree from a given distance matrix.\n" +
        "tree-drawer       Visualize a given tree.\n"+
        "tree-metrics      Obtain summary metrics for given kmer tree.\n" +
        "query             Query kmers from a sketch against a given kmer-tree.\n" +
        "query-serially    Query multiple sketches serially.\n\n"
      )

    if (args.length == 0) {
      println(help)
    } else {
      args(0) match {
        case "sketch" => sketch.BuildSketch.main(args.drop(1))
        case "sketch-cluster" => sketch.BuildSketchUsingCluster.main(args.drop(1))
        case "sketch-serially" => sketch.BuildSketchSerially.main(args.drop(1))
        case "distance" => compare.SketchDistance.main(args.drop(1))
        case "update" => compare.UpdateMatrix.main(args.drop(1))
        case "tree" => reference_population.ConstructPopulationTree.main(args.drop(1))
        case "tree-drawer" => vizualization.TreeTracing.main(args.drop(1))
        case "tree-metrics" => reference_population.TreeMetrics.main(args.drop(1))
        case "query" => compare.SketchQuery.main(args.drop(1))
        case "query-serially" => compare.SketchQuerySerially.main(args.drop(1))
        case "prep-data" => sketch.PrepareData.main(args.drop(1))
        case _ => println(help)
      }
    }
  }
}
