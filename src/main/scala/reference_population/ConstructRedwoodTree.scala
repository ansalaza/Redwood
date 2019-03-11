package reference_population

import java.io.File

import boopickle.Default._
import utilities.FileHandling._
import utilities.KmerTreeUtils.{Kmers, Leaf, Node, Tree}
import utilities.NewickParser.NwkString
import utilities.ReducedKmerTreeUtils.ReducedKmerTree
import utilities.SketchUtils.{RedwoodSketch, loadRedwoodSketch}

import scala.annotation.tailrec

/**
  * Author: Alex N. Salazar
  * Created on 7-3-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object ConstructRedwoodTree {

  case class Config(
                     tree: File = null,
                     sketchesFile: File = null,
                     prefix: String = null,
                     verbose: Boolean = false,
                     outputDir: File = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("kmer-tree") {
      opt[File]('t', "tree") required() action { (x, c) =>
        c.copy(tree = x)
      } text ("Tree in newick format.")
      opt[String]("prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output file.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      note("\nOPTIONAL\n")
      opt[File]("sketches") action { (x, c) =>
        c.copy(sketchesFile = x)
      } text ("File containing path to all sketches, one per line.")
      opt[Unit]("verbose") action { (x, c) =>
        c.copy(verbose = true)
      } text ("Verbose.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.tree)
      if (config.sketchesFile != null) verifyFile(config.sketchesFile)
      constructRedwoodTree(config)
    }
  }

  def constructRedwoodTree(config: Config): Unit = {
    //construct sketch id -> sketch file
    val (sketches, kmer_length, sketch_size) = {
      if (config.sketchesFile == null) (Map[String, File](), 0, 0)
      else {
        println(timeStamp + "Processing sketches")
        //load paths for sketches
        val tmp = openFileWithIterator(config.sketchesFile).toList.map(new File(_))
        //verify individual sketch file
        tmp.foreach(verifyFile(_))
        //construct sketches map and get kmer length and sketch size
        val (map, klength, slength) = loadSketches(tmp)
        println(timeStamp + "Found " + tmp.size + " sketches with universal kmer-length of " + klength + " and minimum " +
          "sektch size of " + slength)
        (map, klength, slength)
      }
    }

    /**
      * Create a new leaf by attempting load sketch from the constructed sketch map
      * @return Leaf
      */
    def loadLeaf: (String, Double) => Leaf = (id, dist) => {
      //leaf is not in sketch map, return empty leaf
      if(!sketches.contains(id)) new Leaf(id, Set(), 0.0, 0, 0)
      else {
        //load sketch
        val sketch = loadRedwoodSketch(sketches(id))
        //create leaf
        new Leaf(id, sketch.sketch.keySet, dist, sketch.genome_size, sketch.sketch.size)
      }
    }

    //set mutable id variable
    var id = -1

    /**
      * Method to recusrively build tree by recursively parsing newick stree
      * @param nwk Current newick string
      * @param tree Accumulating tree
      * @return Tree
      */
    def buildTree(nwk: NwkString, tree: Tree[Kmers]): Tree[Kmers] = {
      /**
        * Function to update kmers in a given leaf with its companion sub-tree. Given leaf A and companion subtree B,
        * the update kmers for leaf A is: A' = A diff B
        * @return Leaf
        */
      def updateLeaf: (Tree[Kmers], Double, Tree[Kmers]) => Leaf = (_leaf, leaf_dist, subtree) => {
        //load as leaf
        val leaf = _leaf.loadAsLeaf()
        //update leaf's kmers as A' = A - B
        new Leaf("L" + leaf.id, leaf.kmers.diff(subtree.loadKmers()), leaf_dist, leaf.genome_size, leaf.kmers.size)
      }
      //current newick tree is a leaf, return loaded leaf
      if (nwk.isLeaf()) loadLeaf(nwk.str, nwk.distance)
      else {
        //get left and right subtrees
        val (ls, rs) = nwk.getNextSubtrees()
        //build left and right subtrees
        val (lt, rt) = (buildTree(ls, tree), buildTree(rs, tree))
        //update id
        id += 1
        //both subtrees are not leafs, create node as: N = A int B
        if(!ls.isLeaf() && !rs.isLeaf()) new Node(id, lt.loadKmers().intersect(rt.loadKmers()), nwk.distance, lt, rt)
        //left subtree is a leaf, update accordingly
        else if (ls.isLeaf() && !rs.isLeaf()){
          //update left leaf
          val ul = updateLeaf(lt, ls.distance, rt)
          //create node using un-updated leaf kmers but using updated leaf in left subtree
          new Node(id, lt.loadKmers().intersect(rt.loadKmers()), nwk.distance, ul, rt)
        }
          //right subtree is a leaf, update accordingly
        else if(!ls.isLeaf() && rs.isLeaf()){
          //update right leaf
          val ur = updateLeaf(rt, rs.distance, lt)
          //create node using un-updated leaf kmers but using updated leaf in right subtree
          new Node(id, lt.loadKmers().intersect(rt.loadKmers()), nwk.distance, lt, ur)
        }
          //both subtrees are leafs, update accordingly
        else {
          //udpate left and right leafs
          val (ul, ur) = (updateLeaf(lt, ls.distance, rt),updateLeaf(rt, rs.distance, lt))
          //create node using un-updated leaf kmers but using updated leaf in left and right subtrees
          new Node(id, lt.loadKmers().intersect(rt.loadKmers()), nwk.distance, ul, ur)
        }
      }
    }
    println(timeStamp + "Constructing kmer-tree")
    //load newick tree string
    val newick = new NwkString(openFileWithIterator(config.tree).toList.head.drop(1).dropRight(2), 0)
    //construct kmer tree
    val ktree = buildTree(newick, null)
    //construct map of old node ID to new node ID
    val old2new = ktree.inOrderTraversal().zipWithIndex.toMap
    println(timeStamp + "Reducing kmer-tree")
    //create reduced tree along with kmer map
    val (rt, km) = ktree.copy2ReducedKmerTree[Kmers](old2new = old2new)
    println(timeStamp + "Constructed reduced kmer-tree of " + rt.getLeafNamesPostOrder().size + " leafs and " +
      rt.getNodeIDsPostOrder().size + " nodes")
    println(timeStamp + "Writing to disk")
    writeSerialized(Pickle.intoBytes(ReducedKmerTree(rt, km, sketch_size, kmer_length)).array(),
      new File(config.outputDir + "/" + config.prefix + ".rdwt"))
    println("Successfully completed!")
  }

  /**
    * Method to obtain a 3-tuple: (map of sketch id -> sketch file path, kmer length, min sketch size)
    * @param sketches List of sketch files
    * @return 3-tuple as described above
    */
    def loadSketches(sketches: List[File]): (Map[String,File], Int, Int) = {
      @tailrec def _loadSketches(remaining: List[File],
                        smap: Map[String,File],
                        kmer_lengths: Set[Int],
                        sketch_sizes: Set[Int]
                       ): (Map[String, File], Int, Int) = {
        remaining match {
          case Nil => {
            //assert single kmer length
            assert(kmer_lengths.size == 1, "Found different kmer lengths in sketches: " + kmer_lengths.mkString(","))
            if(sketch_sizes.size > 1)
              println(timeStamp + "WARNING: different sketch sizes found: " + sketch_sizes.mkString(","))
            (smap, kmer_lengths.head, sketch_sizes.min)
          }
          case (head::tail) => {
            //get sketch identifier
            val sketch = loadRedwoodSketch(head)
            //assert uniqueness
            assert(!smap.contains(sketch.name), "Dual instance of sketch identifier: " + sketch.name)
            //add to set
            _loadSketches(tail, smap + (sketch.name -> head),
              kmer_lengths + sketch.kmer_length, sketch_sizes + sketch.sketch.size)
          }
        }
      }
      _loadSketches(sketches, Map(), Set(), Set())
    }



}
