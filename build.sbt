name := "Redwood"

version := "0.1"

scalaVersion := "2.12.8"

resolvers += Resolver.bintrayRepo("underscoreio", "training")

libraryDependencies ++= Seq(
  "org.scalatest" %% "scalatest" % "3.0.4" % "test",
  "com.github.scopt" % "scopt_2.12" % "3.7.0",
  "com.github.alexandrnikitin" %% "bloom-filter" % "latest.release",
  "underscoreio" %% "doodle" % "0.8.2"
)
