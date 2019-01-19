name := "Redwood"

version := "0.1"

scalaVersion := "2.12.8"

libraryDependencies ++= Seq(
  "org.scalatest" %% "scalatest" % "3.0.4" % "test",
  "com.github.scopt" % "scopt_2.12" % "3.7.0",
  "io.suzaku" %% "boopickle" % "1.3.0",
  "com.github.alexandrnikitin" %% "bloom-filter" % "latest.release"
)