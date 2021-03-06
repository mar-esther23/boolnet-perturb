% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BNP_Helper.R
\name{isGeneInput}
\alias{isGeneInput}
\title{Determine if a node is an input of the network.}
\usage{
isGeneInput(gene, net)
}
\arguments{
\item{gene}{node to asses}

\item{net}{BoolNet network}
}
\value{
Boolean is the node an input in the network.
}
\description{
Determine if a node is an input of the network.
}
\examples{
isGeneInput("CycD", cellcycle)

}
\keyword{internal}
