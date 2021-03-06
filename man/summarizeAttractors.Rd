% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BNP_Summarize.R
\name{summarizeAttractors}
\alias{summarizeAttractors}
\title{Summarize attractors of a dataframe ignoring certain nodes.}
\usage{
summarizeAttractors(attr, genes = NULL, net = NULL, returnBasin = F)
}
\arguments{
\item{attr}{BoolNet attractor object}

\item{genes}{list of gene names to ignore, by default it detects inputs. To perturb multiple nodes at the same time use a vector inside the list.}

\item{net}{net to detect inputs from}

\item{returnBasin}{default FALSE, return a column with the combined basin size}
}
\value{
dataframe, each column corresponds to the number of attractor, state, or node. If Boolean=FALSE return dataframe, each column corresponds to a property of the attractor.
}
\description{
Summarizes the attractors as a dataframe, ignoring nodes that are not relevant for the analysis (input nodes by default). Returns a dataframe in binary format. It can return the combined basin size of the attractors. The function does not summarize cycles to avoid mistakes.
}
\examples{
data(netTh17Treg)
attr <- getAttractors(netTh17Treg)
summarizeAttractors(attr, net=netTh17Treg, returnBasin=T)
#    attractor state IL2 RORGT STAT3 FOXP3 TGFB basinSize
# 1          1     1   0     0     0     0    0         7
# 2          2     1   1     0     0     0    0        21
# 3          3     1   0     0     1     0    0        60
# 4          4     1   0     0     0     0    1        27
# 5          7     1   1     0     0     1    1        41
# 6         14     1   0     1     1     0    1        62
# 7         20     1   0     0     0     0    0        16
# 8         20     2   1     0     1     0    0        NA
# 9         21     1   0     1     0     0    1         6
# 10        21     2   1     0     1     0    1        NA
# 11        22     1   0     1     0     0    1        16
# 12        22     2   1     0     1     0    1        NA

summarizeAttractors(attr, c("IL2","STAT3","IL2e","IL21e"))

}
