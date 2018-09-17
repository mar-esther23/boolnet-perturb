## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE----------------------------------------------------------
#  library(devtools)
#  install_github("mar-esther23/boolnet-perturb", force = TRUE)

## ------------------------------------------------------------------------
library("BoolNetPerturb")

data(netTh17Treg)
netTh17Treg

## ------------------------------------------------------------------------
attr <- getAttractors(netTh17Treg)
attr.df <- attractorToDataframe(attr, Boolean = TRUE)
head(attr.df)

## ------------------------------------------------------------------------
data("labelsTh17Treg")
labelsTh17Treg

## ------------------------------------------------------------------------
labels <- labelAttractors(attr, labelsTh17Treg)
table(labels)

