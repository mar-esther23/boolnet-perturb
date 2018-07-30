library(testthat) 
library(BoolNet)

source("BoolNetPerturb_Dataframe.R")
test_results <- test_dir(".", reporter="summary")

data("cellcycle")
int.state <- 162
bin.state <- c(0,1,0,0,0,1,0,1,0,0)
names(bin.state) <- cellcycle$genes
attr <- getAttractors(cellcycle)

test_that("bin2intState",{
  expect_equal( attractor2dataframe(attr),
                data.frame(involvedStates = c('162','25/785/849/449/389/141/157'),
                           basinSize = c(512,512), stringsAsFactors=FALSE   ) )
})

test_that("bin2intState",{
  expect_equal( attractor2dataframeBin(attr),
                data.frame(attractor=c(1,2,2,2,2,2,2,2), state=c(1,1,2,3,4,5,6,7),
                           CycD=c(0,1,1,1,1,1,1,1), Rb=c(1,0,0,0,0,0,0,0),
                           E2F=c(0,0,0,0,0,1,1,1), CycE=c(0,1,0,0,0,0,1,1),
                           CycA=c(0,1,1,1,0,0,0,1), p27=c(1,0,0,0,0,0,0,0),
                           Cdc20=c(0,0,0,1,1,0,0,0), Cdh1=c(1,0,0,0,1,1,1,1),
                           UbcH10=c(0,0,1,1,1,1,0,0), CycB=c(0,0,1,1,0,0,0,0), stringsAsFactors=FALSE   ))
})
