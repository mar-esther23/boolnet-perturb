library(testthat) 
library(BoolNet)

source("BoolNetPerturb_Helper.R")
test_results <- test_dir(".", reporter="summary")

data("cellcycle")
int.state <- 162
bin.state <- c(0,1,0,0,0,1,0,1,0,0)
names(bin.state) <- cellcycle$genes

test_that("int2binState",{
  expect_equal( int2binState(int.state,cellcycle$genes),bin.state )
})

test_that("bin2intState",{
  expect_equal( bin2intState(bin.state),int.state )
})

test_that("validateState",{
  expect_equal(validateState(int.state, cellcycle$genes),bin.state)
  expect_equal(validateState(as.character(int.state), cellcycle$genes),bin.state)
  expect_equal(validateState(bin.state, cellcycle$genes),bin.state)
  expect_that(validateState("162/162", cellcycle$genes),throws_error() )
})

test_that("isGeneInput",{
  expect_equal(isGeneInput("CycD",cellcycle),T)
  expect_equal(isGeneInput("CycA",cellcycle),F)
  expect_equal(isGeneInput("Cdc20",cellcycle),F)
})

