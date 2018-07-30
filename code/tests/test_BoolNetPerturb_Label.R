library(testthat) 
library(BoolNet)

source("BoolNetPerturb_Label.R")
test_results <- test_dir(".", reporter="summary")

data("cellcycle")
int.state <- 162
bin.state <- c(0,1,0,0,0,1,0,1,0,0)
node.names <- cellcycle$genes
names(bin.state) <- node.names
attr <- getAttractors(cellcycle)
label.rules = data.frame(labels = c('G1','S','G2','M'),
           rules = c('Rb & ! (E2F | CycB | CycA)',
                     '(CycD & E2F) &! (Rb | CycB | CycA)',
                     'E2F','(Rb & CycB & CycA) & ! CycD'),
           stringsAsFactors=FALSE   )


test_that("labelState",{
  expect_equal(labelState(int.state, node.names, label.rules), 'G1')
  expect_equal(labelState(bin.state, node.names, label.rules), 'G1')
})

test_that("labelState",{
  expect_equal(labelAttractors(attr, label.rules), 
               c("G1","NoLabel","NoLabel","NoLabel","NoLabel","SG2","SG2","G2" ) )
})
