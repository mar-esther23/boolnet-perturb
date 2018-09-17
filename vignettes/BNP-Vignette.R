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

## ------------------------------------------------------------------------
mutants <- perturbNetworkFixedNodes(netTh17Treg, label.rules = labelsTh17Treg)

mutants <- mutants[order(-mutants$WT),]  # order by WT in descending order
mutants[mutants==0] <- NA  # replace 0 for NA for plotting
colfunc <- colorRampPalette(c("cyan", "darkblue"))  # color scale
heatmap(data.matrix(mutants),
        Rowv=NA, Colv=NA, 
        col= colfunc(10),scale="none")

## ------------------------------------------------------------------------
data("envTh17Treg")
envTh17Treg

## ------------------------------------------------------------------------
env.attr <- perturbNetworkFixedNodes(netTh17Treg, label.rules = labelsTh17Treg,
                                     genes  = envTh17Treg$nodes, 
                                     values = envTh17Treg$value, 
                                     names  = envTh17Treg$label)

env.attr <- env.attr[order(-env.attr$All),]  # order by WT in descending order
env.attr[env.attr==0] <- NA  # replace 0 for NA for plotting
colfunc <- colorRampPalette(c("cyan", "darkblue"))  # color scale
heatmap(data.matrix(env.attr),
        Rowv=NA, Colv=NA, 
        scale="none", col= colfunc(10))

## ------------------------------------------------------------------------
# Get WT attractors
attr.df <- getAttractors(netTh17Treg)
attr.df <- attractorToDataframe(attr.df)

# Perturb net and get attractors
perturbed.net <- perturbNetwork(netTh17Treg, perturb="functions")
pertubed.attr <- getAttractors(perturbed.net)
pertubed.attr.df <- attractorToDataframe(pertubed.attr)

# Calculate number of attractors not in the original network
diff <- setdiff(attr.df$involvedStates, pertubed.attr.df$involvedStates) 
length(diff)

## ------------------------------------------------------------------------
attr <- getAttractors(netTh17Treg)
verifySyncronousVsAsyncronous(netTh17Treg, attr, label.rules = labelsTh17Treg)

## ------------------------------------------------------------------------
perturbState(netTh17Treg, state = 0, 
             genes=c("IL21e","TGFBe"), 
             result = "trajectory")

## ------------------------------------------------------------------------
cellfate <- cellFateMap(netTh17Treg, label.rules = labelsTh17Treg)
head(cellfate)

## ------------------------------------------------------------------------
res.der <- derrida(netTh17Treg)
plot(x=names(res.der), y=res.der, 
     ylim = c(0,length(res.der)), 
     xlab = "h_t", ylab = "h_t+1")
abline(0,1)

