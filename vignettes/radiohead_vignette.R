## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  # Change the directory to where the package is saved
#  setwd('C:/Users/xyz/R Package')
#  
#  # install the package (devtools package needed)
#  devtools::install('./RADIOHEAD')

## -----------------------------------------------------------------------------
library(RADIOHEAD)

## -----------------------------------------------------------------------------
y = scale(c_pathway_scores[,"EXOCYTOSIS"], center = T, scale = F)

## -----------------------------------------------------------------------------
X = scale(pc_scores, center = T, scale = T)

## -----------------------------------------------------------------------------
groups = pc_groups

## -----------------------------------------------------------------------------
res = groupSS(y, X, groups, Nmcmc=1000)

## -----------------------------------------------------------------------------
# 'associations' outputs the effect size of the significant associations.
associations = fdr_var_selection(res$b, res$x_cnames)
associations

