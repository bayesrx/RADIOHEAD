## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  # install the package (devtools package needed)
#  if(!require(devtools)) install.packages("devtools")
#  devtools::install_github('shariq-mohammed/RADIOHEAD')

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

## ----message=FALSE, warning=FALSE---------------------------------------------
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!requireNamespace("GSVA", quietly = TRUE)) BiocManager::install("GSVA")
library(GSEABase)
library(GSVA)

## ----message=FALSE, warning=FALSE---------------------------------------------
gene.set = getGmt("https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=GO_EXOCYTOSIS&fileType=gmt")

## -----------------------------------------------------------------------------
gene.expr = matrix(rnorm(5*1000), ncol=5)
colnames(gene.expr) = paste0('Subject', 1:5)
rownames(gene.expr) = paste0('Gene', 1:1000)

# Set row names as the gene names from the pathway
gene.ids = gene.set@.Data[[1]]@geneIds
rownames(gene.expr)[1:length(gene.ids)] = gene.ids

## -----------------------------------------------------------------------------
pathwayScores = gsva(gene.expr, gene.set, method='gsva', verbose=FALSE)
pathwayScores

## -----------------------------------------------------------------------------
n = 5 # number of subjects
for(i in 1:5){
  assign(paste0('mri',i), array(rnorm(10*10*10), dim = c(10,10,10)))
  assign(paste0('tumor',i), array(rep(sample(0:3, 4),c(125,125,125,1000-(3*125))), dim = c(10,10,10)))
}

## -----------------------------------------------------------------------------
# sub-region 1
tum.1 = lapply(1:n,
               function(i){
                 mri = eval(parse(text=paste0('mri',i)))
                 tumor = eval(parse(text=paste0('tumor',i)))
                 x = c(mri[tumor==1])
               })

# sub-region 2
tum.2 = lapply(1:n,
               function(i){
                 mri = eval(parse(text=paste0('mri',i)))
                 tumor = eval(parse(text=paste0('tumor',i)))
                 x = c(mri[tumor==2])
               })

# sub-region 3
tum.3 = lapply(1:n,
               function(i){
                 mri = eval(parse(text=paste0('mri',i)))
                 tumor = eval(parse(text=paste0('tumor',i)))
                 x = c(mri[tumor==3])
               })

## -----------------------------------------------------------------------------
pc.scores.1 = geomPCA(tum.1)
pc.scores.2 = geomPCA(tum.2)
pc.scores.3 = geomPCA(tum.3)

## -----------------------------------------------------------------------------
pc.scores = cbind(pc.scores.1, pc.scores.2, pc.scores.3)
pc.groups = rep(c(1,2,3), c(length(pc.scores.1),
                            length(pc.scores.2), length(pc.scores.3)))

