library(tidyverse)
library(msigdbr)
library(GSVA) 
library(GSEABase)
library(Seurat)

genesets <- msigdbr(species = 'Homo sapiens',category = 'H')
genesets <- subset(genesets,select = c('gs_name','gene_symbol')) %>% as.data.frame()

genesets <- split(genesets$gene_symbol, genesets$gs_name)

library('dior')
rna = dior::read_h5(file='./raw/CRC-P19-T.h5', target.object = 'seurat')
rna = NormalizeData(rna)
expr <- as.matrix(rna@assays$RNA@data)
gsva_mat <- gsva(expr = expr,
               gset.idx.list = genesets,
               kcdf = "Gaussian",#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
               verbose = T,
               parallel.sz = 10)

gsva.df <- data.frame(Genesets=rownames(gsva_mat),gsva_mat,check.rows = F)
write.csv(gsva_mat, file = "CRC-P19-T.gsva.csv", row.names = TRUE)
