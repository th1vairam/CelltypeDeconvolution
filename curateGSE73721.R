#### Load Libraries ####
library(CovariateAnalysis)
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)

library(synapseClient)
library(knitr)
library(githubr)

library(GEOquery)
library(limma)

synapseLogin()

#### Get the latest commit of used files from github ####
thisRepo <- githubr::getRepo(repository = 'th1vairam/CelltypeDeconvolution', ref = 'branch', refName = 'datacuration')
thisFile <- githubr::getPermlink(repository = thisRepo, repositoryPath= 'curateGSE73721.R')

# Get data from synapse
covariates = downloadFile('syn9890517')
fpkm = read.table(synGet('syn9796101')@filePath, header = T, sep = ',', fill = NA) %>%
  dplyr::select(Gene, one_of(covariates$ID))

# Write files to synapse
rownames(fpkm) = toupper(fpkm$Gene)
fpkm$Gene = NULL
lfpkm = log2(fpkm)
lfpkm %>%
  rownameToFirstColumn('hgnc_symbol') %>%
  write.table(file = 'lfpkm.tsv', sep = '\t', quote = F, row.names = F)
expr.obj = File('lfpkm.tsv', parentId = 'syn9796099', name = 'Expression (Log FPKM)')
expr.obj = synStore(expr.obj, 
                    activityName = 'Curate GEO data', 
                    used = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE+73721',
                    executed = thisFile)

fpkm %>%
  rownameToFirstColumn('hgnc_symbol') %>%
  write.table(file = 'fpkm.tsv', sep = '\t', quote = F, row.names = F)
expr.obj = File('fpkm.tsv', parentId = 'syn9796099', name = 'Expression (FPKM)')
expr.obj = synStore(expr.obj, 
                    activityName = 'Curate GEO data', 
                    used = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE+73721',
                    executed = thisFile)