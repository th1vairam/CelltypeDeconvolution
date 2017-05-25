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
library(xlsx)

synapseLogin()

#### Get the latest commit of used files from github ####
thisRepo <- githubr::getRepo(repository = 'th1vairam/CelltypeDeconvolution', ref = 'branch', refName = 'datacuration')
thisFile <- githubr::getPermlink(repository = thisRepo, repositoryPath= 'curateGSE52564.R')

#### Get data from synapse ####
geo.gse = getGEO('GSE52564')
pheno.data = pData(geo.gse[[1]]) %>%
  dplyr::select(characteristics_ch1.2, geo_accession, source_name_ch1) %>%
  dplyr::rename(CellType = characteristics_ch1.2, GSE = geo_accession, Tissue = source_name_ch1) %>%
  dplyr::mutate(CellType = gsub('cell type: ', '', CellType)) %>%
  dplyr::filter(CellType != '') %>%
  dplyr::mutate(CellType = factor(CellType, 
                                  levels = c("Astrocyte", "neuron", "oligodendrocyte precursor cells",
                                             "newly formed oligodendrocytes", "myelinating oligodendrocytes",
                                             "microglia", "endothelial cells"),
                                  labels = c("Astrocyte" = "Astrocyte", "endothelial cells" = "endothelial", 
                                             "microglia" = "microglia", "myelinating oligodendrocytes" = "oligodendrocytes",
                                             "neuron" = "neurons", "newly formed oligodendrocytes" = "oligodendrocytes", 
                                             "oligodendrocyte precursor cells" = "OPC"))) %>%
  dplyr::filter(CellType %in% c('oligodendrocytes', 'astrocytes', 'OPC', 'microglia', 'neurons', 'endothelial'))

# Get supplementary files from GEO
tmp = synQuery('select * from file where parentId == "syn9796071"')
ind = grep('.xls', tmp$file.name)
fpkm = lapply(tmp$file.id[ind], function(id){
  obj = synGet(id)
  dat = read.xlsx(obj@filePath, 1)
  colnames(dat) = c('Gene.Names',paste(str_split(obj$properties$name,'')[[1]][1:10], collapse = ''))
  return(dat)
}) %>%
  plyr::join_all(type = 'full') %>%
  dplyr::select(Gene.Names, one_of(levels(droplevels(pheno.data$GSE))))

# Modify covariates
covariates = pheno.data[,c('Tissue', 'CellType')]
covariates = as.data.frame(lapply(covariates, factor))
rownames(covariates) = pheno.data$GSE

# Write files to synapse
rownames(fpkm) = toupper(fpkm$Gene.Names)
fpkm$Gene.Names = NULL
lfpkm = log2(fpkm)
lfpkm %>%
  rownameToFirstColumn('hgnc_symbol') %>%
  write.table(file = 'lfpkm.tsv', sep = '\t', quote = F, row.names = F)
expr.obj = File('lfpkm.tsv', parentId = 'syn9796071', name = 'Expression (Log FPKM)')
expr.obj = synStore(expr.obj, 
                    activityName = 'Curate GEO data', 
                    used = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE+73721',
                    executed = thisFile)

fpkm %>%
  rownameToFirstColumn('hgnc_symbol') %>%
  write.table(file = 'fpkm.tsv', sep = '\t', quote = F, row.names = F)
expr.obj = File('fpkm.tsv', parentId = 'syn9796071', name = 'Expression (FPKM)')
expr.obj = synStore(expr.obj, 
                    activityName = 'Curate GEO data', 
                    used = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE+73721',
                    executed = thisFile)

covariates %>%
  rownameToFirstColumn('SampleID') %>%
  write.table(file = 'covariates.tsv', sep = '\t', quote = F, row.names = F)
covariates.obj = File('covariates.tsv', parentId = 'syn9796071', name = 'Covariates')
covariates.obj = synStore(covariates.obj, 
                          activityName = 'Curate GEO data', 
                          used = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE+52564',
                          executed = thisFile)