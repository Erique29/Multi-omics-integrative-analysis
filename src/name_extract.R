setwd('./MetaboAnalyst/')
library(readxl)
library(dplyr)
allname<-read_excel('./MetaboAnalyst_result_pathway.xlsx',sheet = 'all_name')
genes<-read_excel('./MetaboAnalyst_result_pathway.xlsx',sheet = 'result_genes')
cmpds<-read_excel('./MetaboAnalyst_result_pathway.xlsx',sheet = 'result_cmpds')
gene_name_map <-read_excel('./MetaboAnalyst_result_pathway.xlsx',sheet = 'gene_name_map')
cmpds_name_map <-read_excel('./MetaboAnalyst_result_pathway.xlsx',sheet = 'cmpds_name_map')
Entriz_symbol <-read_excel('./MetaboAnalyst_result_pathway.xlsx',sheet = 'Entriz_symbol')

merge<-merge(allname, genes, by.x = 'Name',by.y= 'ID',all.x = T, sort = F) %>% 
  merge(cmpds, by.x = 'Name',by.y= 'ID',all.x = T, sort = F) %>% 
  merge(gene_name_map, by.x = 'Name',by.y= 'Symbol',all.x = T, sort = F) %>% 
  merge(cmpds_name_map, by.x = 'Name',by.y= 'KEGG',all.x = T, sort = F) %>% 
  merge(Entriz_symbol, by.x = 'Entrez',by.y= 'GeneID',all.x = T, sort = F)
library(rio)
export(merge,"clipboard")
