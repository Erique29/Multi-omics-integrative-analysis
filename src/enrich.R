library(readxl)

load('enrich.RData')

deg <- read_excel("deg.xlsx")

View(deg)
ce_lm_rbh <- readRDS("D:/diske/lm_ph_fl/04deseq/enrichment/rbh_test/ce_lm_rbh.rds")
View(ce_lm_rbh)
ce_lm_rbh=ce_lm_rbh[,-2]
ce_lm_rbh=ce_lm_rbh[,-1]
ce_lm_rbh=unique(ce_lm_rbh)
deg=merge(deg,ce_lm_rbh,by.x = "Geneid",by.y = "lm_gene",all.x = T)
View(deg)
write.csv(deg, "deg_anno.csv")
write.csv(deg,"deg_anno.csv")
library(readxl)
deg_up <- read_excel("D:/diske/xxx/rna_seq/deg.xlsx",
                     sheet = "up")
View(deg_up)
deg_down <- read_excel("D:/diske/xxx/rna_seq/deg.xlsx", sheet = "down")
load("D:/diske/lm_ph_fl/04deseq/enrichment/enrich_dat.RData")
rm(ph46)
library(dplyr)
library(clusterProfiler)
library(stringr)
up_kegg = enricher(deg_up$Geneid, TERM2GENE=keggid2gene, TERM2NAME=pathway2name,pvalueCutoff = 1,pAdjustMethod = "BH",qvalueCutoff = 1)
down_go = enricher(deg_up$Geneid, TERM2GENE=goid2gene, TERM2NAME=term2name,pvalueCutoff = 1,pAdjustMethod = "BH",qvalueCutoff = 1)
up_kegg = enricher(deg_up$Geneid, TERM2GENE=keggid2gene, TERM2NAME=pathway2name,pvalueCutoff = 1,pAdjustMethod = "BH",qvalueCutoff = 1)
up_go = enricher(deg_up$Geneid, TERM2GENE=goid2gene, TERM2NAME=term2name,pvalueCutoff = 1,pAdjustMethod = "BH",qvalueCutoff = 1)
down_kegg = enricher(deg_down$Geneid, TERM2GENE=keggid2gene, TERM2NAME=pathway2name,pvalueCutoff = 1,pAdjustMethod = "BH",qvalueCutoff = 1)
down_go = enricher(deg_down$Geneid, TERM2GENE=goid2gene, TERM2NAME=term2name,pvalueCutoff = 1,pAdjustMethod = "BH",qvalueCutoff = 1)
write.csv(down_kegg@result,"down_kegg.csv")
write.csv(down_go@result,"down_go.csv")
write.csv(up_kegg@result,"up_kegg.csv")
write.csv(up_go@result,"up_go.csv")
