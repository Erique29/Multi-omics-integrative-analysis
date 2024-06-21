term_list <- list()
term_list=strsplit(drawing$geneID,split = "/")
names(term_list)=drawing$KEGGID
expr = all_compare_2023_7_11
rownames(expr) = expr$gene_id
export<-list()
export <-
  lapply(term_list, function(x) {
    as.data.frame(expr[x,])
  })
library(rio)
export(export,xlsx_file <- tempfile(fileext = ".xlsx"))
