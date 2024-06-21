BiocManager::install("enrichplot")
#Pcoa
library(ggalt)
library(MicrobiotaProcess)
library(patchwork)
library(vegan)
library(ggplot2)
library(ggforce)

##=============================


#==================================================
#DEGs enrichment

library(DESeq2)
library(readxl)
#setwd('d:/Transcriptome/')

readcount<-read_excel('counts.xlsx') 
readcount1<-readcount[,-2]

readcount1<-readcount1[,-1]
rownames(readcount1)<-readcount$Geneid

metaData <- read_excel('counts.xlsx', sheet = 'meta') 

dds <- DESeqDataSetFromMatrix(countData=readcount1, 
                              colData=metaData, 
                              design=~bacteria)
dds <- DESeq(dds#,test = 'Wald', fitType = 'mean'
             )

res <- results(dds)
#res2<-subset(res, padj<.05 & abs(log2FoldChange)>1)
#res3 = cbind(as(res2, "data.frame"))

res2<-subset(res, padj<.05 & abs(log2FoldChange)> 0.3785)
res3 = cbind(as(res2, "data.frame"))


#write.csv(res2, file = "deseq2_results_0.3785.csv")


#vol
#reset par
# par(mfrow=c(1,1))
# # Make a basic volcano plot
# with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6,6)))
# 
# # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
# with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# with(subset(res, padj<.05 & abs(log2FoldChange)>0.585), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#PCA
vsdata <- vst(dds, blind=FALSE)
p.pca<-plotPCA(vsdata, intgroup="dex")+
  theme_bw() + 
  scale_color_gradientn(colours = rainbow(5)) + 
  #xlab(paste0('PC1 (',round(xVar,2),' %)')) + 
  #ylab(paste0('PC2 (',round(yVar,2),' %)')) + 
  theme(axis.text = element_text(color = 'black'),
        legend.position="top",
        legend.title=element_blank())+
  scale_colour_manual(name = "Group", values = c('#1BB6AFFF', '#EE6100FF'),
                      #labels=c("Unconditioned\n       media", "   E. coli-\nconditioned\n    media", "   S. algae-\nconditioned\n    media")
  )+
  geom_mark_ellipse(#data = pCoaVecs2, #aes(x = PCo1, y = PCo2),
    expand = 0.02, radius = 0.01,tol = 0.01,linetype=2,
    linewidth=0.2)
p.pca


#PCOA
library(ggalt)
library(MicrobiotaProcess)
library(patchwork)
library(vegan)
library(ggplot2)
library(ggforce)

# deg <- read_excel("deseq2_results_0.585.xlsx")
# deg2<- merge(deg,readcount,by='Geneid',sort=F)
# pcoa_deg<-deg2[,-c(1:8)] %>% t()
# 
# 
# 
# distMatrix <- vegdist(pcoa_deg,method = "bray") 
# pCoa <- cmdscale(distMatrix, eig = T,k = 2 )
# varExp <- (eigenvals(pCoa)/sum(eigenvals(pCoa)))[1:2]
# xVar <- as.numeric(varExp[1]*100)
# yVar <- as.numeric(varExp[2]*100)
# 
# pCoaVecs <- as.data.frame(pCoa$points)
# colnames(pCoaVecs) <- paste0("PCo",c(1:2))
# pCoaVecs$ID <- row.names(pCoaVecs)
# 
# meta<-read_excel('trans_pcoa.xlsx',sheet = 'meta')
# pCoaVecs2<-merge(pCoaVecs,meta,by='ID')
# 
# #PERMANOVA adonis2
# #adonis_result <- adonis2(pcoa_data1 ~ group, data=meta, permutations=999)
# 
# library(paletteer)
# #pCoaVecs2$group <- factor(pCoaVecs2$group, levels = c('Unconditioned media','S.algae-conditioned media','E.coli-conditioned media'))
# 
# p_pcoa<-ggplot(pCoaVecs2,aes(x=PCo1,y=PCo2,color=group)) + 
#   geom_point(size=4) + 
#   theme_bw() + 
#   scale_color_gradientn(colours = rainbow(5)) + 
#   xlab(paste0('PC1 (',round(xVar,2),' %)')) + 
#   ylab(paste0('PC2 (',round(yVar,2),' %)')) + 
#   theme(axis.text = element_text(color = 'black'),
#         legend.position="top",
#         legend.title=element_blank())+
#   scale_colour_manual(name = "Group", values = c('#1BB6AFFF', '#EE6100FF'),
#                       #labels=c("Unconditioned\n       media", "   E. coli-\nconditioned\n    media", "   S. algae-\nconditioned\n    media")
#   )+
#   geom_mark_ellipse(#data = pCoaVecs2, #aes(x = PCo1, y = PCo2),
#     expand = 0.02, radius = 0.01,tol = 0.01,linetype=2,
#     linewidth=0.2)#+
# #annotate("text",x=-0.25,y=0.2,label=paste("p= ", adonis_result$`Pr(>F)`[1]),size=4)
# 
# 
# p_pcoa
ggsave('DEG_pca.pdf',p.pca,height=4.2,width=4)

#vol
library(EnhancedVolcano)
keyvals <- ifelse(
  res$log2FoldChange < -0.3785, '#1BB6AFFF',
  ifelse(res$log2FoldChange > 0.3785, '#EE6100FF',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#EE6100FF'] <- 'Up'
names(keyvals)[keyvals == 'black'] <- 'NS'
names(keyvals)[keyvals == '#1BB6AFFF'] <- 'Down'
p.vol <- EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                colCustom = keyvals,
                selectLab = rownames(res)[which(names(keyvals) %in% c('Up', 'Down'))],
                title = "",
                subtitle = '',
                caption ='',
                pointSize = 2,
                labSize = 2,
                pCutoff = 0.05,
                FCcutoff = 0.3785,
                widthConnectors = 0.5,
                gridlines.major = F,
                gridlines.minor = FALSE,
                border='full',
                xlim = c(-4,4),
                xlab = bquote(~Log[2] ~ "fold change" ~ italic('S. algae') ~ '/' ~ italic('E. coli')),
                ylab = bquote(~-Log[10] ~ italic(P) ~ Value),
                labCol = "black"
                )
p.vol

ggsave('DEG_vol.pdf',p.vol,height=7.2,width=7)



##enrichment

library(readxl)

load('enrich.RData')

deg <- read_excel("deseq2_results_0.3785.xlsx")

#View(deg)
#ce_lm_rbh <- readRDS("D:/diske/lm_ph_fl/04deseq/enrichment/rbh_test/ce_lm_rbh.rds")
#View(ce_lm_rbh)
#ce_lm_rbh=ce_lm_rbh[,-2]
#ce_lm_rbh=ce_lm_rbh[,-1]
#ce_lm_rbh=unique(ce_lm_rbh)
deg=merge(deg,ce_lm_rbh,by.x = "Geneid",by.y = "lm_gene",all.x = T)
#View(deg)
write.csv(deg, "deg_anno.csv")
#write.csv(deg,"deg_anno.csv")

deg_up <- subset(deg,log2FoldChange>0)
deg_down <- subset(deg,log2FoldChange<0)

#load("D:/diske/lm_ph_fl/04deseq/enrichment/enrich_dat.RData")
#rm(ph46)
library(dplyr)
library(clusterProfiler)
library(stringr)

up_kegg = enricher(deg_up$Geneid, TERM2GENE=keggid2gene, TERM2NAME=pathway2name,pvalueCutoff = 1,pAdjustMethod = "BH",qvalueCutoff = 1)
down_kegg = enricher(deg_down$Geneid, TERM2GENE=keggid2gene, TERM2NAME=pathway2name,pvalueCutoff = 1,pAdjustMethod = "BH",qvalueCutoff = 1)

up_go = enricher(deg_up$Geneid, TERM2GENE=goid2gene, TERM2NAME=term2name,pvalueCutoff = 1,pAdjustMethod = "BH",qvalueCutoff = 1)
down_go = enricher(deg_down$Geneid, TERM2GENE=goid2gene, TERM2NAME=term2name,pvalueCutoff = 1,pAdjustMethod = "BH",qvalueCutoff = 1)

write.csv(down_kegg@result,"down_kegg.csv")
write.csv(down_go@result,"down_go.csv")
write.csv(up_kegg@result,"up_kegg.csv")
write.csv(up_go@result,"up_go.csv")

#复制到剪切板
library(rio)
export(down_kegg@result,"clipboard")

#plot
library(enrichplot)

# p cutoff = 0.05
up_kegg = enricher(deg_up$Geneid, TERM2GENE=keggid2gene, TERM2NAME=pathway2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 1)
down_kegg = enricher(deg_down$Geneid, TERM2GENE=keggid2gene, TERM2NAME=pathway2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 1)

up_go = enricher(deg_up$Geneid, TERM2GENE=goid2gene, TERM2NAME=term2name,pvalueCutoff = 1,pAdjustMethod = "BH",qvalueCutoff = 1)
down_go = enricher(deg_down$Geneid, TERM2GENE=goid2gene, TERM2NAME=term2name,pvalueCutoff = 1,pAdjustMethod = "BH",qvalueCutoff = 1)
#add go class
go_class<- read_excel("pcz-GO.xlsx")
up_go2 <- as.data.frame(up_go) %>% merge(go_class,by='ID',all.x=T, sort = F)
down_go2 <- as.data.frame(down_go) %>% merge(go_class,by='ID',all.x=T, sort = F)
#write.csv(up_go2,"up_go2.csv")
#write.csv(down_go2,"down_go2.csv")

#复制到剪切板
library(rio)
export(down_go2,"clipboard")


################
p.kegg_up<-dotplot(up_kegg,  color = "p.adjust", showCategory = 15)
p.kegg_dowm <- dotplot(down_kegg, color = "p.adjust",showCategory = 15)

p.go_up<-  dotplot(up_go2, color = "p.adjust", showCategory = 15)
p.go_down<-dotplot(down_go, color = "p.adjust",showCategory = 15)


#ggplot2 绘制，上面的方法总是缺失通路
library(ggplot2)

#将generatio转为数字
mixedToFloat <- function(x){
  x <- sapply(x, as.character)
  is.integer  <- grepl("^-?\\d+$", x)
  is.fraction <- grepl("^-?\\d+\\/\\d+$", x)
  is.float <- grepl("^-?\\d+\\.\\d+$", x)
  is.mixed    <- grepl("^-?\\d+ \\d+\\/\\d+$", x)
  stopifnot(all(is.integer | is.fraction | is.float | is.mixed))
  
  numbers <- strsplit(x, "[ /]")
  
  ifelse(is.integer,  as.numeric(sapply(numbers, `[`, 1)),
         ifelse(is.float,    as.numeric(sapply(numbers, `[`, 1)),
                ifelse(is.fraction, as.numeric(sapply(numbers, `[`, 1)) /
                         as.numeric(sapply(numbers, `[`, 2)),
                       as.numeric(sapply(numbers, `[`, 1)) +
                         as.numeric(sapply(numbers, `[`, 2)) /
                         as.numeric(sapply(numbers, `[`, 3)))))
}

##########################################################################
up.go<-read_excel('deseq2_results_0.3785.xlsx',sheet = 'up.go.plot')
up.go$GeneRatio=mixedToFloat(up.go$GeneRatio) 
#把GO的description首字母大写
library(R.utils)
up.go$Description<-capitalize(up.go$Description)

"#FFF5EB" "#FEE6CE" "#FDD0A2" "#FDAE6B" "#FD8D3C" "#F16913" "#D94801" "#A63603" "#7F2704"

"#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6" "#4292C6" "#2171B5" "#08519C" "#08306B"
#options(digits=1)
p.go.up <-ggplot(up.go,aes(x = GeneRatio, 
              y = reorder(Description,-p.adjust),# 按照富集度大小排序
              size = Count,
              colour=p.adjust,
              shape = go_ontology)) +
  geom_point() +
  scale_shape_manual(values=c('BP'=16, 'CC'=15, 'MF' = 17))+# 设置点的形状
  labs(x = "GeneRatio", y = "",color='black')+           # 设置x，y轴的名称
  scale_colour_continuous(                    # 设置颜色图例
    name="P-Value",                        # 图例名称
    low="#F16913",                              # 设置颜色范围
    high="#6BAED6")+
  scale_radius(                               # 设置点大小图例
    range=c(2,4),                             # 设置点大小的范围
    name="Count")+                             # 图例名称
  guides(   
    color = guide_colorbar(order = 1),        # 决定图例的位置顺序
    size = guide_legend(order = 2)
  )+
  theme_bw()+                                  # 设置主题
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40) )+  #调整label 长度
  theme(axis.text = element_text(color="black"),
        axis.text.y = element_text(size=10))
p.go.up
#go-down
down.go<-read_excel('deseq2_results_0.3785.xlsx',sheet = 'down.go.plot')
down.go$GeneRatio=mixedToFloat(down.go$GeneRatio)
#把GO的description首字母大写
library(R.utils)
down.go$Description<-capitalize(down.go$Description)

p.go.down <-ggplot(down.go,aes(x = GeneRatio, 
                           y = reorder(Description,-p.adjust),# 按照富集度大小排序
                           size = Count,
                           colour=p.adjust,
                           shape = go_ontology)) +
  geom_point() +
  scale_shape_manual(values=c('BP'=16, 'CC'=15, 'MF' = 17))+# 设置点的形状
  labs(x = "GeneRatio", y = "",color='black')+           # 设置x，y轴的名称
  scale_colour_continuous(                    # 设置颜色图例
    name="P-Value",                        # 图例名称
    low="#F16913",                              # 设置颜色范围
    high="#6BAED6")+
  scale_radius(                               # 设置点大小图例
    range=c(2,4),                             # 设置点大小的范围
    name="Count")+                             # 图例名称
  guides(   
    color = guide_colorbar(order = 1),        # 决定图例的位置顺序
    size = guide_legend(order = 2)
  )+
  theme_bw()+                                  # 设置主题
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40) )+  #调整label 长度
  theme(axis.text = element_text(color="black"),
        axis.text.y = element_text(size=10))
p.go.down

########
#kegg-plot
#up
up_kegg1<-as.data.frame(up_kegg)
up_kegg1$GeneRatio=mixedToFloat(up_kegg1$GeneRatio)
up_kegg1 <- up_kegg1[order(up_kegg1$pvalue,decreasing=F),] %>% subset(pvalue<0.05)

#up_kegg2 <- up_kegg1[1:15,]

p.kegg_up <-ggplot(up_kegg1,aes(x = GeneRatio, 
                            y = reorder(Description,-p.adjust),# 按照富集度大小排序
                            size = Count,
                            colour=p.adjust)) +
  geom_point(shape = 16) +                    # 设置点的形状
  labs(x = "GeneRatio", y = "",color='black')+           # 设置x，y轴的名称
  scale_colour_continuous(                    # 设置颜色图例
    name="P-Value",                        # 图例名称
    low="#F16913",                              # 设置颜色范围
    high="#6BAED6")+
  scale_radius(                               # 设置点大小图例
    range=c(2,4),                             # 设置点大小的范围
    name="Count")+                             # 图例名称
  guides(   
    color = guide_colorbar(order = 1),        # 决定图例的位置顺序
    size = guide_legend(order = 2)
  )+
  theme_bw()+                                  # 设置主题
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40) )+  #调整label 长度
  theme(axis.text = element_text(color="black"),
        axis.text.y = element_text(size=10))
p.kegg_up
#kegg-down
down_kegg1<-as.data.frame(down_kegg)
down_kegg1$GeneRatio=mixedToFloat(down_kegg1$GeneRatio)
down_kegg1 <- down_kegg1[order(down_kegg1$pvalue,decreasing=F),] %>% subset(pvalue<0.05)
#down_kegg2 <- down_kegg1[1:4,]

p.kegg_down<-ggplot(down_kegg1,aes(x = GeneRatio, 
                               y = reorder(Description,-p.adjust),# 按照富集度大小排序
                               size = Count,
                               colour=p.adjust)) +
  geom_point(shape = 16) +                    # 设置点的形状
  labs(x = "GeneRatio", y = "",color='black')+           # 设置x，y轴的名称
  scale_colour_continuous(                    # 设置颜色图例
    name="P-Value",                        # 图例名称
    low="#F16913",                              # 设置颜色范围
    high="#6BAED6")+
  scale_radius(                               # 设置点大小图例
    range=c(2,4),                             # 设置点大小的范围
    name="Count")+                             # 图例名称
  guides(   
    color = guide_colorbar(order = 1),        # 决定图例的位置顺序
    size = guide_legend(order = 2)
  )+
  theme_bw()+                                  # 设置主题
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40) )+  #调整label 长度
  theme(axis.text = element_text(color="black"),
        axis.text.y = element_text(size=10))
p.kegg_down
library(ggpubr)
p.kegg<-ggarrange(p.kegg_up,p.kegg_down,labels = c('A','B'))
p.go<-ggarrange(p.go.up,p.go.down,labels = c('A','B'))


p.kegg <- p.kegg_up+p.kegg_down
p.go <- p.go.up + p.go.down

ggsave('p.kegg.pdf',p.kegg,width = 12,height = 5)
ggsave('p.go.pdf',p.go,width = 13,height = 10)


upsetplot(up_kegg)



#data export
library(dplyr)
library(readxl)
library(rio)
data<-read_excel('Table S1.xlsx',sheet = 1)
deg<-read_excel('Table S1.xlsx',sheet = 2)
#up kegg
cel01040<-read_excel('Table S2.xlsx',sheet = 'cel01040') %>% merge(data,by.x='Gene',by.y='Geneid') %>% merge(deg,by.x='Gene',by.y='Geneid')
cel01212<-read_excel('Table S2.xlsx',sheet = 'cel01212') %>% merge(data,by.x='Gene',by.y='Geneid') %>% merge(deg,by.x='Gene',by.y='Geneid')
cel00062<-read_excel('Table S2.xlsx',sheet = 'cel00062') %>% merge(data,by.x='Gene',by.y='Geneid') %>% merge(deg,by.x='Gene',by.y='Geneid')
cel00260<-read_excel('Table S2.xlsx',sheet = 'cel00260') %>% merge(data,by.x='Gene',by.y='Geneid') %>% merge(deg,by.x='Gene',by.y='Geneid')
cel04146<-read_excel('Table S2.xlsx',sheet = 'cel04146') %>% merge(data,by.x='Gene',by.y='Geneid') %>% merge(deg,by.x='Gene',by.y='Geneid')

write.csv(cel04146,'./up.kegg.table/cel04146.csv')
#down kegg
cel04142<-read_excel('Table S3.xlsx',sheet = 'cel04142') %>% merge(data,by.x='Gene',by.y='Geneid') %>% merge(deg,by.x='Gene',by.y='Geneid')
cel00910<-read_excel('Table S3.xlsx',sheet = 'cel00910') %>% merge(data,by.x='Gene',by.y='Geneid') %>% merge(deg,by.x='Gene',by.y='Geneid')
write.csv(cel00910,'./down.kegg.table/cel00910.csv')

#go 挨个复制查找
GO_0005615<-read_excel('Table S5.xlsx',sheet = 'GO_0005615') %>% merge(data,by.x='Gene',by.y='Geneid') %>% merge(deg,by.x='Gene',by.y='Geneid')
export(GO_0005615,"clipboard")

