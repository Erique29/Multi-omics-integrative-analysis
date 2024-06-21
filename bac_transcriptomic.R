
library(ggalt)
library(MicrobiotaProcess)
library(patchwork)
library(vegan)
library(ggplot2)
library(ggforce)
library(readxl)
library(DESeq2)

#degs
readcount<-read_excel('bac_deg.xlsx') 
readcount1<-readcount[,c(-2:-7)]

readcount1<-readcount1[,-1]
rownames(readcount1)<-readcount$genes

metaData <- read_excel('bac_deg.xlsx', sheet = 'meta') 

dds <- DESeqDataSetFromMatrix(countData=readcount1, 
                              colData=metaData, 
                              design=~bacteria)
dds <- DESeq(dds#,test = 'Wald', fitType = 'mean'
)

res <- results(dds)
#res2<-subset(res, padj<.05 & abs(log2FoldChange)>1)
#res3 = cbind(as(res2, "data.frame"))

res2<-subset(res, padj<.05 & abs(log2FoldChange)> 0.585)
res3 = cbind(as(res2, "data.frame"))


#write.csv(res2, file = "deseq2_results_0.585.csv")


#vol
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="", xlim=c(-9,12)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>0.585), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#PCA
vsdata <- vst(dds, blind=FALSE)
p.pca<-plotPCA(vsdata, intgroup="bacteria")+
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
    linewidth=0.2)+
  coord_fixed(ratio = 10)
p.pca
ggsave('p.pca.pdf',p.pca,height=8,width=8)
#vol
library(EnhancedVolcano)
keyvals <- ifelse(
  res$log2FoldChange < -0.585, '#1BB6AFFF',
  ifelse(res$log2FoldChange > 0.585, '#EE6100FF',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#EE6100FF'] <- 'Up'
names(keyvals)[keyvals == 'black'] <- 'NS'
names(keyvals)[keyvals == '#1BB6AFFF'] <- 'Down'


selectLab1 <- c("ybhH", "adhB", "LBANOLMJ_00976", "ybhI", "LBANOLMJ_02962", "LBANOLMJ_04193", "nupX", "LBANOLMJ_03439", "dnaJ", "paaF", "yhaI", "gspG", "ndhI_2", "cusC_2", "tonB", "ynbA", "gspH", "atoB", "wzb", "fadK", "LBANOLMJ_02109", "pgaC", "gspE", "fliE", "kbaY", "fdnG_1", "gspI", "pilQ", "LBANOLMJ_03770", "atoD", "wcaJ", "fdnH_2", "ynbC", "glxR", "sad", "gspD", "yhcH_1", "LBANOLMJ_02966", "atoE", "gspF", "mngR_1", "fdnG_2", "gspJ", "argO", "glxK", "cyoE", "flgE", "rmf", "bolA", "glgA", "motB", "rpoS", "LBANOLMJ_01498", "glpK", "LBANOLMJ_02523", "aceB", "malP_2", "glpD", "malQ", "glcB", "motA", "glgC", "ahpC", "glgB", "ridA", "aldB", "aer", "flgD", "aroF" )
p.vol <- EnhancedVolcano(res,
                         lab = rownames(res),
                         x = 'log2FoldChange',
                         y = 'padj',
                         colCustom = keyvals,
                         #selectLab = rownames(res)[which(names(keyvals) %in% c('Up', 'Down'))], 
                         selectLab = rownames(res)[which(rownames(res) %in% selectLab1)],
                         title = "",
                         subtitle = '',
                         caption ='',
                         pointSize = 2,
                         labSize = 4,
                         pCutoff = 0.05,
                         FCcutoff = 0.585,
                         gridlines.major = F,
                         gridlines.minor = FALSE,
                         border='full',
                        # xlim = c(-9,12),
                         xlab = bquote(~Log[2] ~ "fold change" ~ "("~italic('S. algae') ~ '/' ~ italic('E. coli')~")"),
                         ylab = bquote(~-Log[10] ~ italic(P) ~ Value),
                         labCol = "black",
                         drawConnectors = T,
                         widthConnectors = 0.01,
                        legendLabels = F
)
p.vol

ggsave('DEG_vol2.pdf',p.vol,height=7.2,width=7)

###vol2
vol_data <- as.data.frame(res) %>% 
  mutate(padj = ifelse(padj == 0, 3.94318123600563E-306*0.1, padj)) %>% 
  mutate(p2=-log10(padj)) %>% mutate(name= row.names(vol_data)) 

cut_off_P =0.05
cut_off_log2FC =0.585

vol_data$Sig = ifelse(vol_data$padj < cut_off_P &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                        abs(vol_data$log2FoldChange) >= cut_off_log2FC,  #abs绝对值
                      ifelse(vol_data$log2FoldChange > cut_off_log2FC ,'Up','Down'),'NS')
vol_data$Sig <- factor(vol_data$Sig, levels = c('Up','NS','Down'))


library(ggtext) # markdown text
p_vol<-ggplot(vol_data, aes(x =log2FoldChange, y=p2, colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#EE6100FF", "#d2dae2","#1BB6AFFF")) + 
  #xlim(c(-8.5, 8.5)) +  #调整点的颜色和x轴的取值范围
  #geom_vline(xintercept=c(-cut_off_log2FC, cut_off_log2FC),lty=3,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = -log10(cut_off_P), lty=3,col="red",lwd=0.585) +  #添加y轴辅助线
  labs(x= 'Log~2~ fold enrichment (*S. algae* / *E. coli*)', y="--log~10~(*P*value)") +  #x、y轴标签
  #ggtitle("单组火山图") + #标题
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        #legend.position="right", 
        legend.title = element_blank(),
        axis.title = element_markdown(size = 11),
        #legend.justification=c(0,1), # 这个参数设置很关键
        #legend.position =c(0,1),
        axis.text = element_text(color = 'black')
        
  )+
  geom_text_repel(data = vol_data %>% 
                    filter(abs(log2FoldChange) >=6), 
                  aes(label = name, x = log2FoldChange, y = p2), 
                  color = 'black', box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines"),hjust= 0.50,force=T)
p_vol
ggsave('DEG_vol.pdf',p_vol,height=8,width=10)




#####=============enrichment   kobas=============================
library(stringr)  #调整label长度
up_data<- read_excel('bac_deg.xlsx',sheet = 'up_plot')
up_kegg<- filter(up_data, Database == 'KEGG PATHWAY')  
up_biocyc<-filter(up_data, Database == 'BioCyc') 


#把up_biocyc的description首字母大写
library(R.utils)
up_biocyc$Term<-capitalize(up_biocyc$Term)


#options(digits=1)
#up kegg
p.up_kegg <-ggplot(up_kegg,aes(x = GeneRatio, 
                           y = reorder(Term,-padj),# 按照富集度大小排序
                           size = `Input number`,
                           colour=padj,
                           #shape = go_ontology
                           )) +
  geom_point() +
  #scale_shape_manual(values=c('BP'=16, 'CC'=15, 'MF' = 17))+# 设置点的形状
  labs(x = "GeneRatio", y = "",color='black')+           # 设置x，y轴的名称
  scale_colour_continuous(                    # 设置颜色图例
    name="P-Value",                        # 图例名称
    low="#f1766d",                              # 设置颜色范围
    high="#a8cae8")+
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
p.up_kegg

#up biocyc
p.up_biocyc <-ggplot(up_biocyc,aes(x = GeneRatio, 
                               y = reorder(Term,-padj),# 按照富集度大小排序
                               size = `Input number`,
                               colour=padj,
                               #shape = go_ontology
)) +
  geom_point() +
  #scale_shape_manual(values=c('BP'=16, 'CC'=15, 'MF' = 17))+# 设置点的形状
  labs(x = "GeneRatio", y = "",color='black')+           # 设置x，y轴的名称
  scale_colour_continuous(                    # 设置颜色图例
    name="P-Value",                        # 图例名称
    low="#f1766d",                              # 设置颜色范围
    high="#a8cae8")+
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
p.up_biocyc

##down
down_data<- read_excel('bac_deg.xlsx',sheet = 'down_plot')
down_kegg<- filter(down_data, Database == 'KEGG PATHWAY')  
down_biocyc<-filter(down_data, Database == 'BioCyc') 


#把down_biocyc的description首字母大写
library(R.utils)
down_biocyc$Term<-capitalize(down_biocyc$Term)


#options(digits=1)
#down kegg
p.down_kegg <-ggplot(down_kegg,aes(x = GeneRatio, 
                               y = reorder(Term,-padj),# 按照富集度大小排序
                               size = `Input number`,
                               colour=padj,
                               #shape = go_ontology
)) +
  geom_point() +
  #scale_shape_manual(values=c('BP'=16, 'CC'=15, 'MF' = 17))+# 设置点的形状
  labs(x = "GeneRatio", y = "",color='black')+           # 设置x，y轴的名称
  scale_colour_continuous(                    # 设置颜色图例
    name="P-Value",                        # 图例名称
    low="#f1766d",                              # 设置颜色范围
    high="#a8cae8")+
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
p.down_kegg

#down biocyc
p.down_biocyc <-ggplot(down_biocyc,aes(x = GeneRatio, 
                                   y = reorder(Term,-padj),# 按照富集度大小排序
                                   size = `Input number`,
                                   colour=padj,
                                   #shape = go_ontology
)) +
  geom_point() +
  #scale_shape_manual(values=c('BP'=16, 'CC'=15, 'MF' = 17))+# 设置点的形状
  labs(x = "GeneRatio", y = "",color='black')+           # 设置x，y轴的名称
  scale_colour_continuous(                    # 设置颜色图例
    name="P-Value",                        # 图例名称
    low="#f1766d",                              # 设置颜色范围
    high="#a8cae8")+
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
p.down_biocyc

p.kegg <- p.up_kegg+p.down_kegg
p.biocyc <- p.up_biocyc + p.down_biocyc
ggsave('p.kegg.pdf',p.kegg,width = 12,height = 5)
ggsave('p.biocyc.pdf',p.biocyc,width = 12,height = 8)



#data export
library(dplyr)
library(readxl)
library(rio)
setwd('../bac_transcrip/')
data<-read_excel('Table S6.xlsx',sheet = 1)

#up kegg 挨个复制查找
kegg<-read_excel('Table S12.xlsx',sheet = 'PWY0-381') %>% merge(data,by.x='Gene',by.y='genes',all.x = T) %>% 
export("clipboard")

#diff plot
setwd('diff_enrich/')
biocyc_data<-read_excel('top_enrichment.xlsx',sheet = 'plot_biocyc')
kegg_data<-read_excel('top_enrichment.xlsx',sheet = 'plot_kegg')

library(ComplexHeatmap)
library(tidyverse)
kegg2 <- kegg_data[,c('S. algae', 'E. coli')]
rownames(kegg2) <- kegg_data$Term
biocyc2 <- biocyc_data[,c('S. algae', 'E. coli')]
rownames(biocyc2) <- biocyc_data$Term

Heatmap(kegg2,col = c('white', '#FF0000'))
Heatmap(biocyc2,col = c('white', '#FF0000'))
        

######enrich new  eggnog 重新注释 =====================
library(stringr)  #调整label长度
library(ggtext)  #markdown
# top 10% different genes
setwd('d:/3.bacs_omics/bac_transcrip/diff_enrich/')
kegg_ecoli_t10<- read_excel('diff_deg.xlsx',sheet = 'kegg_ecoli')
kegg_salgae_t10<- read_excel('diff_deg.xlsx',sheet = 'kegg_salgae')

#plot kegg e coli

p.kegg.ecoli <- ggplot(data=kegg_ecoli_t10, aes(x=reorder(`Term Name`,-`p-value`), y=GeneHitsInSelectedSet , fill = 'orange')) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  #scale_fill_manual(values = 'orange') + 
  theme_test() + 
    xlab("KEGG term") + ylab('Gene number')+
  theme(axis.text=element_text(color="black")) +
  theme(legend.position = "none")+
  ggtitle('*E. coli*')+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_markdown(),
        axis.text.y = element_text(size=10)) 
  
p.kegg.ecoli

#plot kegg s algae

p.kegg.salgae <- ggplot(data=kegg_salgae_t10, aes(x=reorder(`Term Name`,-`p-value`), y=GeneHitsInSelectedSet , fill = 'orange')) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  #scale_fill_manual(values = 'orange') + 
  theme_test() + 
  xlab("KEGG term") + ylab('Gene number')+
  theme(axis.text=element_text(color="black")) +
  theme(legend.position = "none")+
  ggtitle('*S. algae*')+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_markdown(),
        axis.text.y = element_text(size=10)) 

p.kegg.salgae

library(ggpubr)
p.kegg_diff<- ggarrange(p.kegg.ecoli, p.kegg.salgae,labels = c('A','B'))
p.kegg_diff

ggsave("diff_kegg.pdf", p.kegg_diff, width = 10, height = 6)


#GO enrichment
#go_class<- read_excel('pcz-GO.xlsx')
#go_class<- read_excel('amigo.xlsx')
go_ecoli_t10<- read_excel('diff_deg.xlsx',sheet = 'go_ecoli_plot')
go_salgae_t10<- read_excel('diff_deg.xlsx',sheet = 'go_salgae_plot')
#go_ecoli_t10$Class<- factor(go_ecoli_t10$Class, levels = c('BP', 'CC', 'MF'))
#PLOT GO
#把首字母大写
library(R.utils)
go_ecoli_t10$GO_Name<-capitalize(go_ecoli_t10$GO_Name)
go_salgae_t10$GO_Name<-capitalize(go_salgae_t10$GO_Name)

library(RColorBrewer)
color <- c('BP' = "#66C3A5",'CC' = "#8DA1CB",'MF' = "#FD8D62")
  #brewer.pal(3,"Dark2")
colorl <- rep(color)
p.go.ecoli <- ggplot(go_ecoli_t10) +
  aes(x = reorder(GO_Name,order), y = HitsGenesCountsInSelectedSet, fill = Class) +
  geom_bar(stat = "identity") +
  #scale_fill_hue() +
  coord_flip()+
  scale_fill_manual(values =color)+
  ylab('Gene number')+ xlab('GO term')+ggtitle('*E. coli*')+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_markdown(),
        )+
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    #axis.text = element_text(size=12,face="plain",color="black"),
    axis.text.y = element_text(colour = 'black',hjust=1,vjust=0.6,size=10),
    #axis.title.y = element_blank(),
    #legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    legend.direction = "vertical",
    legend.position = c(0.9,0.92),
    legend.background = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = "black"),
    
    plot.background = element_blank()
  )+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50) ) #调整label 长度
p.go.ecoli

p.go.salgae <- ggplot(go_salgae_t10) +
  aes(x = reorder(GO_Name,order), y = HitsGenesCountsInSelectedSet, fill = Class) +
  geom_bar(stat = "identity") +
  #scale_fill_hue() +
  coord_flip()+
  scale_fill_manual(values =color)+
  ylab('Gene number')+ xlab('GO term')+ggtitle('*S. algae*')+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_markdown())+
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    #axis.text = element_text(size=12,face="plain",color="black"),
    axis.text.y = element_text(colour = 'black',hjust=1,vjust=0.6,size=10),
    #axis.title.y = element_blank(),
    #legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    legend.direction = "vertical",
    legend.position = c(0.9,0.92),
    legend.background = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = "black"),
   
    plot.background = element_blank()
  )+
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 40) )+ #调整label 长度
  scale_y_continuous(limits = c(0, 10), breaks = scales::pretty_breaks()) #坐标轴整数


p.go.salgae


p.go_diff<- ggarrange(p.go.salgae, p.go.ecoli, labels = c('A','B'),widths = c(1,1))
p.go_diff

ggsave("diff_go.pdf", p.go_diff, width = 12, height = 9)

#common=================================
setwd('d:/3.bacs_omics/bac_transcrip/bac_trans_common_new/')

library(stringr)  #调整label长度
library(ggtext)  #markdown

#top 10
up_kegg<- read_excel('common.deg.xlsx',sheet = 'up.kegg') %>% arrange('p-value') %>% mutate(GeneRatio=GeneHitsInSelectedSet/AllGenesInSelectedSet)
down_kegg<- read_excel('common.deg.xlsx',sheet = 'down.kegg') %>% arrange('p-value') %>% mutate(GeneRatio=GeneHitsInSelectedSet/AllGenesInSelectedSet)
#增加gene counts

#top 10
up_kegg<- up_kegg[1:10,]
down_kegg<-down_kegg[1:10,]

#kegg-plot
#up
p.kegg_up <-ggplot(up_kegg,aes(x = GeneRatio, 
                               #y = `Term Name`,
                                y = reorder(`Term Name`,-`p-value`),# 按照富集度大小排序
                                size = GeneHitsInSelectedSet,
                                colour=`p-value`)) +
  geom_point(shape = 16) +                    # 设置点的形状
  labs(x = "GeneRatio", y = "",color='black')+           # 设置x，y轴的名称
  scale_colour_continuous(                    # 设置颜色图例
    name="P-Value",                        # 图例名称
    low="#db6968",                              # 设置颜色范围
    high="#0074b3")+
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

p.kegg_down <-ggplot(down_kegg,aes(x = GeneRatio, 
                               #y = `Term Name`,
                               y = reorder(`Term Name`,-`p-value`),# 按照富集度大小排序
                               size = GeneHitsInSelectedSet,
                               colour=`p-value`)) +
  geom_point(shape = 16) +                    # 设置点的形状
  labs(x = "GeneRatio", y = "",color='black')+           # 设置x，y轴的名称
  scale_colour_continuous(                    # 设置颜色图例
    name="P-Value",                        # 图例名称
    low="#db6968",                              # 设置颜色范围
    high="#0074b3")+
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

p.kegg

#GO
up.go<-read_excel('common.deg.xlsx',sheet = 'up.go.plot') %>% arrange('P_value') %>% mutate(GeneRatio=HitsGenesCountsInSelectedSet/AllGenesCountsInSelectedSet)
down.go<-read_excel('common.deg.xlsx',sheet = 'down.go.plot') %>% arrange('P_value') %>% mutate(GeneRatio=HitsGenesCountsInSelectedSet/AllGenesCountsInSelectedSet)


#把GO的description首字母大写
library(R.utils)
up.go$GO_Name<-capitalize(up.go$GO_Name)

"#FFF5EB" "#FEE6CE" "#FDD0A2" "#FDAE6B" "#FD8D3C" "#F16913" "#D94801" "#A63603" "#7F2704"

"#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6" "#4292C6" "#2171B5" "#08519C" "#08306B"
#options(digits=1)
p.go.up <-ggplot(up.go,aes(x = GeneRatio, 
                           y = reorder(GO_Name,-P_value),# 按照富集度大小排序
                           size = HitsGenesCountsInSelectedSet,
                           colour=P_value,
                           shape = Class)) +
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

#down
p.go.down <-ggplot(down.go,aes(x = GeneRatio, 
                           y = reorder(GO_Name,-P_value),# 按照富集度大小排序
                           size = HitsGenesCountsInSelectedSet,
                           colour=P_value,
                           shape = Class)) +
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
p.go<-ggarrange(p.go.up,p.go.down,labels = c('A','B'))
p.go
ggsave('p.kegg.pdf',p.kegg,width = 12,height = 5)
ggsave('p.go.pdf',p.go,width = 13,height = 10)


#################   new name kobas ---common   =======================

library(stringr)  #调整label长度
up_biocyc<- read_excel('./bac_trans_common_new/kobas_enrichment.xlsx',sheet = 'up_biocyc') %>% arrange(`P-Value`) %>% 
  mutate(BgRatio = Input_number/ Background_number)
down_biocyc<- read_excel('./bac_trans_common_new/kobas_enrichment.xlsx',sheet = 'down_biocyc') %>% arrange(`P-Value`) %>% 
  mutate(BgRatio = Input_number/ Background_number)

#top15
up_biocyc <- up_biocyc[1:15,]
down_biocyc <- down_biocyc[1:15,]


#把up_biocyc的description首字母大写
library(R.utils)
up_biocyc$`#Term`<-capitalize(up_biocyc$`#Term`)
down_biocyc$`#Term`<-capitalize(down_biocyc$`#Term`)

#options(digits=1)

#up biocyc
p.up_biocyc <-ggplot(up_biocyc,aes(x = BgRatio, 
                                   y = reorder(`#Term`,-`P-Value`),# 按照富集度大小排序
                                   size = `Input_number`,
                                   colour=`P-Value`,
                                   #shape = go_ontology
)) +
  geom_point() +
  #scale_shape_manual(values=c('BP'=16, 'CC'=15, 'MF' = 17))+# 设置点的形状
  labs(x = "BgRatio", y = "",color='black')+           # 设置x，y轴的名称
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
p.up_biocyc



#down biocyc
p.down_biocyc <-ggplot(down_biocyc,aes(x = BgRatio, 
                                       y = reorder(`#Term`,-`P-Value`),# 按照富集度大小排序
                                       size = `Input_number`,
                                       colour=`P-Value`,
                                       #shape = go_ontology
)) +
  geom_point() +
  #scale_shape_manual(values=c('BP'=16, 'CC'=15, 'MF' = 17))+# 设置点的形状
  labs(x = "BgRatio", y = "",color='black')+           # 设置x，y轴的名称
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
p.down_biocyc

library(ggpubr)
p.biocyc <- ggarrange(p.up_biocyc, p.down_biocyc, labels = c('A','B'))
p.biocyc

ggsave('./bac_trans_common_new/p.biocyc.pdf',p.biocyc,width = 12,height = 6)



#export table
setwd('d:/3.bacs_omics/bac_transcrip/bac_trans_common_new/')

up_biocyc <- read_excel('./kobas_enrichment.xlsx', sheet = 'up_biocyc') %>% arrange(`P-Value`)
down_biocyc <- read_excel('./kobas_enrichment.xlsx', sheet = 'down_biocyc') %>% arrange(`P-Value`)
meta <- read_excel('./kobas_enrichment.xlsx', sheet = 'meta')
#将重复的基因名全部分开
biocyc <- up_biocyc[1:15,]
name <- read_excel("name.xlsx")
library(tidyverse)

df = biocyc %>% separate_rows(Input,sep = "[|]")

df = merge(df,name,by.x = "Input",by.y = "id",all.x = T)

df = df %>% unique() %>% as.data.frame()

df  = df[,-1]

df =df %>% select(1,8)

df_new  <- df %>%
  group_by(`#Term`) %>%
  summarise(Input = toString(genes)) %>%
  ungroup()

df_b = biocyc %>% select(1:7)

df_final  = merge(df_b,df_new,by="#Term")

#export
export(df_final,"clipboard")
drawing <- df_final %>% arrange(`P-Value`)

term_list <- list()
term_list=strsplit(drawing$Input,split = ", ")
names(term_list)=drawing$`#Term`
expr = meta
rownames(expr) = expr$genes
export<-list()
export <-
  lapply(term_list, function(x) {
    as.data.frame(expr[x,])
  })
library(rio)
export(export,xlsx_file <- tempfile(fileext = ".xlsx"))


#export GO table
meta<-read_excel('common.deg.xlsx',sheet = 'meta')
name_s<-read_excel('common.deg.xlsx',sheet = 'trans_name_salgae')
name_e<-read_excel('common.deg.xlsx',sheet = 'trans_name_op50')
meta1<-merge(meta,name_s,sort=F, all.x = T) %>% merge(name_e,,sort=F, all.x = T)

up_go<- read_excel('common.deg.xlsx', sheet = 'up.go.plot')
down_go<- read_excel('common.deg.xlsx', sheet = 'down.go.plot')

#export

drawing <- down_go 

term_list <- list()
term_list=strsplit(drawing$GenesOfSelectedSetInGOterm,split = ",")
names(term_list)=drawing$GO_ID2

expr = meta1
rownames(expr) = expr$E.coli_name
export<-list()
export <-
  lapply(term_list, function(x) {
    as.data.frame(expr[x,])
  })
library(rio)
export(export,xlsx_file <- tempfile(fileext = ".xlsx"))

#kegg
up_kegg<- read_excel('common.deg.xlsx', sheet = 'up.kegg')
down_kegg<- read_excel('common.deg.xlsx', sheet = 'down.kegg')
up_kegg <- up_kegg[1:10,]
down_kegg <- down_kegg[1:10,]

drawing <- down_kegg 
term_list <- list()
term_list=strsplit(drawing$GeneListInSelectedSets,split = ", ")
names(term_list)=drawing$`Term Name`
    
expr = meta1
rownames(expr) = expr$E.coli_name
export<-list()
export <-
  lapply(term_list, function(x) {
    as.data.frame(expr[x,])
  })
library(rio)
export(export,xlsx_file <- tempfile(fileext = ".xlsx"))

#diff top 10% export
library(data.table)
setwd('../diff_enrich/')
salgae_go<- read_excel('diff_deg.xlsx', sheet = 'go_salgae_plot')
ecoli_go<- read_excel('diff_deg.xlsx', sheet = 'go_ecoli_plot')

meta_salgae<- fread('hqb5_fpkm_top_anno.csv') %>% as.data.frame()
meta_ecoli<- fread('op50_fpkm_top_anno.csv') %>% as.data.frame()


drawing <- ecoli_go 
term_list <- list()
term_list=strsplit(drawing$GenesOfSelectedSetInGOterm,split = ",")
names(term_list)=drawing$GO_ID2

expr = meta_ecoli
rownames(expr) = expr$genes2
export<-list()
export <-
  lapply(term_list, function(x) {
    as.data.frame(expr[x,])
  })
library(rio)
export(export,xlsx_file <- tempfile(fileext = ".xlsx"))

#kegg
salgae_kegg<- read_excel('diff_deg.xlsx', sheet = 'kegg_salgae')
ecoli_kegg<- read_excel('diff_deg.xlsx', sheet = 'kegg_ecoli')

meta_salgae<- fread('hqb5_fpkm_top_anno.csv') %>% as.data.frame()
meta_ecoli<- fread('op50_fpkm_top_anno.csv') %>% as.data.frame()


drawing <- ecoli_kegg 
term_list <- list()
term_list=strsplit(drawing$GeneListInSelectedSets,split = ", ")
names(term_list)=drawing$`Term Name`

expr = meta_ecoli
rownames(expr) = expr$genes2
export<-list()
export <-
  lapply(term_list, function(x) {
    as.data.frame(expr[x,])
  })
library(rio)
export(export,xlsx_file <- tempfile(fileext = ".xlsx"))


