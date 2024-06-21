rm(list=ls())

BiocManager::install("ComplexHeatmap ")


library(data.table)
library(readxl)
library(dplyr)
library(ggtext) # markdown text
library(gt)
library(ggrepel)
library(ggpubr)

setwd('d:/Metabolome/RESULTS/')
S.algae_E.coli<-read_excel('HQB5.vs.OP50_all.xlsx')
S.algae_LB<-read_excel('HQB5.vs.LB_all.xlsx')
E.coli_LB<-read_excel('OP50.vs.LB_all.xlsx')

S.algae_LB2 <- subset(S.algae_LB,select = c(Compound_ID, LB_1, LB_2, LB_3, S.algae_vs_LB_FC, S.algae_vs_LB_log2FC, 
                                          S.algae_vs_LB_Pvalue, S.algae_vs_LB_ROC, S.algae_vs_LB_VIP))
E.coli_LB2 <- subset(E.coli_LB, select = c(Compound_ID, E.coli_vs_LB_FC, E.coli_vs_LB_log2FC, 
                                           E.coli_vs_LB_Pvalue, E.coli_vs_LB_ROC, E.coli_vs_LB_VIP))


meta_all_1<-merge(S.algae_E.coli,S.algae_LB2,by='Compound_ID')
meta_all<-merge(meta_all_1, E.coli_LB2,by='Compound_ID')

#add average
meta_all <- mutate(meta_all,S.algae = (S.algae_1+S.algae_2+S.algae_3+S.algae_4+S.algae_5)/5)
meta_all <- mutate(meta_all,E.coli = (E.coli_1+E.coli_2+E.coli_3+E.coli_4+E.coli_5)/5)
meta_all <- mutate(meta_all, LB = (LB_1+LB_2+LB_3)/3)
#save
#openxlsx::write.xlsx(meta_all,'meta_all.xlsx')

#Pcoa
library(ggalt)
library(MicrobiotaProcess)
library(patchwork)
library(vegan)
library(ggplot2)
library(ggforce)

pcoa_data<-read_excel('meta_pcoa.xlsx',sheet = 'pcoa')
pcoa_data1<-pcoa_data[,-1]
row.names(pcoa_data1)<-pcoa_data$Compound_ID

distMatrix <- vegdist(pcoa_data1,method = "cao")  #canberra   clark  bray!  hellinger !  cao!  binomial!
pCoa <- cmdscale(distMatrix, eig = T,k = 2 )
varExp <- (eigenvals(pCoa)/sum(eigenvals(pCoa)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

pCoaVecs <- as.data.frame(pCoa$points)
colnames(pCoaVecs) <- paste0("PCo",c(1:2))
pCoaVecs$ID <- row.names(pCoaVecs)

meta<-read_excel('meta_pcoa.xlsx',sheet = 'meta')
pCoaVecs2<-merge(pCoaVecs,meta,by='ID')

#PERMANOVA adonis2
adonis_result <- adonis2(pcoa_data1 ~ group, data=meta, permutations=999)

library(paletteer)
pCoaVecs2$group <- factor(pCoaVecs2$group, levels = c('Unconditioned media','S.algae-conditioned media','E.coli-conditioned media'))

p_pcoa<-ggplot(pCoaVecs2,aes(x=PCo1,y=PCo2,color=group)) + 
  geom_point(size=4) + 
  theme_bw() + 
  scale_color_gradientn(colours = rainbow(5)) + 
  xlab(paste0('PC1 (',round(xVar,2),' %)')) + 
  ylab(paste0('PC2 (',round(yVar,2),' %)')) + 
  theme(axis.text = element_text(color = 'black'),
        legend.position="top",
        legend.title=element_blank())+
  scale_colour_manual(name = "Group", values = c('black','#1BB6AFFF', '#EE6100FF'),
                      labels=c("Unconditioned\n       media", "   E. coli-\nconditioned\n    media", "   S. algae-\nconditioned\n    media"))+
  geom_mark_ellipse(#data = pCoaVecs2, #aes(x = PCo1, y = PCo2),
                    expand = 0.02, radius = 0.01,tol = 0.01,linetype=2,
                    linewidth=0.2)+
  annotate("text",x=-0.25,y=0.2,label=paste("p= ", adonis_result$`Pr(>F)`[1]),size=4)


p_pcoa
ggsave('pcoa.pdf',p_pcoa,height=5,width=5)

#volcanol
vol_data0 <- read_excel('vol_data.xlsx')
vol_data<- subset(vol_data0, select=c(Compound_ID,Name,Name2, S.algae_vs_E.coli_log2FC,S.algae_vs_E.coli_Pvalue)) %>% 
  mutate(p2=-log10(S.algae_vs_E.coli_Pvalue))
cut_off_P =0.05
cut_off_log2FC =0.585
vol_data$Sig = ifelse(vol_data$S.algae_vs_E.coli_Pvalue < cut_off_P &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                    abs(vol_data$S.algae_vs_E.coli_log2FC) >= cut_off_log2FC,  #abs绝对值
                  ifelse(vol_data$S.algae_vs_E.coli_log2FC > cut_off_log2FC ,'Up','Down'),'NS')
vol_data$Sig <- factor(vol_data$Sig, levels = c('Up','NS','Down'))


cmplist_up<-c('Com_3322_neg', 'Com_404_pos', 'Com_721_neg', 'Com_12413_pos')
cmplist_up1 <- vol_data[vol_data$Compound_ID %in% cmplist_up,]

cmplist_down<-c('Com_1846_neg', 'Com_648_neg', 'Com_13968_neg', 'Com_3415_neg', 'Com_6072_pos', 
                'Com_5448_pos', 'Com_6676_pos', 'Com_13093_neg', 'Com_985_neg', 'Com_14297_pos', 
                'Com_14372_pos', 'Com_4931_neg', 'Com_2610_neg', 'Com_10630_neg', 'Com_2095_neg', 
                'Com_16448_neg', 'Com_499_pos', 'Com_371_neg', 'Com_1393_neg', 'Com_23231_pos')
cmplist_down1 <- vol_data[vol_data$Compound_ID %in% cmplist_down,]


p_vol<-ggplot(vol_data, aes(x =S.algae_vs_E.coli_log2FC, y=p2, colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#EE6100FF", "#d2dae2","#1BB6AFFF")) + xlim(c(-8.5, 8.5)) +  #调整点的颜色和x轴的取值范围
  #geom_vline(xintercept=c(-cut_off_log2FC, cut_off_log2FC),lty=3,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = -log10(cut_off_P), lty=3,col="red",lwd=0.5) +  #添加y轴辅助线
  labs(x= 'Log~2~ fold enrichment (*S. algae* / *E. coli*)', y="--log~10~(*P*value)") +  #x、y轴标签
  #ggtitle("单组火山图") + #标题
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        #legend.position="right", 
        legend.title = element_blank(),
        axis.title = element_markdown(size = 11),
        legend.justification=c(0,1), # 这个参数设置很关键
        legend.position =c(0,1),
        axis.text = element_text(color = 'black')
        
  )+
  geom_text_repel(data = cmplist_up1 #%>% 
                    #filter(abs(S.algae_vs_E.coli_log2FC) >=5)
                  , 
                  aes(label = Name, x = S.algae_vs_E.coli_log2FC, y = p2), 
                  color = 'black', box.padding = unit(0.7, "lines"),hjust= 0.30)+
  geom_text_repel(data = cmplist_down1 #%>% 
                  #filter(abs(S.algae_vs_E.coli_log2FC) >=5)
                  , 
                  aes(label = Name2, x = S.algae_vs_E.coli_log2FC, y = p2), 
                  color = 'black', box.padding = unit(0.7, "lines"),hjust= 0.30)
p_vol


ggsave('vol2.pdf',p_vol,width = 8,height = 6)

top_20 <- bind_rows(   #分别筛选差异显著前10个的上下调基因，并合并两组数值进行绘图
  vol_data %>%
    filter(Sig == 'Up') %>%
    arrange(desc(abs(S.algae_vs_E.coli_log2FC))) %>%
    head(10),
  vol_data %>%
    filter(Sig == 'Down') %>%
    arrange(desc(abs(S.algae_vs_E.coli_log2FC))) %>%
    head(10))

top_20 %>% gt()  #将数据制成表
#绘图——添加基因标签框图#
p_vol +
  geom_label_repel(data = top_20,
                   aes(S.algae_vs_E.coli_log2FC, p2, label = Name), size = 3, fill="#CCFFFF"
                   )
AB<-ggarrange(p_pcoa,p_vol,labels = c('A','B'))                   
AB
ggsave('Fig.meta.pdf',AB,height = 4,width = 8)

#KEGG
#enrich
library("MetaboAnalystR")
setwd('d:/3.bacs_omics/Metabolome/RESULTS/MetaboAnalyst/')
#up
mSet<-InitDataObjects("conc", "msetqea", FALSE)
mSet<-Read.TextData(mSet, "D:/3.bacs_omics/Metabolome/RESULTS/MetaboAnalyst/up.new/up_enrich_data_named.csv", "colu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "pdf", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "pdf", 72, width=NA)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "RaMP_pathway", 2);
mSet<-CalculateGlobalTestScore(mSet)
mSet<-PlotQEA.Overview(mSet, "qea_0_", "net", "pdf", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "qea", "qea_dot_0_", "pdf", 72, width=NA)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
mSet<-CalculateGlobalTestScore(mSet)
mSet<-PlotQEA.Overview(mSet, "qea_1_", "net", "pdf", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "qea", "qea_dot_1_", "pdf", 72, width=9)
mSet<-SaveTransformedData(mSet)

#down
mSet<-InitDataObjects("conc", "msetqea", FALSE)
mSet<-Read.TextData(mSet, "d:/3.bacs_omics/Metabolome/RESULTS/MetaboAnalyst/down.new/down_enrich_data_named.csv", "colu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "pdf", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "pdf", 72, width=NA)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "RaMP_pathway", 2);
mSet<-CalculateGlobalTestScore(mSet)
mSet<-PlotQEA.Overview(mSet, "qea_0_", "net", "pdf", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "qea", "qea_dot_0_", "pdf", 72, width=NA)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
mSet<-CalculateGlobalTestScore(mSet)
mSet<-PlotQEA.Overview(mSet, "qea_1_", "net", "pdf", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "qea", "qea_dot_1_", "pdf", 72, width=9)
mSet<-SaveTransformedData(mSet)

#all
mSet<-InitDataObjects("conc", "msetqea", FALSE)
mSet<-Read.TextData(mSet, "d:/3.bacs_omics/Metabolome/RESULTS/MetaboAnalyst/all.new/all_enrich_data_named.csv", "colu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "pdf", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "pdf", 72, width=NA)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "RaMP_pathway", 2);
mSet<-CalculateGlobalTestScore(mSet)
mSet<-PlotQEA.Overview(mSet, "qea_0_", "net", "pdf", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "qea", "qea_dot_0_", "pdf", 72, width=NA)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
mSet<-CalculateGlobalTestScore(mSet)
mSet<-PlotQEA.Overview(mSet, "qea_1_", "net", "pdf", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "qea", "qea_dot_1_", "pdf", 72, width=9)
mSet<-SaveTransformedData(mSet)


#metID  annotation


#class pie plot
library(readxl)
library(tidyverse)
up_data<-read_excel('./meta_all.xlsx',sheet = 'up_class2')
down_data<-read_excel('./meta_all.xlsx',sheet = 'down_class2')

library(ggplot2)
#up-class
up_data <- up_data %>% 
  arrange(desc(Class_I)) %>%
  mutate(prop = count / sum(up_data$count) *100)%>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
#两位小数
up_data$prop<- round(up_data$prop, 2)

#coloyr
library(paletteer)
paletteer_d("LaCroixColoR::paired")

col <- c('Lipids and lipid-like molecules'= '#C70E7BFF', 'Organic acids and derivatives' = '#FC6882FF', 
         'Organoheterocyclic compounds'='#007BC3FF', 'Nucleosides, nucleotides, and analogues'= '#54BCD1FF',
         'Organic oxygen compounds'='#EF7C12FF', 'Benzenoids'='#F4B95AFF','Phenylpropanoids and polyketides'='#009F3FFF',
         'Organic nitrogen compounds'='#8FDA04FF', 'Hydrocarbons' = '#AF6125FF', 'Alkaloids and derivatives'='#F4E3C7FF')
    #B25D91FF #EFC7E6FF #EF7C12FF #F4B95AFF
p.up <- ggplot(up_data, aes(x="", y=prop, fill=Class_I)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)+
  theme_void()+ 
  #theme(legend.position="none") +
  geom_text(aes(y = ypos, label = prop), color = "black", size=3) +
  scale_fill_manual(values = col, # assign legend colors
                    guide = guide_legend(reverse = TRUE)
  )
p.up

#down-class
down_data <- down_data %>% 
  arrange(desc(Class_I)) %>%
  mutate(prop = count / sum(down_data$count) *100)%>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
#两位小数
down_data$prop<- round(down_data$prop, 2)

#coloyr
library(paletteer)
paletteer_d("LaCroixColoR::paired")

col <- c('Lipids and lipid-like molecules'= '#C70E7BFF', 'Organic acids and derivatives' = '#FC6882FF', 
         'Organoheterocyclic compounds'='#007BC3FF', 'Nucleosides, nucleotides, and analogues'= '#54BCD1FF',
         'Organic oxygen compounds'='#EF7C12FF', 'Benzenoids'='#F4B95AFF','Phenylpropanoids and polyketides'='#009F3FFF',
         'Organic nitrogen compounds'='#8FDA04FF', 'Hydrocarbons' = '#AF6125FF', 'Alkaloids and derivatives'='#F4E3C7FF')
#B25D91FF #EFC7E6FF #EF7C12FF #F4B95AFF
p.down <- ggplot(down_data, aes(x="", y=prop, fill=Class_I)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)+
  theme_void()+ 
  #theme(legend.position="none") +
  geom_text(aes(y = ypos, label = prop), color = "black", size=3) +
  scale_fill_manual(values = col, # assign legend colors
                    guide = guide_legend(reverse = TRUE)
  )
p.down

library(ggpubr)
p.pie <- ggarrange(p.up, p.down, labels = c('A', 'B'), common.legend = T,legend='right'
                   )
p.pie
ggsave('pie.pdf',p.pie,  width = 10, height = 5)


#metscape
library(readxl)
meta<-read_excel('./meta_all.xlsx')
id<-read_excel('Metscape.xlsx')
merge<-merge(id, meta, by.x = 'Query', by.y = 'Name', sort = F)
library(rio)
export(merge,"clipboard")

#trans
setwd('c:/Users/Gemini/Desktop/')
trans<-read_excel('./Metscape_trans.xlsx',sheet = 2)
id<-read_excel('Metscape_trans.xlsx', sheet = 3)
merge<-merge(id, trans, by.x = 'gene_cor', by.y = 'gene', sort = F, all.x = T)




#joint
library(tidyverse)
setwd('./up_down_joint/')
input <-read_excel('./up_all_pathway/jointpa_matched_features.xlsx',sheet = 2)
entrz<-read_excel('./up_all_pathway/jointpa_matched_features.xlsx',sheet = 'entrz')
symbol<- read_excel('./up_all_pathway/jointpa_matched_features.xlsx',sheet = 'symbol')

merge<-merge(input,entrz,by.x='ID',by.y='Symbol',all.x=T, sort=F) %>% merge(symbol, by.x='Entrez',by.y='GeneID',all.x=T,sort=F)
library(rio)
export(merge,"clipboard")
