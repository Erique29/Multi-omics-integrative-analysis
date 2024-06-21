#curve plot
library(ggalt)
library(tidyverse)
library(paletteer)
library(readxl)
#curve_1
curve_1 <-read_excel('growth_data.xlsx', sheet='growthplot')


curve_1 %>% 
  rowwise() %>% 
  mutate(mean_value=mean(c(rep1,rep2,rep3)),
         std_error=plotrix::std.error(c(rep1,rep2,rep3))) %>% 
  ungroup() -> new.curve_1

new.curve_1$strain<-factor(new.curve_1$strain,levels=c('E. coli','S. algae'))


#plot--cur1
p.growth<- ggplot(data=new.curve_1,aes(x=Day,y=mean_value,group=strain, color=color#,linetype=line
))+
  geom_errorbar(aes(ymin=mean_value-std_error,
                    ymax=mean_value+std_error),linetype="solid",
                width=0.1,linewidth=0.4)+
  geom_point(size=0.8)+
  geom_line(linewidth =0.5)+
  #geom_xspline(spline_shape = -0.6,size =6)+
  #scale_linetype_manual(values=c( 'solid',"dotted"))+
  geom_vline(xintercept=5, linewidth= 0.5,color="lightgrey",linetype = "dashed")+
  scale_color_manual(values = c('#1BB6AFFF', '#EE6100FF'),
                     name=NULL)+
  scale_x_continuous(#limits = c(3,10),
    expand = expansion(mult=c(0,0)),
    breaks = seq(10)
  )+
  scale_y_continuous(limits = c(-0.5,100),
                     expand = expansion(mult=c(0,0)),
                     #breaks = seq(),
  )+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(legend.text = element_text(hjust=0),
        # legend.position = c(0.1,0.4)
        axis.text = element_text(colour = "black")
  )+
  labs(x="Day",
       y=expression(paste("L4 proportion of L. marina (%)")))
p.growth


#life curve
library(tidyverse)
library(ggsurvfit)
library(gghighlight)
library(tidycmprsk)
library(paletteer) 
lifespan_data <- read_excel("./growth_data.xlsx",sheet = 'lifespan_plot')
lifespan_data$group<-factor(lifespan_data$group,levels=c('E. coli','S. algae'))


surcical<-
  survfit2(Surv(time, status) ~ group, data = lifespan_data) %>% 
  ggsurvfit(size = 0.5)+
  add_censor_mark()+
  #add_confidence_interval() +
  #add_risktable(theme=theme_test()+
  #                theme(axis.title = element_blank(),
  #                      axis.text.x = element_blank(),
  #                      axis.ticks.x=element_blank(),
  #                      axis.text.y = element_text(color="black",size=10)),
  #              combine_groups=F)+
  add_quantile(color ="grey80",size=0.5,linetype =5)+theme_classic()+
  labs(y = "Percentage Survival",x='Day')+
  add_pvalue(caption = "Log-rank {p.value}",location  = "annotation", x = 20)+
  guides(colour =guide_legend(ncol=1))+
  scale_color_manual(values = c('#1BB6AFFF', '#EE6100FF'))+
  scale_fill_manual(values = c('#1BB6AFFF', '#EE6100FF'))+
  scale_x_continuous(breaks = seq(0, 24, by = 2))+
  theme(axis.text = element_text(color="black"))
#theme(legend.position="none")
surcical

library(ggpubr)
Fig1<- ggarrange(p.growth, surcical,common.legend = T, legend = 'top', labels = c('A','B') )
Fig1
ggsave('Fig. 1.pdf', Fig1, width = 7.5,height = 3.5)
