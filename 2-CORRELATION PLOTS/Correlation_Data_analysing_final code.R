library(readxl)
library(ggplot2)
library(hrbrthemes)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)

##Libraries

data<-read.csv(file.choose())
data<- read_xlsx(file.choose())

#######Data import to excel##############
#data$siginificant6.5= "No significant"
data$siginificant6.5[data$logFC_6.5SS>=1 & data$adjPval_6.5SS<=0.05]="Significant in 6.5"
data$siginificant6.5[data$logFC_6.5SS<=-1 & data$adjPval_6.5SS<=0.05] ="Significant in 6.5"
#data$siginificant7.1= "No significant"
data$siginificant7.1[data$logFC_7.1SS>=1 & data$adjPval_7.1SS<=0.05]="Significant in 7.1"
data$siginificant7.1[data$logFC_7.1SS<=-1 & data$adjPval_7.1SS<=0.05] ="Significant in 7.1"
#data$siginificant6.8= "No significant"
data$siginificant6.8[data$logFC_6.8SS>=1 & data$adjPval_6.8SS<=0.05]="Significant in 6.8"
data$siginificant6.8[data$logFC_6.8SS<=-1 & data$adjPval_6.8SS<=0.05] ="Significant in 6.8"
data$common6.5and7.1[data$logFC_7.1SS>=1&data$adjPval_7.1SS<=0.05&data$logFC_6.5SS>=1&data$adjPval_6.5SS<=0.05]="Significant 6.5&7.1"
data$common6.5and7.1[data$logFC_7.1SS<=-1&data$adjPval_7.1SS<=0.05&data$logFC_6.5SS<=-1&data$adjPval_6.5SS<=0.05]="Significant 6.5&7.1"
data$common6.5and6.8[data$logFC_6.8SS>=1&data$adjPval_6.8SS<=0.05&data$logFC_6.5SS>=1&data$adjPval_6.5SS<=0.05]="Significant 6.5&6.8"
data$common6.5and6.8[data$logFC_6.8SS<=-1&data$adjPval_6.8SS<=0.05&data$logFC_6.5SS<=-1&data$adjPval_6.5SS<=0.05]="Significant 6.5&6.8"
data$common6.8and7.1[data$logFC_7.1SS>=1&data$adjPval_7.1SS<=0.05&data$logFC_6.8SS>=1&data$adjPval_6.8SS<=0.05]="Significant 6.8&7.1"
data$common6.8and7.1[data$logFC_7.1SS<=-1&data$adjPval_7.1SS<=0.05&data$logFC_6.8SS<=-1&data$adjPval_6.8SS<=0.05]="Significant 6.8&7.1"

write.csv(data, "Analysed in terms of significant.csv")
###########################################################################

cor1=cor(x = data$logFC_6.5SS, y=data$logFC_7.1SS,  method = "pearson", use = "complete.obs")
cor1
cols <- c("#0056f7", "#e14141","#808080","#00fff7", "#00ff62", "#f2ff00", "#d000ff")
#grey,blue,red, cyan,green, yellow , purple             

group3<-data$logFC_6.5SS>=1 & data$adjPval_6.5SS<=0.05
group4<-data$logFC_7.1SS>=2.5 & data$adjPval_7.1SS<=0.05
group5<-data$logFC_7.1SS>=1&data$adjPval_7.1SS<=0.05&data$logFC_6.5SS>=1&data$adjPval_6.5SS<=0.05
group6<-data$logFC_6.5SS<=-1 & data$adjPval_6.5SS<=0.05
group7<-data$logFC_7.1SS<=-2.5 & data$adjPval_7.1SS<=0.05
group8<-data$logFC_7.1SS<=-1&data$adjPval_7.1SS<=0.05&data$logFC_6.5SS<=-1&data$adjPval_6.5SS<=0.05

data$colplate="Non-significant"
data$colplate[data$logFC_6.5SS>=1 & data$adjPval_6.5SS<=0.05]="Significant 6.5"
data$colplate[data$logFC_7.1SS>=1 & data$adjPval_7.1SS<=0.05]="Significant 7.1"
data$colplate[data$logFC_7.1SS>=1&data$adjPval_7.1SS<=0.05&data$logFC_6.5SS>=1&data$adjPval_6.5SS<=0.05]="Significant both"
data$colplate[data$logFC_6.5SS<=-1 & data$adjPval_6.5SS<=0.05]="Significant 6.5"
data$colplate[data$logFC_7.1SS<=-1 & data$adjPval_7.1SS<=0.05]="Significant 7.1"
data$colplate[data$logFC_7.1SS<=-1&data$adjPval_7.1SS<=0.05&data$logFC_6.5SS<=-1&data$adjPval_6.5SS<=0.05]="Significant both"

##Study 65 71
k1<-ggplot(data, aes(data$logFC_6.5SS, data$logFC_7.1SS,
                     color=data$colplate, label=data$Gene.names)) +
  geom_point( size = 2)+theme_minimal()+xlim(-5,5)+ylim(-5,5)+  
  geom_hline(yintercept = 0, size=1)+
  geom_vline(xintercept=0)+theme_bw()+ #scale_color_manual(values = cols)+
  #geom_text(aes(label=ifelse(group3|group4|group5|group6|group7|group8,as.character(data$Gene.names),'')))+
  annotate("rect", xmin=-Inf, xmax=0, ymin=0, ymax=Inf, alpha=0.2, fill="blue")+
  annotate("rect", xmin=0, xmax=Inf, ymin=-Inf, ymax=0, alpha=0.2, fill="blue")+
  labs(title="Stiff_vs_Soft under pH:6.5 vs 7.1 (R=0.6178575)", y="Log2(Stiff/Soft_pH:7.1)", x="Log2(Stiff/Soft_pH:6.5)", caption="")+
  labs(color="Significance")+
  scale_color_manual(values = c("gray", "blue3", "brown1", "#009E73", "#009E73"))+
  geom_text_repel(aes(label=ifelse(group3|group4|group5|group6|group7|group8,as.character(data$Gene.names),'')),size = 3.5)+
  coord_fixed(ratio = 1) +theme(text = element_text(size = 20))

k1
plot_grid(k1, k2, k3, labels = c('A', 'B','C'))

###############################################################
##6.5vs6.8
cor(x = data$logFC_6.5SS, y=data$logFC_6.8SS,  method = "pearson", use = "complete.obs")

group3<-data$logFC_6.5SS>=1 & data$adjPval_6.5SS<=0.05
group4<-data$logFC_6.8SS>=2.5 & data$adjPval_6.8SS<=0.05
group5<-data$logFC_6.8SS>=1&data$adjPval_6.8SS<=0.05&data$logFC_6.5SS>=1&data$adjPval_6.5SS<=0.05
group6<-data$logFC_6.5SS<=-1 & data$adjPval_6.5SS<=0.05
group7<-data$logFC_6.8SS<=-2.5 & data$adjPval_6.8SS<=0.05
group8<-data$logFC_6.8SS<=-1&data$adjPval_6.8SS<=0.05&data$logFC_6.5SS<=-1&data$adjPval_6.5SS<=0.05

data$colplate="Non-significant"
data$colplate[data$logFC_6.5SS>=1 & data$adjPval_6.5SS<=0.05]="Significant 6.5,logFC>=1"
data$colplate[data$logFC_6.8SS>=1 & data$adjPval_6.8SS<=0.05]="Significant 6.8,logFC>=1"
data$colplate[data$logFC_6.8SS>=1&data$adjPval_6.8SS<=0.05&data$logFC_6.5SS>=1&data$adjPval_6.5SS<=0.05]="Significant 6.5_6.8,logFC>=1"
data$colplate[data$logFC_6.5SS<=-1 & data$adjPval_6.5SS<=0.05]="Significant 6.5, logFC<=-1"
data$colplate[data$logFC_6.8SS<=-1 & data$adjPval_6.8SS<=0.05]="Significant 6.8, logFC<=-1"
data$colplate[data$logFC_6.8SS<=-1&data$adjPval_6.8SS<=0.05&data$logFC_6.5SS<=-1&data$adjPval_6.5SS<=0.05]="Significant 6.5_6.8,logFC<=-1"

##Graphing
k1<-ggplot(data, aes(data$logFC_6.5SS, data$logFC_6.8SS,
                     color=data$colplate, label=data$Gene.names)) +
  geom_point( size = 2)+theme_minimal()+xlim(-5,5)+ylim(-5,5)+  
  geom_hline(yintercept = 0, size=1)+
  geom_vline(xintercept=0)+theme_bw()+ #scale_color_manual(values = cols)+
  #geom_text(aes(label=ifelse(group3|group4|group5|group6|group7|group8,as.character(data$Gene.names),'')))+
  annotate("rect", xmin=-Inf, xmax=0, ymin=0, ymax=Inf, alpha=0.2, fill="blue")+
  annotate("rect", xmin=0, xmax=Inf, ymin=-Inf, ymax=0, alpha=0.2, fill="blue")+
  labs(title="Stiff/Soft-pH:6.5_vs_6.8 (R=0.6547192)", y="Stiff/Soft_pH:6.8", x="Stiff/Soft_pH:6.5", caption="Source: HG_YHU")+
  labs(color="Significance")+
  #geom_text_repel(aes(data$logFC_6.5SS, data$logFC_6.8SS, label = data$Gene.names), size = 3)+
  scale_color_manual(values = c("gray", "blue3", "brown1", "#009E73", "#009E73"))+
  geom_text_repel(aes(label=ifelse(group3|group4|group5|group6|group7|group8,as.character(data$Gene.names),'')),size = 3.5)+
  coord_fixed(ratio = 1)
#geom_hline(yintercept=1, linetype="dashed", 
#          color = "green", size=1)+ geom_hline(yintercept=-1, linetype="dashed", 
#                                               color = "green", size=1)+
#geom_vline(xintercept=1, linetype="dashed", 
#         color = "green", size=1)+ geom_vline(xintercept=-1, linetype="dashed", 
#                                             color = "green", size=1)


k1
plot_grid(k1, k2, k3, labels = c('A', 'B','C'))

###########################################################################
#6.8 and 7.1
cor(x = data$logFC_6.8SS, y=data$logFC_7.1SS,  method = "pearson", use = "complete.obs")
      
group3<-data$logFC_6.8SS>=1 & data$adjPval_6.8SS<=0.05
group4<-data$logFC_7.1SS>=2.5 & data$adjPval_7.1SS<=0.05
group5<-data$logFC_7.1SS>=1&data$adjPval_7.1SS<=0.05&data$logFC_6.8SS>=1&data$adjPval_6.8SS<=0.05
group6<-data$logFC_6.8SS<=-1 & data$adjPval_6.8SS<=0.05
group7<-data$logFC_7.1SS<=-2.5 & data$adjPval_7.1SS<=0.05
group8<-data$logFC_7.1SS<=-1&data$adjPval_7.1SS<=0.05&data$logFC_6.8SS<=-1&data$adjPval_6.8SS<=0.05

data$colplate="Non-significant"
data$colplate[data$logFC_6.8SS>=1 & data$adjPval_6.8SS<=0.05]="Significant 6.8"
data$colplate[data$logFC_7.1SS>=1 & data$adjPval_7.1SS<=0.05]="Significant 7.1"
data$colplate[data$logFC_7.1SS>=1&data$adjPval_7.1SS<=0.05&data$logFC_6.8SS>=1&data$adjPval_6.8SS<=0.05]="Significant both"
data$colplate[data$logFC_6.8SS<=-1 & data$adjPval_6.8SS<=0.05]="Significant 6.8"
data$colplate[data$logFC_7.1SS<=-1 & data$adjPval_7.1SS<=0.05]="Significant 7.1"
data$colplate[data$logFC_7.1SS<=-1&data$adjPval_7.1SS<=0.05&data$logFC_6.8SS<=-1&data$adjPval_6.8SS<=0.05]="Significant both"


k1<-ggplot(data, aes(data$logFC_6.8SS, data$logFC_7.1SS,
                     color=data$colplate, label=data$Gene.names)) +
  geom_point( size = 2)+theme_minimal()+xlim(-5,5)+ylim(-5,5)+  
  geom_hline(yintercept = 0, size=1)+
  geom_vline(xintercept=0)+theme_bw()+ #scale_color_manual(values = cols)+
  #geom_text(aes(label=ifelse(group3|group4|group5|group6|group7|group8,as.character(data$Gene.names),'')))+
  annotate("rect", xmin=-Inf, xmax=0, ymin=0, ymax=Inf, alpha=0.2, fill="blue")+
  annotate("rect", xmin=0, xmax=Inf, ymin=-Inf, ymax=0, alpha=0.2, fill="blue")+
  labs(title="Stiff_vs_Soft under pH:6.8 vs 7.1 (R=0.9586469)", y="Log2(Stiff/Soft_pH:7.1)", x="Log2(Stiff/Soft_pH:6.8)", caption="")+
  labs(color="Significance")+
  scale_color_manual(values = c("gray", "blue3", "brown1", "#009E73", "#009E73"))+
  geom_text_repel(aes(label=ifelse(group3|group4|group5|group6|group7|group8,as.character(data$Gene.names),'')),size = 3.5)+
  coord_fixed(ratio = 1)+theme(text = element_text(size = 20)) 
        


k1
plot_grid(k1, k2, k3, labels = c('A', 'B','C'))

