library(ggplot2)
library(ggpubr)
library(cowplot)
library(wesanderson)

data<-read.csv(file.choose())

SS6.5 <- data.frame(length = data$logFC_6.5SS)
SS7.1 <- data.frame(length = data$logFC_7.1SS)
SS6.8 <- data.frame(length = data$logFC_6.8SS)

SS6.5$Condition <- 'Stiff_vs_Soft:6.5'
SS7.1$Condition <- 'Stiff_vs_Soft:7.1'
SS6.8$Condition <- 'Stiff_vs_Soft:6.8'
vegLengths <- rbind(SS6.5, SS7.1,SS6.8)
grp.mean=mean(length)

ggplot(vegLengths, aes(length, fill = Condition)) + geom_density(alpha = 0.5, position="identity")+
  xlim(-1, 1)+ xlab("LogFC") + ylab("Protein Density")+
  theme_bw()+theme(text = element_text(size = 15)) #Density

gghistogram(vegLengths, x="length", fill="Condition",alpha =0.5,
            linetype = "solid",add="mean", rug =1,
            position = position_identity(), bins = 35)+xlim(-1,1)+ylim(0,500)+
  xlab("Log2FC") + ylab("Protein Numbers")+theme(text = element_text(size = 18))+
  theme(legend.position="right")

gghistogram(vegLengths, x="length", fill="Condition",alpha =0.5,
            position = position_identity(), bins = 35)+xlim(-1,1)+ylim(0,500)+
  xlab("Log2FC") + ylab("Protein Numbers")+
  theme(legend.position="right")+theme_bw()+theme(text = element_text(size =20))


##
phist<-gghistogram(vegLengths, x="length",y="..count..", fill="Condition",alpha = 0.5,bins = 35,
                   rug = 1)+
  xlim(-1,1)+ylim(0,500)+
  xlab("Log2FC") + ylab("Number of proteins")+theme(text = element_text(size = 15))+
  theme(legend.position="bottom")
  
pdensity<-ggdensity(vegLengths, x="length", color="Condition", alpha = 0, size = 1)+xlim(-1,1)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "right")  +
  theme_half_open(11, rel_small = 1) +ylab("Density of proteins")+
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend")+theme(text = element_text(size = 15))

aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")

ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

