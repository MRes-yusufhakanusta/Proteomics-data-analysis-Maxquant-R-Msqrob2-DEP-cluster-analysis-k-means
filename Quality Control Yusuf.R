options(java.parameters = "- Xmx1024m")
library(eulerr)
library(dplyr)
library(msqrob2)
library(ggplot2)
library(plotly)
library(gridExtra)
library(pheatmap)
library(ReactomePA)
library(Biobase)
library(MSnbase)
library(limma)
library(readxl)
library(reshape2)
library(cowplot)
library(ggfortify)
library(IRdisplay)
library(xlsx)
library(xlsxjars)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(htmlwidgets)
#options(warn = defaultW)
library(Matrix)
library(kableExtra)
library(ComplexHeatmap)
library(ggrepel)
library(ExploreModelMatrix)
library(enrichplot)
##Libraries
#library(rJava)
library(ggrepel)
library(tidyverse)



peptidesFile <- "C:/Users/YHU/Desktop/29.09.2022-9606/peptides.txt" 

ecols <- grep(
  "Intensity\\.",
  names(read.delim(peptidesFile))
)
ecols

pe <- readQFeatures(
  table = peptidesFile,
  fnames = 1,
  ecol = ecols,
  name = "peptideRaw", sep="\t")
pe

##Experimental Design 
colData(pe)$matrix <- gsub("\\.", "", gsub(".*_(\\w+\\.*)_\\d.*", "\\1", 
                                           colnames(pe[["peptideRaw"]]))) %>%
  unlist %>%
  as.factor

colData(pe)$pH <- gsub(".*\\.pH_(\\S+)_\\w+\\.*_\\d.*", "\\1", 
                       colnames(pe[["peptideRaw"]])) %>%
  unlist %>%
  as.factor

colData(pe)$matrix_pH_group <-  colData(pe)$matrix:colData(pe)$pH%>%
  unlist %>%
  as.factor

colData(pe)$Donor <- c("WH094","WH094","WH094","WH094","WH094","WH094","WH312",
                       "WH312","WH312","WH312","WH312","WH312","WH314","WH314",
                       "WH314","WH314","WH314","WH314","WH330","WH330")%>%
  as.factor

design<-as.data.frame(colData(pe))
design$Sample<-gsub("Intensity.", "", rownames(design))
design[, c(5, 1:4)] %>%
  kbl(row.names = FALSE) %>%
  kable_classic("hover", full_width = F, font_size = 20)

##Missing value
rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA
MSnbase::plotNA(assay(pe[["peptideRaw"]]))+
  xlab("Peptide index (ordered by data completeness)")+
  theme_linedraw()

##Pre-processing
pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
melt_assay<-melt(assay(pe[["peptideLog"]]))

ggplot(melt(assay(pe[["peptideLog"]])))+
  geom_density(aes(value, colour=Var2))+
  scale_x_continuous("Log2 Raw Peptide Intensity")+
  scale_color_discrete("")+
  theme_linedraw()

##Filtering
pe <- filterFeatures(pe,
                     ~ Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins)
)
pe <- filterFeatures(pe, ~ Potential.contaminant != "+")
pe <- filterFeatures(pe, ~ Reverse != "+")
pe <- filterFeatures(pe, ~ nNonZero >= 2)
nrow(pe[["peptideLog"]])

##Normalizing Data
pe <- normalize(pe, i = "peptideLog", method = "quantiles", name = "peptideNorm")

ggplot(melt(assay(pe[["peptideNorm"]])))+
  geom_density(aes(value, colour=Var2))+
  scale_x_continuous("Log2 Normalised Peptide Intensity")+
  scale_color_discrete("")+
  theme_linedraw()

ggplot(melt(assay(pe[["peptideNorm"]])))+
  geom_boxplot(aes(value, fill=Var2))+
  scale_x_continuous("Log2 Normalised Peptide Intensity")+
  scale_color_discrete("")+
  theme_linedraw()

##Multiscale drawing
mdsplot<-limma::plotMDS(assay(pe[["peptideNorm"]]),
                        col = colData(pe)$matrix:colData(pe)$pH %>%
                          as.numeric,
                        labels = gsub("Intensity.", "",rownames(colData(pe))), plot=FALSE)

mds_data<-data.frame(Dim1=mdsplot$x,Dim2=mdsplot$y,
                     matrix=design$matrix, pH=design$pH, matrix_ph=design$matrix_pH_group, donor=design$Donor,
                     sample=gsub("Intensity.", "", design$Sample))

ggplot(mds_data)+
  geom_text(aes(Dim1, Dim2, colour=donor, label=sample))+
  scale_x_continuous(paste("Dim1 ", round(100*mdsplot$var.explained[1], 2), "%", sep=""), limits=c(-2.5,3))+
  scale_y_continuous(paste("Dim1 ", round(100*mdsplot$var.explained[2], 2), "%", sep=""), limits=c(-1.5,2.5))+
  theme_linedraw()


##Summarize the protein level
pe <- aggregateFeatures(pe,
                        i = "peptideNorm",
                        fcol = "Proteins",
                        na.rm = TRUE,
                        name = "proteinRobust",
                        fun = MsCoreUtils::robustSummary)


##Multiscale summarised plot
mdsplot<-limma::plotMDS(assay(pe[["proteinRobust"]]),
                        col = colData(pe)$matrix:colData(pe)$pH %>%
                          as.numeric,
                        labels = gsub("Intensity.", "",rownames(colData(pe))), plot=FALSE)

mds_data<-data.frame(Dim1=mdsplot$x,Dim2=mdsplot$y,
                     matrix=design$matrix, pH=design$pH, matrix_ph=design$matrix_pH_group, donor=design$Donor, 
                     sample=gsub("Intensity.", "", design$Sample))

ggplot(mds_data)+
  geom_text(aes(Dim1, Dim2, colour=donor, label=sample))+
  scale_x_continuous(paste("Dim1 ", round(100*mdsplot$var.explained[1], 2), "%", sep=""), limits=c(-2.5,2.5))+
  scale_y_continuous(paste("Dim1 ", round(100*mdsplot$var.explained[2], 2), "%", sep=""), limits=c(-1.5,2.5))+
  theme_linedraw()

##Deleting unwanted data
ph7.1_Soft<-as.data.frame(assay(pe[["proteinRobust"]])[,grep("7.1_Soft", colnames(assay(pe[["proteinRobust"]])))])
colnames(ph7.1_Soft)<-gsub("Intensity.", "", colnames(ph7.1_Soft))

for(i in 1:ncol(ph7.1_Soft)){ph7.1_Soft[,i]<-ifelse(is.na(ph7.1_Soft[,i]), FALSE, TRUE)}
eufit <- euler(as.matrix(ph7.1_Soft), shape = "ellipse")
plot(eufit,
     quantities = TRUE,
     # fill = "transparent",
     lty = 1:4,
     labels = list(font = 4))

##Data Analysis
pe<-msqrob(object=pe, i="proteinRobust", 
           formula= ~matrix*pH+(1|Donor), ridge=TRUE, overwrite=TRUE)

VisualizeDesign(colData(pe), designFormula = ~matrix*pH)$plotlist

##Make a contrast

L2 <-makeContrast(
  c(#Stiff/Soft @pH6.5
    "ridgematrixStiff=0",
    
    #Stiff/Soft @pH6.8
    "ridgematrixStiff+ridgematrixStiff:pH6.8=0",
    
    #Stiff/Soft @pH7.1
    "ridgematrixStiff+ridgematrixStiff:pH7.1=0",
    
    #pH6.8/ph6.5 in soft
    "ridgepH6.8=0",
    
    #pH6.8/ph6.5 in stiff
    "ridgepH6.8+ridgematrixStiff:pH6.8=0",
    
    #pH7.1/ph6.5 in soft
    "ridgepH7.1=0",
    
    #pH7.1/ph6.5 in stiff
    "ridgepH7.1+ridgematrixStiff:pH7.1=0",
    
    #pH7.1/ph6.8 in soft
    "ridgepH7.1-ridgepH6.8=0",
    
    #pH7.1/ph6.8 in stiff
    "ridgepH7.1+ridgematrixStiff:pH7.1-ridgepH6.8+ridgematrixStiff:pH6.8=0",
    
    #Stiff/Soft
    "ridgematrixStiff+(ridgematrixStiff:pH6.8/3)+(ridgematrixStiff:pH7.1/3)=0",
    
    #pH6.8/pH6.5
    "ridgepH6.8+0.5*ridgematrixStiff:pH6.8=0",
    
    #pH7.1/pH6.5
    "ridgepH7.1+0.5*ridgematrixStiff:pH7.1=0",
    
    #pH7.1/ph6.8
    "ridgepH7.1+0.5*ridgematrixStiff:pH7.1-ridgepH6.8+0.5*ridgematrixStiff:pH6.8=0",
    
    #Interaction of Stiff/Soft with pH6.8/pH6.5
    "ridgematrixStiff:pH6.8=0", 
    
    #Interaction of Stiff/Soft with pH7.1/pH6.5
    "ridgematrixStiff:pH7.1=0", 
    
    #Interaction of Stiff/Soft with pH7.1/pH6.8
    "ridgematrixStiff:pH7.1-ridgematrixStiff:pH6.8=0"

  ),
  parameterNames = rowData(pe[["proteinRobust"]])$msqrobModels[[1]] %>%
    getCoef %>%
    names
)

##Hypothesis test
pe <- hypothesisTest(object = pe, i = "proteinRobust", contrast = L2, overwrite=TRUE)

##Significance Summary

sigsummary<-data.frame(Contrast=colnames(rowData(pe[["proteinRobust"]]))[13:28])
sigsummary$SignificantProteins=0
for(i in 13:28){
  tt<-rowData(pe[["proteinRobust"]])[,i]
  tt<-tt[!is.na(tt$adjPval),]
  sigsummary$SignificantProteins[i-12]<-nrow(tt[tt$adjPval<=0.05, ])
}
sigsummary$Description<-c("Stiff/Soft @pH6.5", 
                          "Stiff/Soft @pH6.8",
                          "Stiff/Soft @pH7.1",
                          "pH6.8/ph6.5 in Soft",
                          "pH6.8/ph6.5 in Stiff", 
                          "pH7.1/ph6.5 in Soft",
                          "pH7.1/ph6.5 in Stiff",
                          "pH7.1/ph6.8 in Soft",
                          "pH7.1/ph6.8 in Stiff",
                          "Stiff/Soft",
                          "pH6.8/pH6.5",
                          "pH7.1/pH6.5",
                          "pH7.1/ph6.8",
                          "Interaction of Stiff/Soft with pH6.8/pH6.5",
                          "Interaction of Stiff/Soft with pH7.1/pH6.5",
                          "Interaction of Stiff/Soft with pH7.1/pH6.8"
                          
)

sigsummary[, c(1,3,2)] %>%
  kbl(row.names = FALSE) %>%
  kable_classic("hover", full_width = F, font_size = 20)

write.csv(sigsummary,"sigsummary.csv")


##Protein Overlaps
genesLookup<-rowData(pe[["proteinRobust"]])$Gene.names
names(genesLookup)<-rowData(pe[["proteinRobust"]])$Proteins

res_stiff_6.5<-rowData(pe[["proteinRobust"]])[,13]
res_stiff_6.5<-res_stiff_6.5[!is.na(res_stiff_6.5$adjPval) & res_stiff_6.5$adjPval<=0.05, ]
res_stiff_6.8<-rowData(pe[["proteinRobust"]])[,14]
res_stiff_6.8<-res_stiff_6.8[!is.na(res_stiff_6.8$adjPval) & res_stiff_6.8$adjPval<=0.05, ]
res_stiff_7.1<-rowData(pe[["proteinRobust"]])[,15]
res_stiff_7.1<-res_stiff_7.1[!is.na(res_stiff_7.1$adjPval) & res_stiff_7.1$adjPval<=0.05, ]

res_stiff_all<-  rowData(pe[["proteinRobust"]])[,22]
res_stiff_all<-res_stiff_all[!is.na(res_stiff_all$adjPval) & res_stiff_all$adjPval<=0.05, ]

sigProts<-unique(c(rownames(res_stiff_6.5), rownames(res_stiff_6.8), rownames(res_stiff_7.1), rownames(res_stiff_all)))

sigProts<-unique(c(rownames(res_stiff_6.5), rownames(res_stiff_6.8), rownames(res_stiff_7.1)))
eudata<-data.frame(Prots=sigProts, Genes=genesLookup[sigProts])

eudata$pH6.5<-ifelse(eudata$Prots %in% rownames(res_stiff_6.5), TRUE, FALSE)
eudata$pH6.8<-ifelse(eudata$Prots %in% rownames(res_stiff_6.8), TRUE, FALSE)
eudata$pH7.1<-ifelse(eudata$Prots %in% rownames(res_stiff_7.1), TRUE, FALSE)
eudata$all<-ifelse(eudata$Prots %in% rownames(res_stiff_all), TRUE, FALSE)

eufit <- euler(as.matrix(eudata[, 3:5]), shape = "ellipse")
plot(eufit,
     quantities = TRUE,
     # fill = "transparent",
     lty = 1:4,
     labels = list(font = 2),col=c("red", "blue", "black", "purple"))


head(eudata[, 3:5])

write.csv(eudata,"eudata.csv")



##Volcano Plots 
comparisonLookup<-sigsummary$Description
names(comparisonLookup)<-sigsummary$Contrast


for(i in names(comparisonLookup)){
  # i<-names(comparisonLookup)[1]
  compData<-rowData(pe[["proteinRobust"]])[, i]
  compData<-compData[!is.na(compData$adjPval),]
  compData$Gene.names<-genesLookup[rownames(compData)]
  
  volcano <- ggplot(
    compData,
    aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) +
    geom_point(cex = 2) +
    scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
    theme_bw() +
    ggtitle(as.character(comparisonLookup[i]))+
    geom_text_repel(data=compData[compData$adjPval<=0.05,], aes(logFC, -log10(pval),label=Gene.names), 
                    box.padding = 1, max.overlaps = 50,nudge_x = .15, segment.curvature = -1e-20)
  # volcano
  title<-gsub(" ", "_", as.character(gsub("Interaction of ", "", gsub(" in", "", gsub("\\/", "_", gsub("\\@", "", comparisonLookup[i]))))))
  ggsave(plot=volcano, paste("C:/Users/YHU/Desktop/29.09.2022-9606/",title,"_Volcano.pdf", sep=""), height=10, width=10)
}

### new volcano

for(i in names(comparisonLookup)){
  # i<-names(comparisonLookup)[1]
  compData<-rowData(pe[["proteinRobust"]])[, i]
  compData<-compData[!is.na(compData$adjPval),]
  compData$Gene.names<-genesLookup[rownames(compData)]
  
  volcano <- ggplot(
    compData,
    aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) +
    geom_point(cex = 4) +
    scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
    theme_bw() +
    ggtitle(as.character(comparisonLookup[i]))+
    geom_text_repel(data=compData[compData$adjPval<=0.05&&logFC>=1&&logFC<=-1,], aes(logFC, -log10(pval),label=Gene.names), 
                    box.padding = 1, max.overlaps = Inf,nudge_x = .15, segment.curvature = -1e-20)+
    geom_vline(xintercept = 0,
               linetype = "dashed")+
    theme(text = element_text(size = 15)) 
  

  # volcano
  title<-gsub(" ", "_", as.character(gsub("Interaction of ", "", gsub(" in", "", gsub("\\/", "_", gsub("\\@", "", comparisonLookup[i]))))))
  ggsave(plot=volcano, paste("C:/Users/YHU/Desktop/29.09.2022-9606/",title,"_Volcano.pdf", sep=""), height=10, width=10)
}

##writeResults

wb<-createWorkbook()
for(i in names(comparisonLookup)){
  df<-rowData(pe[["proteinRobust"]])[[i]]
  df$Protein<-rownames(df)
  df$Gene.names<-genesLookup[rownames(df)]
  colnames(df)
  addDataFrame(df[, c(7:8, 1:6)], createSheet(wb, sheetName = as.character(gsub("Interaction of ", "", gsub(" in", "", gsub("\\/", " ", gsub("\\@", "", comparisonLookup[i])))))), row.names = FALSE)
  # print(as.character(gsub("Interaction of ", "", gsub(" in", "", gsub("\\/", " ", gsub("\\@", "", comparisonLookup[i]))))))
}
saveWorkbook(wb, "C:/Users/YHU/Desktop/29.09.2022-9606.xlsx")




BP<-list()
CC<-list()
MF<-list()
REACTOME<-list()

for(i in names(comparisonLookup)){
  # i = names(comparisonLookup)[1]
  print(i)
  print(comparisonLookup[i])
  compData<-rowData(pe[["proteinRobust"]])[, i]
  compData<-compData[!is.na(compData$adjPval),]
  compData$Gene.names<-genesLookup[rownames(compData)]
  
  up<-unlist(strsplit(rownames(compData[compData$adjPval<=0.05 & compData$logFC>0, ]), ";"))
  down<-unlist(strsplit(rownames(compData[compData$adjPval<=0.05 & compData$logFC<0, ]), ";"))
  bg<-unlist(strsplit(rownames(compData), ";"))
  
  BP[[paste(i, "up", sep=":")]]<-enrichGO(gene=up,universe = bg,OrgDb = "org.Hs.eg.db", ont = "BP", readable = T, keyType = "UNIPROT")
  BP[[paste(i, "down", sep=":")]]<-enrichGO(gene=down,universe = bg, OrgDb = "org.Hs.eg.db", ont = "BP", readable = T, keyType = "UNIPROT")
  
  CC[[paste(i, "up", sep=":")]]<-enrichGO(gene=up,universe = bg, OrgDb = "org.Hs.eg.db", ont = "CC", readable = T, keyType = "UNIPROT")
  CC[[paste(i, "down", sep=":")]]<-enrichGO(gene=down,universe = bg, OrgDb = "org.Hs.eg.db", ont = "CC", readable = T, keyType = "UNIPROT")
  
  MF[[paste(i, "up", sep=":")]]<-enrichGO(gene=up,universe = bg, OrgDb = "org.Hs.eg.db", ont = "MF", readable = T, keyType = "UNIPROT")
  MF[[paste(i, "down", sep=":")]]<-enrichGO(gene=down,universe = bg, OrgDb = "org.Hs.eg.db", ont = "MF", readable = T, keyType = "UNIPROT")
  
  df_up<-bitr(up, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  df_down<-bitr(down, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  df_bg<-bitr(bg, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  
  REACTOME[[paste(i, "up", sep=":")]]<-enrichPathway(gene=unique(df_up$ENTREZID), universe = unique(df_bg$ENTREZID), readable = TRUE, organism = "human")
  REACTOME[[paste(i, "down", sep=":")]]<-enrichPathway(gene=unique(df_down$ENTREZID), universe = unique(df_bg$ENTREZID), readable = TRUE, organism = "human")
}
#####
#Create report
allList<-list("BP"=BP, 
              "MF"=MF, 
              "CC"=CC, 
              "REACTOME"=REACTOME)


for(i in names(comparisonLookup)){
  wb<-createWorkbook()
  for(n in names(allList)){
    funcRes<-allList[[n]]
    if( is.null(funcRes[[paste(i, "up", sep=":")]])){
      next;
    }else{
      tt<-funcRes[[paste(i, "up", sep=":")]]@result
      addDataFrame(tt, createSheet(wb, sheetName = paste(n,"_UP", sep="")), row.names = FALSE)
      
      if(nrow(tt[tt$p.adjust<=0.05, ])>0){
        # print(paste("something here:", n))
        dd<-barplot(funcRes[[paste(i, "up", sep=":")]], showCategory = 10)
        filename<-paste("C:/Users/YHU/Desktop/New folder/",gsub("\\W", "_", comparisonLookup[i]),"_", n, "_up.pdf", sep="")
        if(nrow(tt[tt$p.adjust<=0.05, ])>10){
          ggsave(plot=dd, filename, height=8.75, width=6)
        }else{
          ggsave(plot=dd, filename, height=1+(nrow(tt[tt$p.adjust<=0.05, ])*0.75), width=6)
        }
      }
    }
    
    if( is.null(funcRes[[paste(i, "down", sep=":")]])){
      next;
    }else{
      tt<-funcRes[[paste(i, "down", sep=":")]]@result
      addDataFrame(tt, createSheet(wb, sheetName = paste(n, "_DOWN", sep="")), row.names = FALSE)
      if(nrow(tt[tt$p.adjust<=0.05, ])>0){
        # print(paste("something here:", n))
        dd<-barplot(funcRes[[paste(i, "down", sep=":")]], showCategory = 10)
        filename<-paste("C:/Users/YHU/Desktop/New folder/",gsub("\\W", "_", comparisonLookup[i]),"_", n, "_down.pdf", sep="")
        if(nrow(tt[tt$p.adjust<=0.05, ])>10){
          ggsave(plot=dd, filename, height=8.75, width=6)
        }else{
          ggsave(plot=dd, filename, height=1+(nrow(tt[tt$p.adjust<=0.05, ])*0.75), width=6)
        }
      }
    }
  }
  xlsx_filename<-paste("C:/Users/YHU/Desktop/New folder/",gsub("\\W", "_", comparisonLookup[i]), ".xlsx", sep="")
  saveWorkbook(wb, xlsx_filename)
}


###Dot plot
allList<-list("BP"=BP, 
              "MF"=MF, 
              "CC"=CC, 
              "REACTOME"=REACTOME)


for(i in names(comparisonLookup)){
  wb<-createWorkbook()
  for(n in names(allList)){
    funcRes<-allList[[n]]
    if( is.null(funcRes[[paste(i, "up", sep=":")]])){
      next;
    }else{
      tt<-funcRes[[paste(i, "up", sep=":")]]@result
      addDataFrame(tt, createSheet(wb, sheetName = paste(n,"_UP", sep="")), row.names = FALSE)
      
      if(nrow(tt[tt$p.adjust<=0.05, ])>0){
        # print(paste("something here:", n))
        dd<-dotplot(funcRes[[paste(i, "up", sep=":")]], showCategory = 10)
        filename<-paste("C:/Users/YHU/Desktop/New folder/",gsub("\\W", "_", comparisonLookup[i]),"_", n, "_up.pdf", sep="")
        if(nrow(tt[tt$p.adjust<=0.05, ])>10){
          ggsave(plot=dd, filename, height=8.75, width=6)
        }else{
          ggsave(plot=dd, filename, height=1+(nrow(tt[tt$p.adjust<=0.05, ])*0.75), width=6)
        }
      }
    }
    
    if( is.null(funcRes[[paste(i, "down", sep=":")]])){
      next;
    }else{
      tt<-funcRes[[paste(i, "down", sep=":")]]@result
      addDataFrame(tt, createSheet(wb, sheetName = paste(n, "_DOWN", sep="")), row.names = FALSE)
      if(nrow(tt[tt$p.adjust<=0.05, ])>0){
        # print(paste("something here:", n))
        dd<-dotplot(funcRes[[paste(i, "down", sep=":")]], showCategory = 10)
        filename<-paste("C:/Users/YHU/Desktop/New folder/",gsub("\\W", "_", comparisonLookup[i]),"_", n, "_down.pdf", sep="")
        if(nrow(tt[tt$p.adjust<=0.05, ])>10){
          ggsave(plot=dd, filename, height=8.75, width=6)
        }else{
          ggsave(plot=dd, filename, height=1+(nrow(tt[tt$p.adjust<=0.05, ])*0.75), width=6)
        }
      }
    }
  }
  xlsx_filename<-paste("C:/Users/YHU/Desktop/New folder/",gsub("\\W", "_", comparisonLookup[i]), ".xlsx", sep="")
  saveWorkbook(wb, xlsx_filename)
}


