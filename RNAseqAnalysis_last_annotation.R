#!/usr/bin/Rscript



#################################
#######  RNA-seq Analysis #######
#################################

# Clean environnement
ls()
rm(list=ls())
ls()



###########################################
# Installing packages bioconductor DEseq2 #
###########################################

#install.packages("pheatmap")
#install.packages("xlsx")
#install.packages("RColorBrewer")
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 BiocManager::install("IHW") 
 #BiocManager::install("IHW", version = "3.8")
source("https://bioconductor.org/biocLite.R")
biocLite("IHW")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2") 
#BiocManager::install("DESeq2", version = "3.8") 

# BiocManager::install("apeglm")
#install.packages("rggobi", type = "source")



####################
# Loading packages #
####################

library("plyr")

library("RCurl")
library("XML")
library("RSQLite")
library("DESeq2")

library("ggplot2")
library("xlsx")
library("pheatmap")
library("RColorBrewer")
library("IHW")
library("reshape2")
library("dplyr")
#library("mclust")
library("cdparcoord")
library(Rmisc)
library(GGally)
library(MASS)
library(colorspace) # get nice colors
library(lattice)
library("apeglm")
library(tidyr)
library(readr)
library(knitr)
library(kableExtra)



# working directory
setwd("~/RNAseqAnalysis/Dev_larvaire_1_new_annotation") 



################
# Loading data #
################

coldata <- read.table("design.txt", sep="\t",header=T,na.string="NA")
rownames(coldata) <- coldata$Ind
coldata <- coldata[order(rownames(coldata)),] 
coldata$Ind <- NULL
rownames(coldata)
coldata
    #             Condition
    #  J13-LU-TL1       umbo
    #  J13-LU-TL2       umbo
    #  J13-LU-TL3       umbo
    #  J22-LO-TL1a   oeillee
    #  J22-LO-TL1b   oeillee
    #  J22-LO-TL1c   oeillee
    #  J4-LD1              D
    #  J4-LD2              D
    #  J4-LD3              D
    #  J8-LV-TL1     velyger
    #  J8-LV-TL2     velyger
    #  J8-LV-TL3     velyger

cts<-read.table("join_devlarve_last_annotation_ID.txt", header=T,check.names=FALSE)
rownames(cts)<-cts$genes
cts$genes <- NULL
colnames(cts)
cts <- cts[,order(colnames(cts))] 
cts <- cts[,colnames(cts) %in% rownames(coldata) ]
cts <- cts[,order(as.numeric(colnames(cts)))]
rownames(coldata)
colnames(cts)
head(cts)

# J13-LU-TL1 J13-LU-TL2 J13-LU-TL3 J22-LO-TL1a J22-LO-TL1b J22-LO-TL1c J4-LD1 J4-LD2 J4-LD3 J8-LV-TL1
# evm.TU.scaffold10000size69197.1         19         15         35          21          26          29     23     20     39        14
# evm.TU.scaffold10000size69197.2         11          7         16           7          12          14     13     11     14         7
# evm.TU.scaffold10000size69197.3          0          2          0           4           6           4      2      0      1         5
# evm.TU.scaffold10001size69087.1         20          2         16          15          19          14     17     15     24        17
# evm.TU.scaffold10001size69087.2         40         23         51          20          40          54     43     32     44        18
# evm.TU.scaffold10002size59123.1          4          0          5           1           3           5      5      8      7         7

# Modify colnames cts
#colnames(cts) = gsub("_gsnap*.*","", colnames(cts))
#colnames(cts) = gsub("*.*---NEBNext","NEBNext", colnames(cts))



###########
## CHECK ##
###########

# check : verifie si le nombre de colonnes de design est le meme que celui de join_devlarve.txt
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
coldata$Condition<-as.factor(coldata$Condition)
str(cts)



##############
### DEseq2 ###
##############

# DESeqDataSet is a subclass of RangedSummarizedExperiment, used to store the input values, intermediate calculations and results of an analysis of differential expression
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Condition)
# Name_filtered <-rownames(dds)
colSums(cts)

# filtering
dds <- estimateSizeFactors(dds) # this function estimates the size factors using the "median ratio method" described by Equation 5 in Anders and Huber (2010)
dds$sizeFactor



#################################
### Normalisation des donnees ###
#################################

# Pour chaque colonne, v?rifie si le nombre de reads est >=10 + il faut que cela soit vrai au moins 3fois (par ligne).
idx <- rowSums(counts(dds,normalized=TRUE) >= 10 ) >= 3   
dds <- dds[idx,]                                          # Filtre sur les valeurs pr?c?dentes
dim(dds)

# Matrix log
norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 2)
#write.table(log.norm.counts,file="log_matrix.txt",quote=F)
norm.counts<-as.data.frame(norm.counts)



###########
##  PCA  ##
###########

# Data transformation : quickly estimate dispersion trend and apply a variance stabilizing transformation : remove the dependence of the variance on the mean
vsd.fast <- vst(dds, fitType='local',blind=FALSE)

### Sample PCA plot for transformed data
plotData<-plotPCA(vsd.fast,intgroup=c("Condition"), returnData=TRUE)
percentVar <- round(100 * attr(plotData, "percentVar"))
colnames(plotData)[3] = "Stages"
plotData$Stages = sub("oeillee", "Eye-spot", plotData$Stages) # Rename larval stage for the PCA plot
plotData$Stages = sub("D", "D-shaped", plotData$Stages)
plotData$Stages = sub("velyger", "Velyger", plotData$Stages)
plotData$Stages = sub("umbo", "Umbo", plotData$Stages)


PCA<-ggplot(plotData, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill=factor(Stages)),colour="black", shape=21, alpha=0.9, stroke=1, size=4) +
  scale_fill_manual(values=c('#abd9e9','#d7191c', '#fdae61', '#2c7bb6'), breaks=c("D-shaped", "Velyger", "Umbo", "Eye-spot")) +
  xlab(paste0("PC1 (",percentVar[1],"%)")) +
  ylab(paste0("PC2 (",percentVar[2],"%)")) +
  theme(axis.text.x = element_text(size=14,color="black"), 
        axis.text.y =element_text(size=14,color="black"), 
        axis.title.y = element_text(face="bold", size = 15), 
        axis.title.x = element_text(face="bold", size = 15), 
        panel.background = element_rect(fill = "white",color="black",size=1), 
        legend.position=c(0.9,0.15), 
        legend.background = element_rect(size = 0.5, linetype = "solid", color = "black", fill = "white")) +
  labs(fill = "Larval stage") + 
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed")
PCA

# Colors used :
# abd9e9 blue
# 2c7bb6 darkblue
# fdae61 orange
# d7191c red

theme_glob<-theme(axis.text.x=element_text(colour="black",size=11),
                  axis.text.y=element_text(colour="black",size=11),
                  axis.title.x=element_text(colour="black",size=12,face="bold"),
                  axis.title.y=element_text(colour="black",size=12,face="bold"),
                  panel.background = element_blank(),
                  panel.border=element_rect(fill=NA),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.margin=unit(c(1,1,1,1),"line"),
                  legend.title=element_blank(),
                  aspect.ratio=1,
                  legend.background = element_rect(size=0.3, linetype="solid", colour ="black"),
                  legend.position = c(0.86,0.13),
                  legend.key = element_rect(fill=NA),
                  legend.text = element_text(colour="black",size=12,face="bold")) 
couleurs=c("green","orange2","blue","red") 
#tiff(file = "03_results/rda.tiff", width = 20, height = 15, units = "cm", res = 300)
set.seed(1)
graph.pca<-PCA +theme_glob + scale_color_manual(values=couleurs,labels=c("D","oeillee","Umbo","velyger")) 
graph.pca
#dev.off()


##############
## PCAtools ##
##############
#install.packages("PCAtools") # Not available for R version 3.5.1.
#devtools::install_github('kevinblighe/PCAtools')
#BiocManager::install('PCAtools')

#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')
#
#BiocManager::install('PCAtools')
library(PCAtools) # Require R version 4.0.3 which work fine with DEseq2 too


vsd.assay<-assay(vst(dds, blind=TRUE))
write.table(vsd.assay, file ="vsd_assay_devlarve.txt", sep= "\t" )

mat.tot<-(as.matrix(vsd.assay))

p <- pca(mat.tot, metadata = coldata,removeVar = 0.1)
biplot(p,x='PC1',y='PC2',lab = NULL,legendPosition = 'right')
pairsplot(p)


options(ggrepel.max.overlaps = 120)
biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)



#increase max overlaps:
options(ggrepel.max.overlaps = 120)
plotloadings(p, labSize = 3)

plotloadings(p,
             rangeRetain = 0.01,
             labSize = 3.0,
             title = 'Loadings plot',
             subtitle = 'PC1, PC2, PC3, PC4, PC5',
             caption = 'Top 1% variables',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)

#eigencorplot(p,metavars = c('pH','OmCA','Length','ShellWidth','DFW30'))

horn <- parallelPCA(mat.tot)
horn$n
horn$rotated



elbow <- findElbowPoint(p$variance)
elbow


##########
# Biplot #
##########

biplot(p)

pairsplot(p)

plotloadings(p, labSize = 3)

############
# Sreeplot #
############

screeplot(p,
          components = getComponents(p, 1:20),
          vline = c(horn$n, elbow),
          hline = c(80)) +
  
  geom_label(aes(x = horn$n + 0.5, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 40,
                 label = 'Elbow method', vjust = -1, size = 8)) +
  geom_label(aes(x = 4.75, y = 69,
                 label = '80 % explain variation', vjust = -1, size = 8))

# Number of PCs that contributes to a pre-selected total of explained variation
which(cumsum(p$variance) > 80)[1] # 2PCs account for >80% explained variation.





#check unbalanced gene representation
# Most expressed sequences per sample
#'
#' Proportion of reads associated with the three most expressed sequences per sample
#'
#' @param counts \code{matrix} of counts
#' @param n number of most expressed sequences to return
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the bars (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A \code{matrix} with the percentage of reads of the three most expressed sequences and a file named majSeq.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

majSequences <- function(counts, n=3, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){
  
  seqnames <- apply(counts, 2, function(x){x <- sort(x, decreasing=TRUE); names(x)[1:n]})
  seqnames <- unique(unlist(as.character(seqnames)))
  
  sum <- apply(counts,2,sum)
  counts <- counts[seqnames,]
  sum <- matrix(sum,nrow(counts),ncol(counts),byrow=TRUE)
  p <- round(100*counts/sum,digits=3)
  
  if (outfile) png(filename="Check_Unbalanced_Gene_Representation.png",width=min(3600,1800+800*ncol(counts)/10),height=1800,res=300)
  maj <- apply(p, 2, max)
  seqname <- rownames(p)[apply(p, 2, which.max)]
  x <- barplot(maj, col=col[as.integer(group)], main="Percentage of reads from most expressed sequence",
               ylim=c(0, max(maj)*1.2), las=2, ylab="Percentage of reads")
  legend("topright", levels(group), fill=col[1:nlevels(group)], bty="n")
  for (i in 1:length(seqname)) text(x[i], maj[i]/2, seqname[i], cex=0.8, srt=90, adj=0)
  if (outfile) dev.off()
  
  return(invisible(p))
}

majSequences(counts(dds), group=coldata$Condition)



######################
## Diff. expression ##
######################

## Estimates differences
dds <- DESeq(dds, test="LRT", reduced=~1) # Likelihood ratio test
dds
res <- results(dds,alpha=0.01)
res <- res[order(res$padj),]
head(res)
resultsNames(dds)


##############################
### Export data for revigo ###
##############################

# D vs velyger
resDvel<- results(dds, contrast=c("Condition","D","velyger"),alpha=0.01) # FDR of 0.01
summary(resDvel)
resOrderedDvel <- resDvel[order(resDvel$pvalue),]     # Order by pvalue
fdrDEDvel<-subset(resOrderedDvel, padj<0.01)          # Keep only padj<0.01
finalDEDvel<-subset(fdrDEDvel,abs(log2FoldChange)>2)  # keep gene with differential expression x4
summary(finalDEDvel)
final_DV = data.frame(rownames(finalDEDvel))              # Retrieve rownames : genes names
final_DV[,2] = data.frame(finalDEDvel[,"log2FoldChange"]) # Retrieve log2FoldChange
colnames(final_DV)[1]="Gene"                              # Call first column gene
colnames(final_DV)[2]=c("log2FoldChange")                 # Call second column log2FoldChange
final_DV$Stade = "DV"                                     # Call third column by contrast : here "DV"
#write.csv(finalDEDvel,file="Larve_D_vs_Velyger.csv",quote=F)
#write.table(data.frame(rownames(finalDEDvel)),file="Larve_D_vs_Velyger.txt",quote=F,row.names=FALSE,col.names =FALSE)

# D vs Umbo
resDvel<- results(dds, contrast=c("Condition","D","umbo"),alpha=0.01) 
summary(resDvel)
resOrderedDvel <- resDvel[order(resDvel$pvalue),]
fdrDEDvel<-subset(resOrderedDvel, padj<0.01)
finalDEDvel<-subset(fdrDEDvel,abs(log2FoldChange)>2) 
summary(finalDEDvel)
final_DU = data.frame(rownames(finalDEDvel))
final_DU[,2] = data.frame(finalDEDvel[,"log2FoldChange"])
colnames(final_DU)[1]="Gene"
colnames(final_DU)[2]=c("log2FoldChange")
final_DU$Stade = "DU"
#write.csv(finalDEDvel,file="Larve_D_vs_Umbo.csv",quote=F)
#write.table(data.frame(rownames(finalDEDvel)),file="Larve_D_vs_Umbo.txt",quote=F,row.names=FALSE,col.names =FALSE)

# D vs oeillee
resDvel<- results(dds, contrast=c("Condition","D","oeillee"),alpha=0.01) 
summary(resDvel)
resOrderedDvel <- resDvel[order(resDvel$pvalue),]
fdrDEDvel<-subset(resOrderedDvel, padj<0.01)
finalDEDvel<-subset(fdrDEDvel,abs(log2FoldChange)>2) 
summary(finalDEDvel)
final_DO = data.frame(rownames(finalDEDvel))
final_DO[,2] = data.frame(finalDEDvel[,"log2FoldChange"])
colnames(final_DO)[1]="Gene"
colnames(final_DO)[2]=c("log2FoldChange")
final_DO$Stade = "DO"
#write.csv(finalDEDvel,file="Larve_D_vs_Oeillee.csv",quote=F)
#write.table(data.frame(rownames(finalDEDvel)),file="Larve_D_vs_Oeillee.txt",quote=F,row.names=FALSE,col.names =FALSE)


# Veliger vs Umbo
resDvel<- results(dds, contrast=c("Condition","velyger","umbo"),alpha=0.01) 
summary(resDvel)
resOrderedDvel <- resDvel[order(resDvel$pvalue),]
fdrDEDvel<-subset(resOrderedDvel, padj<0.01)
finalDEDvel<-subset(fdrDEDvel,abs(log2FoldChange)>2)
summary(finalDEDvel)
final_VU = data.frame(rownames(finalDEDvel))
final_VU[,2] = data.frame(finalDEDvel[,"log2FoldChange"])
colnames(final_VU)[1]="Gene"
colnames(final_VU)[2]=c("log2FoldChange")
final_VU$Stade = "VU"
#write.csv(finalDEDvel,file="Larve_velyger_vs_Umbo.csv",quote=F)
#write.table(data.frame(rownames(finalDEDvel)),file="Larve_velyger_vs_Umbo.txt",quote=F,row.names=FALSE,col.names =FALSE)

# Veliger vs oeillee
resDvel<- results(dds, contrast=c("Condition","velyger","oeillee"),alpha=0.01)
summary(resDvel)
resOrderedDvel <- resDvel[order(resDvel$pvalue),]
fdrDEDvel<-subset(resOrderedDvel, padj<0.01)
finalDEDvel<-subset(fdrDEDvel,abs(log2FoldChange)>2)
summary(finalDEDvel)
final_VO = data.frame(rownames(finalDEDvel))
final_VO[,2] = data.frame(finalDEDvel[,"log2FoldChange"])
colnames(final_VO)[1]="Gene"
colnames(final_VO)[2]=c("log2FoldChange")
final_VO$Stade = "VO"
#write.csv(finalDEDvel,file="Larve_velyger_vs_Oiellee.csv",quote=F)
#write.table(data.frame(rownames(finalDEDvel)),file="Larve_velyger_vs_Oiellee.txt",quote=F,row.names=FALSE,col.names =FALSE)

# Umbo vs oeillee
resDvel<- results(dds, contrast=c("Condition","umbo","oeillee"),alpha=0.01)
summary(resDvel)
resOrderedDvel <- resDvel[order(resDvel$pvalue),]
fdrDEDvel<-subset(resOrderedDvel, padj<0.01)
finalDEDvel<-subset(fdrDEDvel,abs(log2FoldChange)>2)
summary(finalDEDvel)
final_UO = data.frame(rownames(finalDEDvel))
final_UO[,2] = data.frame(finalDEDvel[,"log2FoldChange"])
colnames(final_UO)[1]="Gene"
colnames(final_UO)[2]=c("log2FoldChange")
final_UO$Stade = "UO"
#write.csv(finalDEDvel,file="Larve_Umbo_vs_Oiellee.csv",quote=F)
#write.table(data.frame(rownames(finalDEDvel)),file="Larve_Umbo_vs_Oiellee.txt",quote=F,row.names=FALSE,col.names =FALSE)


######################################
## Regroupe info in one data frame  ##
######################################

## Heatmap gene diff expr entre les differents stades
# Concatene tous les data.frame
AllCondtion = data.frame(rbind(final_DO, final_DU, final_DV, final_UO, final_VO, final_VU)) # concatene l'information de chaque comparaison en 1 seul tableau : AllCondition
Unique_scaffold_DEGs=data.frame(unique(AllCondtion$Gene))
colnames(Unique_scaffold_DEGs)="Genes_DEGs"
#write.table(Unique_scaffold_DEGs, file ="Unique_scaffold_DEGs.txt", sep= "\t" )

# data.frame lisible par heatmap
AllCondlog2FCsup2                           = data.frame(dcast(AllCondtion, Gene ~ Stade, value.var = "log2FoldChange"))
AllCondlog2FCsup2[is.na(AllCondlog2FCsup2)] = 0
rownames(AllCondlog2FCsup2)                 = AllCondlog2FCsup2[,1]
AllCondlog2FCsup2                           = data.frame(AllCondlog2FCsup2[,c("DO", "DU", "DV", "UO", "VO", "VU")]) # Keep this data frame for EDG (ExprDiffGenes)



# Retrieve matrix log 
norm.counts     <- counts(dds, normalized=TRUE)   # Retrieve normalized countcounts
log.norm.counts <- log2(norm.counts + 2)          # Log2 on normalized count

# norm.counts<-as.data.frame(norm.counts)
logncsubset = log.norm.counts[which(rownames(log.norm.counts) %in% rownames(AllCondlog2FCsup2)),] # Retrieve DEG (ExprDiffGenes) with AllCondlog2FCsup2
colnames(logncsubset) = c("Umbo_1", "Umbo_2", "Umbo_3", "Eye-spot_1", "Eye-spot_2", "Eye-spot_3", "D-shaped_1", "D-shaped_2", "D-shaped_3", "Velyger_1", "Velyger_2", "Velyger_3")

logncsubset_minRowMeans = logncsubset - rowMeans(logncsubset) # En soustrayant la moyenne au log2
subj                    = colnames(logncsubset_minRowMeans)   # Retrieve colnames of logncsubset_minRowMeans for aka2
aka2 = data.frame(Condition = factor(c("Umbo", "Umbo", "Umbo", "Eye-spot", "Eye-spot", "Eye-spot", "D-shaped", "D-shaped", "D-shaped", "Velyger", "Velyger", "Velyger")))
rownames(aka2)          = subj                                # retrieve colnames for aka2
my_colour= list(Condition = c('D-shaped'="#abd9e9", Velyger="#2c7bb6", Umbo="#fdae61", 'Eye-spot'="#d7191c"))

## Colors
# "#2c7bb6" : bleu fonc?
# "#abd9e9" : bleu clair
# "#fdae61" : orange
# "#d7191c" : rouge


##################################
# Define best number of clusters #
##################################

# Bayesian approach
#d_clust <- Mclust(as.matrix(logncsubset_minRowMeans), G=1:15, modelNames = mclust.options("emModelNames"))
#d_clust$BIC
#plot(d_clust)



################
###  HEATMAP  ##
################

set.seed(3) # to fix the random number used in heatmap and always keep the same cluster after a restart cession.
# K-mean fix to 8
out<-pheatmap(logncsubset_minRowMeans,kmeans_k=8,
               annotation_col = aka2,
               annotation_colors = my_colour,
               annotation_names_col = FALSE,
               cellheight = 20,
               cellwidth = 20,
               angle_col = 45
              )
#ggsave("Rplot_Heatmap_larval_stages.pdf", plot = out, device = "pdf", width = 8, height = 12)




#################
## VIOLIN PLOT ##
#################
mycluster               = data.frame(out$kmeans$cluster)  # Retrieve cluster from pheatmap
colnames(mycluster)[1]  = "cluster"                       # Rename first column cluster
mycluster$Gene          = rownames(mycluster)             # Add rownames (Gene info) to the data.frame

tmp=data.frame(mycluster[order(mycluster$cluster),])
tmp2=data.frame(tmp[,1])
rownames(tmp2)=rownames(tmp)
colnames(tmp2)="cluster"

pheatmap(logncsubset_minRowMeans,
         annotation_row = tmp2[,],
         cutree_rows = 2
)





logdata_EachStadLarv      = data.frame(logncsubset_minRowMeans) # change name
logdata_EachStadLarv$Gene = rownames(logdata_EachStadLarv)      # Add rownames as columns (genes) 


p<-list()

pdf("Rplot_violinPerCluster.pdf")
  for (i in 1:max(mycluster$cluster)) {       # For each cluster : add data from log2_ncount
    print("################")
    Cluster_name = paste("Cluster", i, sep = "_") # Generate Cluster_name : keep cluster number for plotting
    print(Cluster_name)
    Cluster <- join(data.frame(mycluster[which(mycluster==i),]), logdata_EachStadLarv) # Add logdata_EachStadLarv (log2-mean) to each cluster => Cluster
    print("Dim :")
    # join logdata_EachStadLarv (log2 - mean) for each cluster and put it in variable Cluster_[i]
    print(dim(assign(paste("Cluster", i, sep = "_"), join(data.frame(mycluster[which(mycluster==i),]), logdata_EachStadLarv)))) 
    
    # Format data 
    nm_Cluster = Cluster
    nm_Cluster[,3:length(colnames(nm_Cluster))] = round(nm_Cluster[,3:length(colnames(nm_Cluster))], digits = 4) # round data with 4 digits
    tmp = nm_Cluster[,-1]                 # Delete first column
    tmp = tmp[,-1]                        # Delete first column again ... 
    tmp = tmp[,order(colnames(tmp))]      # order tmp's colnames 
    tmp$'D-shaped' = rowMeans(tmp[,1:3])  # Mean of triplicats for stade D
    tmp$'Eye-spot' = rowMeans(tmp[,4:6])  # Mean of triplicats for stade Oeillee
    tmp$Umbo    = rowMeans(tmp[,7:9])     # Mean of triplicats for stade umbo
    tmp$Velyger = rowMeans(tmp[,10:12])   # Mean of triplicats for stade velyger
    tmp2        = tmp[,c("D-shaped", "Velyger", "Umbo", "Eye-spot")]    # Order new columns "means" from triplicats into a new data.Frame tmp2
    nm_Cluster_mean               = cbind(tmp2, nm_Cluster$Gene)        # Add Gene information to tmp2
    colnames(nm_Cluster_mean)[5]  = "Gene"                              # Rename last column into "Gene"
    nm_Cluster_mean$Genes         = as.character(nm_Cluster_mean$Gene)  # Change Gene column (factor) into Genes column (character)
    nm_Cluster_mean               = nm_Cluster_mean[,-5]                # Delete Gene column (factor)

    write.table(x = nm_Cluster_mean, paste("Cluster", i, "_Genes_logFc.txt", sep = ""), quote = F, row.names = FALSE, col.names = FALSE)
    m_data = melt(nm_Cluster_mean)
    print(ggplot(m_data, aes(x = variable, y = value, color=variable)) + 
            geom_violin(trim=FALSE) + ylab("Log2 - mean") + xlab("") + 
            ggtitle(paste(Cluster_name, paste(nrow(nm_Cluster_mean), "genes", sep = " "), "Boxplot for each larval stages on genes differentially expressed", sep = " / ")) + 
            theme(plot.title = element_text(size=10, face = "bold")) + 
            theme(axis.text.x = element_text(size=15,color="black"), axis.text.y =element_text(size=15,color="black"), axis.title.y = element_text(face="plain",size=15), panel.background = element_rect(fill = "white",color="black",size=1), panel.grid = element_blank(), legend.position="none") + 
            stat_summary(fun.y=mean, geom = "point", color="black", size=3))
    
    
    p[[i]]<-ggplot(m_data, aes(x = variable, y = value, fill = variable)) + 
            geom_violin(trim=FALSE) + 
            geom_boxplot(width=0.1, fill="white") +
            ylab("Log2 - mean") + 
            xlab("") + 
            ggtitle(paste("Cluster ", i, " (", nrow(nm_Cluster_mean), ")", sep = "")) + 
            theme(plot.title = element_text(size=14, face = "bold", vjust = -9, hjust= 0.025)) + 
            theme(axis.text.x = element_text(size=15,color="black"), axis.text.y =element_text(size=15,color="black"), axis.title.y = element_text(face="plain",size=15), panel.background = element_rect(fill = "white",color="black",size=1), panel.grid = element_blank(), legend.position="none") +
            scale_fill_manual(values=c("#abd9e9", "#2c7bb6", "#fdae61", "#d7191c")) +
            stat_summary(fun.y=mean, geom = "point", color="red", size=1)
              
  }
dev.off()

layout = matrix(c(1,2,3,4,5,6,7,8), nrow = 2, byrow = TRUE)
multiplot(p[[1]],p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], p[[7]], p[[8]], layout = layout)


##################################################
## MAKE PRETTY TABLE FOR DIPLAYING DEGs RESULTS ##
##################################################
install.packages("kableExtra")
library(kableExtra)



################
### TABLEAUX ###
################

# Figure nombre de reads / reads count
tmp = read.delim(file = "clipboard", stringsAsFactors = F, check.names = F, header = T)
colnames(tmp)=c("Sample", "Row reads fastq", "Trim reads fastq", "Mapped reads")

knitr::kable(tmp, align = 'r') %>%
      kable_classic(full_width = F, html_font = "Cambria")



###


































######################################################################################################################

                                                ### ZONE DE TEST  ###

######################################################################################################################



#######################

###  KEGG PATHWAY  ###

######################

source("https://bioconductor.org/biocLite.R")
biocLite(c("topGO", "KEGGREST", "org.At.tair.db", "Rgraphviz"))
library(topGO)
library(KEGGREST)


infile <- "I5_v_C_time6.txt"
tmp <- read.delim(infile)

# OR
#tmp <- read.delim("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2018-June-RNA-Seq-Workshop/master/friday/I5_v_C_time6.txt")

geneList <- tmp$P.Value
names(geneList) <- tmp$Gene
# Create topGOData object
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSelectionFun = function(x)x,
              annot = annFUN.org, mapping = "org.At.tair.db")

# Kolmogorov-Smirnov testing
resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
head(tab, 15)


rna.meth.terms <- genesInTerm(GOdata)[["GO:0001510"]] # get genes associated with term
p.values.in <- geneList[names(geneList) %in% rna.meth.terms]
p.values.out <- geneList[!(names(geneList) %in% rna.meth.terms)]
plot.ecdf(p.values.in, verticals = T, do.points = F, col = "red", lwd = 2, xlim = c(0,1),
          main = "Empirical Distribution of DE P-Values by Annotation with 'RNA Methylation'",
          cex.main = 0.9, xlab = "p", ylab = "Probabilty(P-Value < p)")
ecdf.out <- ecdf(p.values.out)
xx <- unique(sort(c(seq(0, 1, length = 201), knots(ecdf.out))))
lines(xx, ecdf.out(xx), col = "black", lwd = 2)
legend("bottomright", legend = c("Genes Annotated with 'RNA Methylation'", "Genes Not Annotated with 'RNA Methylation'"), lwd = 2, col = 2:1, cex = 0.9)

par(cex = 0.3)
showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 3, useInfo = "def")
par(cex = 1)




##################???
## 4 KEGG Pathway #
###################


# Pull all pathways for AT  
pathways.list <- keggList("pathway")
head(pathways.list)

# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list)) 
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           }
)

head(genes.by.pathway)

infile <- "I5_v_C_time6.txt"
DE.table <- read.delim(infile, stringsAsFactors = F)

# OR

#DE.table <- read.delim("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2018-June-RNA-Seq-Workshop/master/friday/I5_v_C_time6.txt")

geneList <- DE.table$P.Value
names(geneList) <- DE.table$Gene
head(geneList)

# Wilcoxon test for each pathway
pVals.by.pathway <- t(sapply(names(genes.by.pathway),
                             function(pathway) {
                               pathway.genes <- genes.by.pathway[[pathway]]
                               list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
                               list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
                               scores.in.pathway <- geneList[list.genes.in.pathway]
                               scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
                               if (length(scores.in.pathway) > 0){
                                 p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
                               } else{
                                 p.value <- NA
                               }
                               return(c(p.value = p.value, Annotated = length(list.genes.in.pathway)))
                             }
))

# Assemble output table
outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
outdat$pathway.name <- pathways.list[outdat$pathway.code]
outdat$p.value <- pVals.by.pathway[,"p.value"]
outdat$Annotated <- pVals.by.pathway[,"Annotated"]
outdat <- outdat[order(outdat$p.value),]
head(outdat)

sessionInfo()







tmp = nm_Cluster_mean
tmp = (nm_Cluster_mean[,-5])
tmp2 = melt(tmp)

fit <- aov(value ~ variable, data=tmp2) # y est la variable numérique et A indique les groupes
summary(fit)


pdf("Rplot_BowplotPerCluster.pdf")
for (i in 1:max(mycluster$cluster)) {       # For each cluster : add data from log2_ncount
  print("################")
  Cluster_name = paste("Cluster", i, sep = "_") # Generate Cluster_name : keep cluster number for plotting
  print(Cluster_name)
  Cluster <- join(data.frame(mycluster[which(mycluster==i),]), logdata_EachStadLarv) # Add logdata_EachStadLarv (log2-mean) to each cluster => Cluster
  print("Dim :")
  # join logdata_EachStadLarv (log2 - mean) for each cluster and put it in variable Cluster_[i]
  print(dim(assign(paste("Cluster", i, sep = "_"), join(data.frame(mycluster[which(mycluster==i),]), logdata_EachStadLarv)))) 
  
  # Format data 
  nm_Cluster = Cluster
  nm_Cluster[,3:length(colnames(nm_Cluster))] = round(nm_Cluster[,3:length(colnames(nm_Cluster))], digits = 4) # round data with 4 digits
  tmp = nm_Cluster[,-1]                 # Delete first column
  tmp = tmp[,-1]                        # Delete first column again ... 
  tmp = tmp[,order(colnames(tmp))]      # order tmp's colnames 
  tmp$'D-shaped' = rowMeans(tmp[,1:3])  # Mean of triplicats for stade D
  tmp$'Eye-spot' = rowMeans(tmp[,4:6])  # Mean of triplicats for stade Oeillee
  tmp$Umbo    = rowMeans(tmp[,7:9])     # Mean of triplicats for stade umbo
  tmp$Velyger = rowMeans(tmp[,10:12])   # Mean of triplicats for stade velyger
  tmp2        = tmp[,c("D-shaped", "Velyger", "Umbo", "Eye-spot")]    # Order new columns "means" from triplicats into a new data.Frame tmp2
  nm_Cluster_mean               = cbind(tmp2, nm_Cluster$Gene)        # Add Gene information to tmp2
  colnames(nm_Cluster_mean)[5]  = "Gene"                              # Rename last column into "Gene"
  nm_Cluster_mean$Genes         = as.character(nm_Cluster_mean$Gene)  # Change Gene column (factor) into Genes column (character)
  nm_Cluster_mean               = nm_Cluster_mean[,-5]                # Delete Gene column (factor)
  
  write.table(x = nm_Cluster_mean, paste("Cluster", i, "_Genes_logFc.txt", sep = ""), quote = F, row.names = FALSE, col.names = FALSE)
  m_data = melt(nm_Cluster_mean)
  print(ggplot(m_data, aes(x = variable, y = value, color=variable)) + 
          geom_boxplot(trim=FALSE) + ylab("Log2 - mean") + xlab("") + 
          ggtitle(paste(Cluster_name, paste(nrow(nm_Cluster_mean), "genes", sep = " "), "Boxplot for each larval stages on genes differentially expressed", sep = " / ")) + 
          theme(plot.title = element_text(size=10, face = "bold")) + 
          theme(axis.text.x = element_text(size=15,color="black"), axis.text.y =element_text(size=15,color="black"), axis.title.y = element_text(face="plain",size=15), panel.background = element_rect(fill = "white",color="black",size=1), panel.grid = element_blank(), legend.position="none") + 
          stat_summary(fun.y=mean, geom = "point", color="black", size=3))
  
  
  p[[i]]<-ggplot(m_data, aes(x = variable, y = value, fill = variable)) + 
    geom_boxplot(trim=FALSE) +
    ylab("Log2 - mean") + 
    xlab("") + 
    ggtitle(paste("Cluster ", i, " (", nrow(nm_Cluster_mean), ")", sep = "")) + 
    theme(plot.title = element_text(size=14, face = "bold", vjust = -9, hjust= 0.025)) + 
    theme(axis.text.x = element_text(size=15,color="black"), axis.text.y =element_text(size=15,color="black"), axis.title.y = element_text(face="plain",size=15), panel.background = element_rect(fill = "white",color="black",size=1), panel.grid = element_blank(), legend.position="none") +
    scale_fill_manual(values=c("#abd9e9", "#2c7bb6", "#fdae61", "#d7191c")) +
    stat_summary(fun.y=mean, geom = "point", color="red", size=1)
  
}
dev.off()

layout = matrix(c(1,2,3,4,5,6,7,8), nrow = 2, byrow = TRUE)
multiplot(p[[1]],p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], p[[7]], p[[8]], layout = layout)


#############################



########################################################
# cluster : sortir les genes pour chacuns des clusters #
########################################################

head(mycluster)

for (i in 1:max(mycluster$cluster)) {
  print(i)
  ClusterFile = data.frame(rownames(mycluster[which(mycluster$cluster==i),]))
  colnames(ClusterFile)[1] = "Genes"
  write.table(x = ClusterFile, paste("GenesCluster", i, ".txt", sep = ""), quote = F, row.names = FALSE, col.names = FALSE)
}




#################################
# Work on gene list from Jeremy #
#################################

# retrieve info from annotation.txt with a list of genes.
res_annotation  = read.table("annotation.txt", sep="\t", header = T, na.strings = "NA")
My318DEG_biomin = read.table("Fichier_GO_plus_genes_F_Biomin.txt", sep="\t", header = F, na.strings = "NA")
My1944DEGs      = data.frame(rownames(AllCondlog2FCsup2))
colnames(My318DEG_biomin) = "Transcripts"
colnames(My1944DEGs)      = "Transcripts"
My1944DEGs_complete = res_annotation[which(res_annotation$Transcripts %in% My1944DEGs$Transcripts),]
My1944DEGs_complete$Cgigas = gsub("\\|", "",gsub(".*ref|", "", My1944DEGs_complete$Cgigas))
My1944DEGs_complete$Best_hit= gsub("\\|", "",gsub(".*ref|", "", My1944DEGs_complete$Best_hit))
df=data.frame(t(data.frame(strsplit(as.character(My1944DEGs_complete$Sprot), '\\|'))))
tmp=cbind(My1944DEGs_complete, df)
colnames(tmp)[6:8]=c("sp", "Uniprot_Accession", "Uniprot_ID")
tmp = tmp[,c("Transcripts", "Uniprot_Accession", "Uniprot_ID", "Pfucata", "Cgigas",  "Best_hit")]
tmp$Biomin_gene = 0 # Add a new
tmp[which(tmp$Transcripts %in% My318DEG_biomin$Transcripts),"Biomin_gene"]=1

Genes_biomin_Jeremy = read.delim(file = "Gene_Biomin_Jeremy.txt", stringsAsFactors = F, check.names = F, header = F)
colnames(Genes_biomin_Jeremy)= "Transcripts"
Genes_biomin_Jeremy$Biomin_gene_jeremy = 1
Genes_biomin_1944DEGs_withJeremy = join(tmp, Genes_biomin_Jeremy)
Genes_biomin_1944DEGs_withJeremy[which(is.na(Genes_biomin_1944DEGs_withJeremy$Biomin_gene_jeremy)), "Biomin_gene_jeremy"] = 0
Genes_biomin_1944DEGs_withJeremy[which(Genes_biomin_1944DEGs_withJeremy$Biomin_gene_jeremy==1),]




######################################################
##   Work with Gene associated with selected genes  ##
######################################################

# retrieve info from annotation.txt with a list of genes.
My1944DEGs_name = data.frame(rownames(AllCondlog2FCsup2)) 
res_annotation  = read.table("annotation.txt", sep="\t", header = T, na.strings = "NA")
My318DEG_Biomin         = read.table("Fichier_GO_plus_genes_F_Biomin.txt", sep="\t", header = F, na.strings = "NA")       # Retrieve genes names associated with biomin obtain with rervigo
My318DEG_ImmuneSystem   = read.table("Fichier_GO_plus_genes_F_ImmuneSytem.txt", sep="\t", header = F, na.strings = "NA")  # Retrieve genes names associated with Immune system obtain with rervigo
My318DEG_Perception     = read.table("Fichier_GO_plus_genes_F_perception.txt", sep="\t", header = F, na.strings = "NA")   # Retrieve genes names associated with Perception obtain with rervigo
colnames(My1944DEGs_name)       = "Transcripts"
colnames(My318DEG_Biomin)       = "Transcripts"
colnames(My318DEG_ImmuneSystem) = "Transcripts"
colnames(My318DEG_Perception)   = "Transcripts"


# Format data to keep ID / accession from uniprot and Pfucata, Cgigas and Best_hit
My1944DEGs_annotation = res_annotation[which(res_annotation$Transcripts %in% My1944DEGs_name$Transcripts),] # Ajout de l'information contenu dans le fichier annotation (annotation genome)
#My1944DEGs_annotation$Cgigas = gsub("\\|", "",gsub(".*ref|", "", My1944DEGs_annotation$Cgigas))            # Elimine de l'info pour la colonne Cgigas
#My1944DEGs_annotation$Best_hit= gsub("\\|", "",gsub(".*ref|", "", My1944DEGs_annotation$Best_hit))         # Elimine de l'info pour la colonne Best_hit
df=data.frame(t(data.frame(strsplit(as.character(My1944DEGs_annotation$Sprot), '\\|'))))                    # split column Sprot to retrieve uniprot_Accession and uniprot_ID
tmp=cbind(My1944DEGs_annotation, df)                                                                        # Merge cette info avec notre data_frame contenant l'annotation du g?nome
colnames(tmp)[6:8]=c("sp", "Uniprot_Accession", "Uniprot_ID")                                               # Renome les colonnes split?es 
tmp = tmp[,c("Transcripts", "Uniprot_Accession", "Uniprot_ID", "Sprot", "Pfucata", "Cgigas",  "Best_hit")]  # R?ordonne et r?cup?re les colonnes d'interet
tmp$Biomin_gene       = 0                                                                                   # Add a new column for Biomin genes avec pour valeur 0
tmp$ImmuneSystem_gene = 0                                                                                   # Add a new column for Immune system genes avec pour valeur 0
tmp$Perception_gene   = 0                                                                                   # Add a new column for Perception genes avec pour valeur 0
tmp[which(tmp$Transcripts %in% My318DEG_Biomin$Transcripts),"Biomin_gene"]=1                                # Ajout de la valeur 1 si associ? ? la biomin
tmp[which(tmp$Transcripts %in% My318DEG_ImmuneSystem$Transcripts),"ImmuneSystem_gene"]=1                    # Ajout de la valeur 1 si associ? au systeme immunitaire
tmp[which(tmp$Transcripts %in% My318DEG_Perception$Transcripts),"Perception_gene"]=1                        # Ajout de la valeur 1 si associ? ? la perception
rownames(tmp)=NULL                                                                                          # R?initiallise le nom des lignes
nrow(tmp[(which(tmp$Biomin_gene==1)),])                                                                     # Compte le nombre de genes associ? ? la biomin
nrow(tmp[(which(tmp$ImmuneSystem_gene==1)),])                                                               # Compte le nombre de genes associ? au systeme immunitaire
nrow(tmp[(which(tmp$Perception_gene==1)),])                                                                 # Compte le nombre de genes associ? ? la perception
My1944DEGs_annotation = tmp                                                                                 # Sauvegarde toute cette info dans My1994DEGs_annotation



###################################################
# Select genes associated with biomineralization: #
###################################################

# We add gene associated with biomin in P.margaratifera from jeremy list
Protein_Biomin_selected = read.delim(file = "35proteins_Biomin.txt", stringsAsFactors = F, check.names = F, header = F) # selectionnees ? la main.
Protein_Biomin_selected$BiominSelected = 1                                      # Rajoute une colonne avec les genes selectionn?s qui correpondent ? 1
colnames(Protein_Biomin_selected)[1]="Uniprot_Accession"                        # Renomme la 1ere colonne
My1944DEGs_annotation = join(My1944DEGs_annotation, Protein_Biomin_selected)    # Merge/join 2 data.frame
nrow(My1944DEGs_annotation[which(My1944DEGs_annotation$BiominSelected==1), ])   # check number of genes selected
# Duplicated (one protein match with 2 genes TRINITY in our annotation file...)
My1944DEGs_annotation[which(My1944DEGs_annotation$BiominSelected==1), ][which(duplicated(My1944DEGs_annotation[which(My1944DEGs_annotation$BiominSelected==1), "Uniprot_Accession"])),]
nrow(My1944DEGs_annotation[which(My1944DEGs_annotation$BiominSelected==1), ][which(duplicated(My1944DEGs_annotation[which(My1944DEGs_annotation$BiominSelected==1), "Uniprot_Accession"])),])
#Nombre de genes en elevant les duplicated :
nrow(My1944DEGs_annotation[which(My1944DEGs_annotation$BiominSelected==1), ]) - nrow(My1944DEGs_annotation[which(My1944DEGs_annotation$BiominSelected==1), ][which(duplicated(My1944DEGs_annotation[which(My1944DEGs_annotation$BiominSelected==1), "Uniprot_Accession"])),])
# Fichier contenant les genes uniques:
My1944DEGs_annotation[which(My1944DEGs_annotation$BiominSelected==1), ][which(duplicated(My1944DEGs_annotation[which(My1944DEGs_annotation$BiominSelected==1), "Uniprot_Accession"])),]$BiominSelected = 0



###############################################
# Select genes associated with Immune system: #
###############################################

Protein_Immune_system = read.delim(file = "57proteins_ImmuneSystem.txt", stringsAsFactors = F,check.names = F, header = F)
Protein_Immune_system$ImmuneSystemSelected = 1
colnames(Protein_Immune_system)[1]="Uniprot_Accession"
My1944DEGs_annotation = join(My1944DEGs_annotation, Protein_Immune_system)
nrow(My1944DEGs_annotation[which(My1944DEGs_annotation$ImmuneSystemSelected==1), ])
# Duplicated (one protein match with 2 genes TRINITY in our annotation file...)
My1944DEGs_annotation[which(My1944DEGs_annotation$ImmuneSystemSelected==1), ][which(duplicated(My1944DEGs_annotation[which(My1944DEGs_annotation$ImmuneSystemSelected==1), "Uniprot_Accession"])),]
nrow(My1944DEGs_annotation[which(My1944DEGs_annotation$ImmuneSystemSelected==1), ][which(duplicated(My1944DEGs_annotation[which(My1944DEGs_annotation$ImmuneSystemSelected==1), "Uniprot_Accession"])),])
# Nombre de genes en enlevant les duplicated :
nrow(My1944DEGs_annotation[which(My1944DEGs_annotation$ImmuneSystemSelected==1), ]) - nrow(My1944DEGs_annotation[which(My1944DEGs_annotation$ImmuneSystemSelected==1), ][which(duplicated(My1944DEGs_annotation[which(My1944DEGs_annotation$ImmuneSystemSelected==1), "Uniprot_Accession"])),])
# Fichier contenant les genes uniques:
My1944DEGs_annotation[which(My1944DEGs_annotation$ImmuneSystemSelected==1), ][which(duplicated(My1944DEGs_annotation[which(My1944DEGs_annotation$ImmuneSystemSelected==1), "Uniprot_Accession"])),]$ImmuneSystemSelected = 0



############################################
# Select genes associated with Perception: #
############################################

Protein_Perception = read.delim(file = "52Protein_Perception.txt", stringsAsFactors = F,check.names = F, header = F)
Protein_Perception$PerceptionSelected = 1
colnames(Protein_Perception)[1]="Uniprot_Accession"
My1944DEGs_annotation = join(My1944DEGs_annotation, Protein_Perception)
nrow(My1944DEGs_annotation[which(My1944DEGs_annotation$PerceptionSelected==1), ])
# Duplicated (one protein match with 2 genes TRINITY in our annotation file...)
My1944DEGs_annotation[which(My1944DEGs_annotation$PerceptionSelected==1), ][which(duplicated(My1944DEGs_annotation[which(My1944DEGs_annotation$PerceptionSelected==1), "Uniprot_Accession"])),]
nrow(My1944DEGs_annotation[which(My1944DEGs_annotation$PerceptionSelected==1), ][which(duplicated(My1944DEGs_annotation[which(My1944DEGs_annotation$PerceptionSelected==1), "Uniprot_Accession"])),])
# Nombre de genes en enlevant les duplicated :
nrow(My1944DEGs_annotation[which(My1944DEGs_annotation$PerceptionSelected==1), ]) - nrow(My1944DEGs_annotation[which(My1944DEGs_annotation$PerceptionSelected==1), ][which(duplicated(My1944DEGs_annotation[which(My1944DEGs_annotation$PerceptionSelected==1), "Uniprot_Accession"])),])
# Fichier contenant les genes uniques:
My1944DEGs_annotation[which(My1944DEGs_annotation$PerceptionSelected==1), ][which(duplicated(My1944DEGs_annotation[which(My1944DEGs_annotation$PerceptionSelected==1), "Uniprot_Accession"])),]$PerceptionSelected = 0

# Recapitulatif:
nrow(My1944DEGs_annotation[which(My1944DEGs_annotation$BiominSelected==1),])
nrow(My1944DEGs_annotation[which(My1944DEGs_annotation$ImmuneSystemSelected==1),])
nrow(My1944DEGs_annotation[which(My1944DEGs_annotation$PerceptionSelected==1),])
head(My1944DEGs_annotation)



################################################################
##  Retrieve gene of interest and place them in the cluster 9 ##
################################################################

# Biomin
my8cluster = mycluster

mycluster = my8cluster
mycluster$cluster=0
tmp = as.character(My1944DEGs_annotation[which(My1944DEGs_annotation$BiominSelected==1),"Transcripts"])
mycluster[which(mycluster$Gene %in% tmp),"cluster"]=2
nrow(mycluster[which(mycluster$cluster==2),])
mycluster_biomin = mycluster

# Immune system
mycluster = my8cluster
mycluster$cluster=0
tmp = as.character(My1944DEGs_annotation[which(My1944DEGs_annotation$ImmuneSystemSelected==1),"Transcripts"])
mycluster[which(mycluster$Gene %in% tmp),"cluster"]=2
nrow(mycluster[which(mycluster$cluster==2),])
mycluster_ImmuneSystem = mycluster

# Perception
mycluster = my8cluster
mycluster$cluster=0
tmp = as.character(My1944DEGs_annotation[which(My1944DEGs_annotation$PerceptionSelected==1),"Transcripts"])
mycluster[which(mycluster$Gene %in% tmp),"cluster"]=2
nrow(mycluster[which(mycluster$cluster==2),])
mycluster_Perception = mycluster

my3cluster = list(mycluster_biomin, mycluster_ImmuneSystem, mycluster_Perception)
names(my3cluster) = c("mycluster_biomin", "mycluster_ImmuneSystem", "mycluster_Perception")



#########################################################################################################
## Expression of interest genes during larval development focus on Biomin, immune system and perception #
#########################################################################################################

for (i in 1:3) {
  name_cluster = print(names(my3cluster)[i])
  name_cluster = gsub("mycluster_", "", name_cluster)
  my3cluster[[i]]
  print(which(my3cluster[[i]]$cluster==2))
  print(my3cluster[[i]][which(my3cluster[[i]]$cluster==2),])
  
  Cluster <- join(data.frame(my3cluster[[i]][which(my3cluster[[i]]$cluster==2),]), logdata_EachStadLarv) # Add logdata_EachStadLarv (log2-mean) to each cluster => Cluster
    
  # Format data 
  nm_Cluster = Cluster
  nm_Cluster[,3:length(colnames(nm_Cluster))] = round(nm_Cluster[,3:length(colnames(nm_Cluster))], digits = 4) # round data with 4 digits
  tmp = nm_Cluster[,-1]                     # Delete first column
  tmp = tmp[,-1]                            # Delete first column again ... 
  tmp = tmp[,order(colnames(tmp))]          # order tmp's colnames 
  tmp$'D-shaped'  = rowMeans(tmp[,1:3])     # Mean of triplicats for stade D
  tmp$'Eye-spot'  = rowMeans(tmp[,4:6])     # Mean of triplicats for stade Oeillee
  tmp$Umbo        = rowMeans(tmp[,7:9])     # Mean of triplicats for stade umbo
  tmp$Velyger     = rowMeans(tmp[,10:12])   # Mean of triplicats for stade velyger
  tmp2            = tmp[,c("D-shaped", "Velyger", "Umbo", "Eye-spot")]    # Order new columns "means" from triplicats into a new data.Frame tmp2
  nm_Cluster_mean               = cbind(tmp2, nm_Cluster$Gene)            # Add Gene information to tmp2
  colnames(nm_Cluster_mean)[5]  = "Gene"                                  # Rename last column into "Gene"
  nm_Cluster_mean$Genes         = as.character(nm_Cluster_mean$Gene)      # Change Gene column (factor) into Genes column (character)
  nm_Cluster_mean               = nm_Cluster_mean[,-5]                    # Delete Gene column (factor) 
  
  # Rearrange tmp pour faire une heatmap et determiner le nombre de cluster/group parmis les 29 gene de biomin
  tmp=nm_Cluster_mean
  rownames(tmp)=tmp$Genes
  
  set.seed(1) # to fix the random number used in heatmap and always keep the same cluster after a restart cession.
  out_biomin<-pheatmap(tmp[1:4],kmeans_k = 3,
                       cellheight = 20,
                       cellwidth = 20,
                       angle_col = 45
  )
  
  mycluster_biomin        = data.frame(out_biomin$kmeans$cluster)  # Retrieve cluster from pheatmap
  colnames(mycluster_biomin)[1]  = "cluster"                       # Rename first column cluster
  mycluster_biomin$Gene          = rownames(mycluster_biomin)      # Add rownames (Gene info) to the data.frame
  tmp = cbind(tmp,mycluster_biomin)
  
  levels <-  tmp %>%
    dplyr::select(Gene, "D-shaped", Velyger, Umbo, "Eye-spot", cluster) %>%
    gather(key = "stade", value = "expr.level", -Gene, -cluster) %>%
    mutate(num.stade = case_when(stade == "D-shaped" ~ 1,
                                 stade == "Velyger" ~ 2,
                                 stade == "Umbo" ~ 3,
                                 stade == "Eye-spot" ~ 4))
  p[[i]]<-ggplot(levels, aes(x = num.stade, y = expr.level,fill = Gene, color = as.factor(cluster))) + geom_line() + 
    labs(color="Cluster") + scale_x_discrete(name=paste("Cluster", name_cluster, sep = "_"), limits=c("D-shaped", "Veliger", "Umbo", "Eye-spot")) +
    scale_y_continuous(name = "Expression level", limits = c(-5,5))
  
  m_data = melt(nm_Cluster_mean)
  ggplot(m_data, aes(x = variable, y = value, fill = variable)) + 
    geom_violin(trim=FALSE) + 
    geom_boxplot(width=0.1, fill="white") +
    ylab("Log2 - mean") + 
    xlab("") + 
    ggtitle(paste("Cluster ", i, " (", nrow(nm_Cluster_mean), ")", sep = "")) + 
    theme(plot.title = element_text(size=14, face = "bold", vjust = -9, hjust= 0.025)) + 
    theme(axis.text.x = element_text(size=15,color="black"), axis.text.y =element_text(size=15,color="black"), axis.title.y = element_text(face="plain",size=15), panel.background = element_rect(fill = "white",color="black",size=1), panel.grid = element_blank(), legend.position="none") +
    scale_fill_manual(values=c("#abd9e9", "#2c7bb6", "#fdae61", "#d7191c")) +
    stat_summary(fun.y=mean, geom = "point", color="red", size=1)

}
dev.off()
layout = matrix(c(1,2,3), nrow = 1, byrow = FALSE)
multiplot(p[[1]],p[[2]], p[[3]], layout = layout)
# Save as pdf 6x20

## Same goal but for each cluster separated:
# Plot only the genes of cluster 1 :
print(ggplot(levels[which(levels$cluster==1),], aes(x = num.stade, y = expr.level,fill = Gene, color = as.factor(cluster))) + 
        geom_line(color = "#cc4c02") + 
        labs(color="Cluster") + scale_x_discrete(name="Development stage", limits=c("D-shaped", "Veliger", "Umbo", "Eye-spot")) +
        scale_y_continuous(name = "Expression level", limits = c(-5,5)))

# Plot only the genes of cluster 2 :
ggplot(levels[which(levels$cluster==2),], aes(x = num.stade, y = expr.level,fill = Gene, color = as.factor(cluster))) + 
  geom_line(color = "#2ca25f") + 
  labs(color="Cluster") + scale_x_discrete(name="Development stage", limits=c("D-shaped", "Veliger", "Umbo", "Eye-spot")) +
  scale_y_continuous(name = "Expression level", limits = c(-5,5))

# Plot only the genes of cluster 3 :
ggplot(levels[which(levels$cluster==3),], aes(x = num.stade, y = expr.level,fill = Gene, color = as.factor(cluster))) + 
  geom_line(color = "#0570b0") + 
  labs(color="Cluster") + scale_x_discrete(name="Development stage", limits=c("D-shaped", "Veliger", "Umbo", "Eye-spot")) +
  scale_y_continuous(name = "Expression level", limits = c(-5,5))



############################################
## Mise en forme pour tableau format html ##
############################################

colnames(tmp)[7]="Transcripts"
Genes_biomin_Tab = join(tmp, My1944DEGs_annotation)
Genes_biomin_Tab = Genes_biomin_Tab[,c("Transcripts", "Uniprot_Accession", "Uniprot_ID", "Pfucata", "Cgigas", "Best_hit", "cluster")]
rownames(Genes_biomin_Tab) <- NULL
kable(Genes_biomin_Tab, align = 'r') %>%
  kable_styling("striped", "condensed", full_width = F) %>%
  #save_kable(file = "Rplot_table_29GenesBiomin_v1.html", self_contained = T)
  #save_kable(file = "Rplot_table_65GenesImmuneSystem_v1.html", self_contained = T)
  save_kable(file = "Rplot_table_65Genesperception_v1.html", self_contained = T)

Genes_biomin_Tab = cbind(tmp, My1944DEGs_annotation)
Genes_biomin_Tab = Genes_biomin_Tab[,c("Transcripts", "Uniprot_Accession", "Uniprot_ID", "cluster")]
rownames(Genes_biomin_Tab) <- NULL
kable(Genes_biomin_Tab, align = 'r') %>%
  kable_styling("striped", "condensed", full_width = F) %>%
  #save_kable(file = "Rplot_table_29GenesBiomin_v2.html", self_contained = T)
  #save_kable(file = "Rplot_table_65GenesImmuneSystem_v2.html", self_contained = T)
  save_kable(file = "Rplot_table_65Genesperception_v2.html", self_contained = T)

rownames(tmp)=NULL
colnames(tmp)[5]="Transcripts"
rownames(Genes_Biomin_selected)=NULL
Genes_Biomin_TableFull = join(tmp, My1944DEGs_annotation)
kable(Genes_Biomin_TableFull, align = 'r') %>%
  kable_styling("striped", "condensed", full_width = F) %>%
  #save_kable(file = "Rplot_table_29GenesBiomin_v3.html", self_contained = T)
  #save_kable(file = "Rplot_table_65GenesImmuneSystem_v3.html", self_contained = T)
  save_kable(file = "Rplot_table_65Genesperception_v3.html", self_contained = T)


# Work on cluster 3 only. Genes that decrease over stages.
#rownames(mycluster_biomin)=NULL
#mycluster_biomin_cluster3 = mycluster_biomin[which(mycluster_biomin$cluster==3),]
#mycluster_biomin_cluster3 = Genes_biomin_Tab[which(Genes_biomin_Tab$Transcripts %in% mycluster_biomin_cluster3$Gene),]


























########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################


                                                    #############################
                                                    ######   ZONE DE TEST #######
                                                    #############################


########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################



#############################################
###  Reassemble 2 tableaux d informations ###
#############################################

uniprotxlsx  = read.xlsx(file = "uniprot-318Genes_Biomin.xlsx", 1, stringsAsFactors = F, check.names = F, header = TRUE)
colnames(uniprotxlsx)[1]="Uniprot_Accession"

## Mise en forme pour tableau format html
colnames(tmp)[7]="Transcripts"
Genes_express_info = join(tmp, My1944DEGs_annotation)

Genes_express_info$Uniprot_Accession = as.character(Genes_express_info$Uniprot_Accession)
Genes_express_info$Uniprot_ID = as.character(Genes_express_info$Uniprot_ID)
Genes_express_info$Pfucata = as.character(Genes_express_info$Pfucata)
Genes_express_info = Genes_express_info[, c("D-shaped", "Velyger", "Umbo", "Eye-spot", "Genes", "cluster", "Transcripts", "Uniprot_Accession", "Uniprot_ID", "Pfucata", "Cgigas", "Best_hit")]

Uniprot_Genes_express_info = join(Genes_express_info, uniprotxlsx, by="Uniprot_Accession")
write.xlsx2(Uniprot_Genes_express_info, "Uniprot_Genes_Biomin_Tab.xlsx")




#################################
# Work on gene list from Jeremy #
#################################

# retrieve info from annotation.txt with a list of genes.
res_annotation  = read.table("annotation.txt", sep="\t", header = T, na.strings = "NA")
My318DEG_biomin = read.table("Fichier_GO_plus_genes_tmp.txt", sep="\t", header = F, na.strings = "NA")
My1944DEGs      = data.frame(rownames(AllCondlog2FCsup2))
colnames(My318DEG_biomin) = "Transcripts"
colnames(My1944DEGs)      = "Transcripts"

My1944DEGs_complete = res_annotation[which(res_annotation$Transcripts %in% My1944DEGs$Transcripts),]
My1944DEGs_complete$Cgigas = gsub("\\|", "",gsub(".*ref|", "", My1944DEGs_complete$Cgigas))
My1944DEGs_complete$Best_hit= gsub("\\|", "",gsub(".*ref|", "", My1944DEGs_complete$Best_hit))
df=data.frame(t(data.frame(strsplit(as.character(My1944DEGs_complete$Sprot), '\\|'))))
tmp=cbind(My1944DEGs_complete, df)
colnames(tmp)[6:8]=c("sp", "Uniprot_Accession", "Uniprot_ID")
tmp = tmp[,c("Transcripts", "Uniprot_Accession", "Uniprot_ID", "Pfucata", "Cgigas",  "Best_hit")]
tmp$Biomin_gene = 0 # Add a new
tmp[which(tmp$Transcripts %in% My318DEG_biomin$Transcripts),"Biomin_gene"]=1

Genes_biomin_Jeremy = read.delim(file = "Gene_Biomin_Jeremy.txt", stringsAsFactors = F, check.names = F, header = F)
colnames(Genes_biomin_Jeremy)= "Transcripts"
Genes_biomin_Jeremy$Biomin_gene_jeremy = 1
Genes_biomin_1944DEGs_withJeremy = join(tmp, Genes_biomin_Jeremy)
Genes_biomin_1944DEGs_withJeremy[which(is.na(Genes_biomin_1944DEGs_withJeremy$Biomin_gene_jeremy)), "Biomin_gene_jeremy"] = 0
Genes_biomin_1944DEGs_withJeremy[which(Genes_biomin_1944DEGs_withJeremy$Biomin_gene_jeremy==1),]





###################################################################################################################


res = read.delim(file = "clipboard", stringsAsFactors = F, check.names = F, header = TRUE)
colnames(res)=c("Regulation", "D-shaped / Veliger", "Veliger / Umbo", "Umbo / Eye-spot")
levels = res %>% gather(key = "stade", value = "DEGs", -Regulation) %>% mutate(num.stade = case_when(stade == "D-shaped / Veliger" ~ 1,stade == "Veliger / Umbo" ~ 2,stade == "Umbo / Eye-spot" ~ 3))
ggplot(data=levels, aes(x=num.stade, y=DEGs, fill=Regulation)) + geom_bar(stat="identity", position=position_dodge()) + 
  geom_text(aes(label=DEGs), vjust=-0.5, position = position_dodge(0.9), size=3.5) +
  scale_x_discrete(name="Larval stages", limits = c("D-shaped / Veliger", "Veliger / Umbo", "Umbo / Eye-spot")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(name = "Differential expression genes")



