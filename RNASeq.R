library(tidyverse)
library(pheatmap)
library(ggplot2)
library(plyr)
library(devtools)
install_github("raivokolde/pheatmap")

#import EdgeR files downloaded from Shushi DE=differential expression over het control
DE_D164A    <- read.delim("~/mnt_rstudio/Coco/BetacatProject/RNASeq/result--D164A--over--het.txt", header=T)
DE_DeltaC   <- read.delim("~/mnt_rstudio/Coco/BetacatProject/RNASeq/result--DeltaC--over--het.txt", header=T)
DE_dm       <- read.delim("~/mnt_rstudio/Coco/BetacatProject/RNASeq/result--dm--over--het.txt", header=T)
DE_KO       <- read.delim("~/mnt_rstudio/Coco/BetacatProject/RNASeq/result--floxflox--over--het.txt", header=T)

####PCA######
#extract FPKM values from one of the EdgeR result dataframes (here DE_KO) (TV are sample names from the sequencing)
all_genes   <-  data.frame("Gene"=DE_KO[,"gene_name"],
                             "het_1"=DE_KO[, "TV6..FPKM."],
                             "het_2"=DE_KO[, "TV7..FPKM."],
                             "het_3"=DE_KO[, "TV8..FPKM."],
                             "het_4"=DE_KO[, "TV10..FPKM."],
                             "floxflox_1"=DE_KO[, "TV21..FPKM."],
                             "floxflox_2"=DE_KO[, "TV23..FPKM."],
                             "dm_1"=DE_KO[, "TV18..FPKM."],
                             "dm_2"=DE_KO[, "TV19..FPKM."],
                             "dm_3"=DE_KO[, "TV20..FPKM."],
                             "DeltaC_1"=DE_KO[, "TV13..FPKM."],
                             "DeltaC_2"=DE_KO[, "TV14..FPKM."],
                             "D164A_1"=DE_KO[, "TV15..FPKM."],
                             "D164A_2"=DE_KO[, "TV16..FPKM."],
                             "D164A_3"=DE_KO[, "TV17..FPKM."], stringsAsFactors = T)
all_genes[1] <- NULL

#average FPKM all genes
het     <- rowMeans(all_genes[,1:4])
KO      <- rowMeans(all_genes[,5:6])
dm      <- rowMeans(all_genes[,7:9])
DeltaC  <- rowMeans(all_genes[,10:11])
D164A   <- rowMeans(all_genes[,12:14])
all_genes_average           <- cbind(het, KO, dm, DeltaC, D164A)
colnames(all_genes_average) <- c("control", "KO", "dm", "DeltaC", "D164A")

#pca
data_for_PCA                <- t(all_genes_average)
dim(data_for_PCA)

# calculate MDS (matrix of dissimilarities)
mds <- cmdscale(dist(data_for_PCA), k=3, eig=TRUE)  

#The variable mds$eig provides the Eigen values for the first 8 principal components:
mds$eig

# transform the Eigen values into percentage
eig_pc <- mds$eig * 100 / sum(mds$eig)

# plot the PCA
barplot(eig_pc,
        las=1,
        xlab="Dimensions", 
        ylab="Proportion of explained variance (%)", y.axis=NULL,
        col="darkgrey")

pca = prcomp(all_genes_average, center=T)
plot(pca$rotation[,1],pca$rotation[,2], xlab = "PC1 (94.8%)", ylab = "PC2 (3.0%)",
     xlim = c(-0.72,-0.3), ylim = c(-0.8,0.6),
     cex=2, pch=21, bg="darkred", col="black")
text(pca$rotation[,1],pca$rotation[,2], row.names(pca$rotation), cex=1, pos=2)
summary(pca)

#########SIGNIFICANT DEGS############
threshold_pvalue = 0.05
threshold_log2fc = 2

positions     <- which(DE_D164A$pValue < threshold_pvalue & (DE_D164A$log2.Ratio < -threshold_log2fc | DE_D164A$log2.Ratio > threshold_log2fc))
sigDEG_D164A  <- data.frame("Gene"=DE_D164A[positions,"gene_name"], "log2ratio"=DE_D164A[positions,"log2.Ratio"], "pValue"=DE_D164A[positions,"pValue"], stringsAsFactors = FALSE)
sigDEG_D164A  %>% remove_rownames %>% column_to_rownames(var="Gene")

positions     <- which(DE_DeltaC$pValue < threshold_pvalue & (DE_DeltaC$log2.Ratio < -threshold_log2fc | DE_DeltaC$log2.Ratio > threshold_log2fc))
sigDEG_DeltaC <- data.frame("Gene"=DE_DeltaC[positions,"gene_name"], "log2ratio"=DE_DeltaC[positions,"log2.Ratio"], "pValue"=DE_DeltaC[positions,"pValue"], stringsAsFactors = FALSE)
sigDEG_DeltaC %>% remove_rownames %>% column_to_rownames(var="Gene")

positions     <- which(DE_dm$pValue < threshold_pvalue & (DE_dm$log2.Ratio < -threshold_log2fc | DE_dm$log2.Ratio > threshold_log2fc))
sigDEG_dm     <- data.frame("Gene"=DE_dm[positions,"gene_name"], "log2ratio"=DE_dm[positions,"log2.Ratio"], "pValue"=DE_dm[positions,"pValue"], stringsAsFactors = FALSE)
sigDEG_dm     %>% remove_rownames %>% column_to_rownames(var="Gene")

positions     <- which(DE_KO$pValue < threshold_pvalue & (DE_KO$log2.Ratio < -threshold_log2fc | DE_KO$log2.Ratio > threshold_log2fc))
sigDEG_KO     <- data.frame("Gene"=DE_KO[positions,"gene_name"],"log2ratio"=DE_KO[positions,"log2.Ratio"], "pValue"=DE_KO[positions,"pValue"], stringsAsFactors = FALSE)
sigDEG_KO     %>% remove_rownames %>% column_to_rownames(var="Gene")

#use these gene lists for EnrichR

#############UPSET PLOT#################
library(UpSetR)
sigDEGs   <- list("KO" = sigDEG_KO$Gene,
                  "dm" = sigDEG_dm$Gene,
                  "DeltaC" = sigDEG_DeltaC$Gene,
                   "D164A" = sigDEG_D164A$Gene)

upset(fromList(sigDEGs), 
      order.by = "freq", 
      nintersects = NA,  
      main.bar.color = "darkblue", 
      keep.order = T, 
      matrix.color = "darkred", 
      mainbar.y.label = "number of intersecting DEGs", 
      sets.x.label= "number of DEGs",
      sets.bar.color = "lightgrey", 
      point.size = 10, 
      empty.intersections=T, 
      text.scale = 4)
 
 
#############MARKER HEATMAP#################                  
#extract FPKM value for marker gene heatmap
markergenes <- data.frame("Axin2", "Myc", "Cd44", "Mmp7", "Sox9", "Fzd2",  #Wnt target
                          "Lgr5", "Olfm4", "Ascl2", "Ephb2", "Ephb3", "Gemin4", "Wdr12", "Tnfrsf19", "Nr2e3", "Slc14a1",  #Stem cell markers
                          "Mki67","Foxm1", "Ccnd1", "Ccnd2", #proliferation
                          stringsAsFactors = FALSE)

markergenesFPKM <- data.frame("Gene"=NA,"floxflox_1"=NA, "floxflox_2"=NA, 
                              "dm_1"=NA, "dm_2"=NA, "dm_3"=NA,  "DeltaC_1"=NA,  "DeltaC_2"= NA,   "D164A_1"=NA ,  "D164A_2"= NA,   "D164A_3"=NA)
row=1
for (i in 1:ncol(markergenes))
{
  name=markergenes[1,i]
  match1 <- which(DE_KO$gene_name == name)
  markergenesFPKM[row,] <- c(name, DE_KO[match1,"TV21..FPKM."], 
                             DE_KO[match1,"TV23..FPKM."], 
                             DE_KO[match1,"TV18..FPKM."],  
                             DE_KO[match1,"TV19..FPKM."], 
                             DE_KO[match1,"TV20..FPKM."], 
                             DE_KO[match1,"TV13..FPKM."], 
                             DE_KO[match1,"TV14..FPKM."], 
                             DE_KO[match1,"TV15..FPKM."],
                             DE_KO[match1,"TV16..FPKM."],
                             DE_KO[match1,"TV17..FPKM."], stringsAsFactors = T)
  row=row+1
}

#make numeric
markergenesFPKMnum            <-  data.matrix(markergenesFPKM[, c(2:11)], rownames.force = NA)
rownames(markergenesFPKMnum)  <- markergenesFPKM$Gene
colnames(markergenesFPKMnum)  <- c("KO_1", "KO_2","dm_1", "dm_2", "dm_3",  "DeltaC_1",  "DeltaC_2",   "D164A_1",  "D164A_2",   "D164A_3")

#prepare palette for pheatmap
paletteLength   <- 50
myColor         <- colorRampPalette(c("blue", "white", "darkorange"))(paletteLength)
breaksList      = seq(-2, 2, by = 0.04)

#prepare annotation for pheatmap
annotation_rows             <- data.frame(markers = rep(c("Wnt signalling", "Stem cells", "Proliferation"), c(6,10,4)))
rownames(annotation_rows)   <- rownames(markergenesFPKMnum)
annotation_rows$markers     <- factor(annotation_rows$markers, levels = c("Wnt signalling", "Stem cells", "Proliferation"))

pheatmap(markergenesFPKMnum, scale="row", 
         color = colorRampPalette(c("blue", "white", "darkorange"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         cluster_rows = F, cluster_cols = T, 
         border_color = "black", 
         legend_breaks = -2:2, 
         cellwidth = 40, cellheight = 30,
         angle_col = "45", 
         annotation_row = annotation_rows,
         fontsize = 20)

###########################GSEA##############
library(fgsea)
library(msigdbr)
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)

###GSEA ON DIFFERENTIAL EXPRESSION ##
#select species and set
msigdbr_show_species()
m_df      <- msigdbr(species = "Mus musculus", category = "H")
Hallmarks <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
m_df      <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
CGP       <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
m_df      <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
KEGG      <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
m_df      <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
REACTOME  <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
m_df      <- msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT")
TFT       <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
m_df      <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
BP        <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
m_df      <- msigdbr(species = "Mus musculus", category = "C6")
oncSig    <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#on significant DEGs in KO
ranks <- sigDEG_KO %>% 
  na.omit()%>%
  mutate(ranking=-log10(pValue)/sign(log2ratio))
ranks <- ranks$ranking
names(ranks) <- sigDEG_KO$Gene
head(ranks, 10)

Hallmarks_KO <- fgsea(pathways = Hallmarks, 
                 stats = ranks,
                 minSize=10,
                 maxSize=500,
                 nperm=1000000)

KO <- ggplot(Hallmarks_KO %>% filter(abs(NES)>1.5 & pval<0.05) %>% head(n= 100), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=pval, width=3)) +
  scale_size_area(max_size = 10)+
  scale_color_gradientn(colours = brewer.pal(11,"RdYlBu")[1:11]) +
  coord_flip() +
  labs(x=" ", y="Normalized Enrichment Score",
       title=" ", cols="black") + 
  theme_classic()+
  scale_y_reverse(limits = c(-1,-3.5))+
  labs(y="")+   
  theme(axis.text.y=element_text(size=15), axis.text.x=element_text(size=15))

#on significant DEGs in dm
ranks <- sigDEG_dm %>% 
  na.omit()%>%
  mutate(ranking=-log10(pValue)/sign(log2ratio))
ranks <- ranks$ranking
names(ranks) <- sigDEG_dm$Gene
head(ranks, 10)

Hallmarks_dm <- fgsea(pathways = Hallmarks, 
               stats = ranks,
               minSize=10,
               maxSize=500,
               nperm=1000000000)

dm <- ggplot(Hallmarks_dm %>% filter(abs(NES)>1.5 & pval<0.05) %>% head(n= 100), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=pval, width=3)) +
  scale_size_area(max_size = 10)+
  scale_color_gradientn(colours = brewer.pal(11,"RdYlBu")[1:11]) +
  coord_flip() +
  labs(x=" ", y="Normalized Enrichment Score",
       title=" ", cols="black") + 
  theme_classic()+
  scale_y_reverse(limits = c(-1,-3.5))+
  labs(y="")+   
  theme(axis.text.y=element_text(size=15), axis.text.x=element_text(size=15))

#on significant DEGs in DeltaC
ranks <- sigDEG_DeltaC %>% 
  na.omit()%>%
  mutate(ranking=-log10(pValue)/sign(log2ratio))
ranks <- ranks$ranking
names(ranks) <- sigDEG_DeltaC$Gene
head(ranks, 10)

Hallmarks_DeltaC <- fgsea(pathways = Hallmarks, 
               stats = ranks,
               minSize=10,
               maxSize=500,
               nperm=1000000)

View(Hallmarks_DeltaC %>% filter(abs(NES)>1.5 & pval<0.05))

DeltaC <- ggplot(Hallmarks_DeltaC %>% filter(abs(NES)>1.5 & pval<0.05) %>% head(n= 100), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=pval, width=3)) +
  scale_size_area(max_size = 10)+
  scale_color_gradientn(colours = brewer.pal(11,"RdYlBu")[1:11]) +
  coord_flip() +
  labs(x=" ", y="NORMALIZED ENRICHMENT SCORE",
       title=" ", cols="black") + 
  scale_y_reverse(limits = c(-1,-3.5))+
  theme_classic()+  
  theme(axis.text.y=element_text(size=15), axis.text.x=element_text(size=15), axis.title=element_text(size=15, face="bold"))

#make one panel
library(ggpubr)
theme_set(theme_pubr())
figure <- ggarrange(KO, dm, DeltaC,
                    common.legend = TRUE, legend = "bottom",
                    labels = c("enriched in KO", "enriched in dm", "enriched in DeltaC"),
                    font.label =list(size = 17, color = "black", face = "bold"),
                    label.x = c(0.3, 0.3, 0.27),
                    ncol = 1, nrow = 3)
figure

#on significant DEGs in D164A
ranks <- sigDEG_D164A %>% 
  na.omit()%>%
  mutate(ranking=-log10(pValue)/sign(log2ratio))
ranks <- ranks$ranking
names(ranks) <- sigDEG_D164A$Gene
head(ranks, 10)

Hallmarks_D164A <- fgsea(pathways = Hallmarks, 
                         stats = ranks,
                         minSize=10,
                         maxSize=500,
                         nperm=1000000)
View(Hallmarks_D164A %>% filter(abs(NES)>1.5 & pval<0.05))

ggplot(Hallmarks_D164A %>% filter(abs(NES)>1.5 & pval<0.05) %>% head(n= 100), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=pval)) +
  scale_color_gradientn(colours = brewer.pal(11,"RdYlBu")[1:11]) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=" ") + 
  theme_classic()+
  labs(y="")+   
  theme(axis.text.y=element_text(size=30))






