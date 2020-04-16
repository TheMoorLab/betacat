library("Rsubread")
library("edgeR")
library("dplyr")
library(RColorBrewer) 
library(wesanderson)
library(devtools)
library(Seurat)
library(fgsea)
library(msigdbr)
library(BiocManager)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(ggraph)
library(cowplot)
library(SeuratWrappers)
library(conos)
library(SeuratData)
library('psupertime')
library('SingleCellExperiment')
library(magrittr)
library(data.table)
library(stringr)
library(ggrepel)

###########IMPORT FILES##########
#raw bam files downloaded from SUSHI after running AtacENCODEApp pipeline (Davis et al, 2018), which includes: 
# - adapter trimming with cutadapt
# - alignment on GRCh38.95 using Bowtie2 with default mapping parameters (Landmead and Salzberg, 2012),
# - duplicate filtering with Picard 
# - peak calling with MACS2 (p<0.01)

#The peak files were merged in BEDtools v2.29.2 (Quinlan and Hall, 2010) and converted into a saf annotation (samples.saf)
# 1. converting narrowPeak files (macs2 with 0.01 pvalue output) to bed "cut -f 1-6 ATAC4.narrowPeak > ATAC4.bed"
# 2. concatenate "cat ATAC1.bed ATAC2.bed ATAC3.bed ATAC4.bed ATAC5.bed ATAC6.bed > samples.bed"
# 3. sort by chromosome and starting position "sort -k1,1 -k2,2n samples.bed > samples.sorted.bed"
# 4. merge "bedtools merge -i samples.sorted.bed > samples.merged.bed"
# 5. convert to saf "awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' samples.merged.bed > samples.saf"

###########EDGER##########
#nraw reads were mapped to the saf annotation with the FeatureCounts function in the package rsubread (Liao et al, 2019)
all_files = featureCounts(files = c( "control1.bam",
                                     "control2.bam",
                                     "control3.bam",
                                     "D164A1.bam",
                                     "D164A2.bam",
                                     "D164A3.bam"), 
                          allowMultiOverlap=TRUE, 
                          isPairedEnd = T, 
                          annot.ext="samples.saf")

head(all_files$counts)
colnames(all_files$counts) = c("control_1", "control_2", "control_3", "D164A_1", "D164A_2", "D164A_3")
head(all_files$counts)
all_counts = as.matrix(all_files$counts)
head(all_counts)
write.table(all_counts, file = "rawCounts_all.txt") 

#Create DGEList object
group <- factor(c(1,1,1, 2,2,2))
y     <- DGEList(counts=all_counts,group=group)
y$samples
barplot(y$samples$lib.size,names=colnames(y),las=2, main = "Barplot of library sizes")

# filter count matrix for peaks with low coverage (minimum read count of 10 for each sample)
keep    <- filterByExpr(y, min.count = 10)
y       <- y[keep, , keep.lib.sizes=FALSE]

#normalize the library and estimate dispersion
y       <- calcNormFactors(y)
y$samples
design  <- model.matrix(~group)
y       <- estimateDisp(y,design)

#plot MDS
plotMDS(y)
logcpm  <- edgeR::cpm(y, log=TRUE)
head(logcpm)

# perform exact test for differential gene expression
et2v1   <- exactTest(y, pair=c("1","2"))
topTags(et2v1)

plotMD(et2v1)
summary(decideTests(et2v1))

# Export all peaks
all_et2v1  = topTags(et2v1, n = Inf)
dim(all_et2v1)
head(all_et2v1$table, 100)
write.table(all_et2v1$table, file = "controlvsD164A.txt")
r = rownames(all_et2v1$table)
all_et2v1$table$GeneID = r
head(all_et2v1$table)

###########ANNOTATION IN HOMER##########
#merge the edgeR output to the peak genomic locations (samples.saf), to generate a pre-BED file 
saf = read.csv("samples.saf", sep = "\t", header = F)
head(saf)
colnames(saf) = c("GeneID", "V2", "V3", "V4", "V5")

#prepare for Homer peak annotation
bed_controlvsD164A = merge(all_et2v1$table, saf, by = "GeneID")
head(bed_controlvsD164A)

## re-arrange the columns according the bed format either after exporting or in R
controlvsD164A_bed <- bed_controlvsD164A[, c("V2", "V3", "V4","V5", "GeneID","logFC", "logCPM", "PValue", "FDR")]
head(controlvsD164A_bed)
colnames(controlvsD164A_bed) <- c('chrom', 'chromStart', 'chromEnd', 'strand', 'GeneID', "logFC", "logCPM", "PValue", "FDR")
write.table(controlvsD164A_bed, file = "controlvsD164A.txt", quote=FALSE, sep="\t", row.names=F)

#add PeakID column in Excel, export as tab delimited file, then annotate peaks in HOMER:
  # 1) changeNewLine.pl controlvsD164A.txt  
  # 2) annotatePeaks.pl controlvsD164A.txt mm10 > annotation.txt (previously installed mm10)
#import annotation file: cp /IMCR_shares/Moorlab/Coco/ATACSeq/annotation.txt /home/moorlab/mnt_rstudio/Coco/BetacatProject/ATACSeq

annotation <- read.delim("annotation.txt")
View(annotation)

#merge by GeneID with result of EdgeR: now you have annotated and differentially accessible peaks
annotateddiffpeaks = merge(all_et2v1$table, annotation, by = "GeneID")
View(annotateddiffpeaks)
nrow(annotateddiffpeaks)

#plot
ggplot(annotateddiffpeaks, aes(x=logFC, y=-log10(PValue), label = Gene.Name)) + geom_text()

#####VOLCANO PLOT######
annotateddiffpeaks$Significant <- ifelse((annotateddiffpeaks$PValue < 0.01 & abs(annotateddiffpeaks$logFC) > 1 ), "Significant", "Not Sig")
ggplot(data=annotateddiffpeaks, aes(x=logFC, y=-log10(PValue), fill=factor(Significant))) + 
  theme_classic()+ 
  xlim(-4,4)+
  ylim(0, 9)+
  annotate("rect", xmin = -4, xmax = -1, ymin = 2, ymax = 9, alpha = 1, fill= "gray72")+
  annotate("rect", xmin = 1, xmax = 4, ymin = 2, ymax = 9, alpha = .5, fill= "brown")+
  annotate("text", x = -2.5, y = 8.5, label = "lost peaks", size=6) +
  annotate("text", x = 2.5, y = 8.5, label = "gained peaks", size=6) +
  geom_point(shape = 21,color="black")+
  scale_fill_manual(values = c("grey", "darkred")) +
  xlab("log2 fold change") + expand_limits(x = 4)+
  ylab("-log10 adjusted p-value") + 
  theme(legend.position = "none", text = element_text(size=20),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))+
  annotate("point", x = -3.263947728, y = -log10(0.0005581391), shape = 21,color="black", fill="firebrick3", size=4)+
  annotate("label", x = -3.2639477280-.6, y = -log10(0.0005581391), label = "Sp5", size=7) +
  annotate("point", x = -2.259145681, y = -log10(0.002873867), shape = 21,color="black", fill="firebrick3", size=4)+
  annotate("label", x = -2.259145681-0.6, y = -log10(0.002873867), label = "Fzd7", size=7) +
  annotate("point", x = -1.2501949, y = -log10(0.008376697), shape = 21,color="black", fill="firebrick3", size=4)+
  annotate("label", x = -1.2501949+0.6, y = -log10(0.008376697), label = "Axin2", size=7)+
  annotate("point", x = -1.8253572, y = -log10(0.003506013), shape = 21,color="black", fill="firebrick3", size=4)+
  annotate("label", x = -1.8253572+0.6, y = -log10(0.003506013)+.2, label = "Cenpf", size=7)+
  annotate("point", x = -1.8015801, y = -log10(0.006348062), shape = 21,color="black", fill="firebrick3", size=4)+
  annotate("label", x = -1.8015801-0.6, y = -log10(0.006348062)-.4, label = "Cebpe", size=7)+
  annotate("point", x = -1.8829891, y = -log10(0.0004397144), shape = 21,color="black", fill="firebrick3", size=4)+
  annotate("label", x = -1.8829891-0.6, y = -log10(0.0004397144), label = "E2f1", size=7)+
  annotate("point", x = 2.559396, y = -log10(6.514718e-05), shape = 21,color="black", fill="firebrick3", size=4)+
  annotate("label", x = 2.559396+1, y = -log10(6.514718e-05), label = "Mapk8/JNK", size=7)+
  annotate("point", x = 1.386477, y = -log10(3.499661e-03), shape = 21,color="black", fill="firebrick3", size=4)+
  annotate("label", x = 1.386477+.8, y = -log10(3.499661e-03)-.4, label = "Mapk13/SAPK4", size=7)+
  annotate("point", x = 1.151880, y = -log10(0.000231533), shape = 21,color="black", fill="firebrick3", size=4)+
  annotate("label", x = 1.151880+.35, y = -log10(0.000231533), label = "Fos", size=7)+
  annotate("point", x = 1.341831, y = -log10(0.002943785), shape = 21,color="black", fill="firebrick3", size=4)+
  annotate("label", x = 1.341831+.6, y = -log10(0.002943785)+0.3, label = "Muc13", size=7) #is JNK target: https://cancerres.aacrjournals.org/content/69/3/765.long
  
#select differentially accessible peaks and export for HOMER motif analysis ./bin/findMotifsGenome.pl Diffpeaks.txt mm10 HOMER/ 
positions         <- which(abs(annotateddiffpeaks$logFC) > 1 & annotateddiffpeaks$PValue < 0.01)
DiffpeaksHOMER    <- (annotateddiffpeaks[positions, c("logFC", "PValue", "Gene.Name")])
write.table(DiffpeaksHOMER, file = "DiffpeaksHOMER.txt", quote=FALSE, sep="\t", row.names=F)

positionsUP       <- which(annotateddiffpeaks$logFC > 1 & annotateddiffpeaks$PValue < 0.01)
DiffpeaksHOMERUP  <- annotateddiffpeaks[positionsUP, c(6:10)]
write.table(DiffpeaksHOMERUP, file = "DiffpeaksHOMERUP.txt", quote=FALSE, sep="\t", row.names=F)

positionsDN       <- which(annotateddiffpeaks$logFC < -1 & annotateddiffpeaks$PValue < 0.01)
DiffpeaksHOMERDN  <- annotateddiffpeaks[positionsDN, c(6:10)]
write.table(DiffpeaksHOMERDN, file = "DiffpeaksHOMERDN.txt", quote=FALSE, sep="\t", row.names=F)

positionsUP       <- which(annotateddiffpeaks$logFC > 1 & annotateddiffpeaks$PValue < 0.01)
DiffpeaksUP       <- annotateddiffpeaks[positionsUP, c("GeneID","Gene.Name","Detailed.Annotation", "logFC", "PValue")]
View(DiffpeaksUP)

positionsDN       <- which(annotateddiffpeaks$logFC < -1 & annotateddiffpeaks$PValue < 0.01)
DiffpeaksDN       <- annotateddiffpeaks[positionsDN, c("GeneID","Gene.Name","Detailed.Annotation", "logFC", "PValue")]
nrow(DiffpeaksDN)

#find gained peaks under control of AP1 complex TFTs
#command line HOMER: ./bin/findMotifsGenome.pl DiffpeaksHOMERUP.txt mm10 HOMER/ -find ./bin/motifs/ap1.motif > AP1.txt
AP1      <- read.delim("~/mnt_rstudio/Coco/BetacatProject/ATACSeq/AP1.txt")
fra1    <- read.delim("~/mnt_rstudio/Coco/BetacatProject/ATACSeq/fra1.txt")
fra2    <- read.delim("~/mnt_rstudio/Coco/BetacatProject/ATACSeq/fra2.txt")
fos     <- read.delim("~/mnt_rstudio/Coco/BetacatProject/ATACSeq/fos.txt")
atf3    <- read.delim("~/mnt_rstudio/Coco/BetacatProject/ATACSeq/atf3.txt")
junb    <- read.delim("~/mnt_rstudio/Coco/BetacatProject/ATACSeq/junb.txt")
batf    <- read.delim("~/mnt_rstudio/Coco/BetacatProject/ATACSeq/batf.txt")
junap1  <- read.delim("~/mnt_rstudio/Coco/BetacatProject/ATACSeq/junap1.txt")
spdef   <- read.delim("~/mnt_rstudio/Coco/BetacatProject/ATACSeq/spdef.txt")

#merge genomic position with annotation and edgeR result
annotateddiffpeaks$PositionID <- tolower(annotateddiffpeaks$PositionID)
AP1targets    = merge(annotateddiffpeaks[positionsUP, c("PositionID", "Gene.Name", "logFC", "PValue", "GeneID")], AP1, by = "PositionID")
fra1targets   = merge(annotateddiffpeaks[positionsUP, c("PositionID", "Gene.Name", "logFC", "PValue", "GeneID")], fra1, by = "PositionID")
fra2targets   = merge(annotateddiffpeaks[positionsUP, c("PositionID", "Gene.Name", "logFC", "PValue", "GeneID")], fra2, by = "PositionID")
fostargets    = merge(annotateddiffpeaks[positionsUP, c("PositionID", "Gene.Name", "logFC", "PValue", "GeneID")], fos, by = "PositionID")
atf3targets   = merge(annotateddiffpeaks[positionsUP, c("PositionID", "Gene.Name", "logFC", "PValue", "GeneID")], atf3, by = "PositionID")
junbtargets   = merge(annotateddiffpeaks[positionsUP, c("PositionID", "Gene.Name", "logFC", "PValue", "GeneID")], junb, by = "PositionID")
batftargets   = merge(annotateddiffpeaks[positionsUP, c("PositionID", "Gene.Name", "logFC", "PValue", "GeneID")], batf, by = "PositionID")
junap1targets = merge(annotateddiffpeaks[positionsUP, c("PositionID", "Gene.Name", "logFC", "PValue", "GeneID")], junap1, by = "PositionID")
spdeftargets  = merge(annotateddiffpeaks[positionsUP, c("PositionID", "Gene.Name", "logFC", "PValue", "GeneID")], spdef, by = "PositionID")

#get unique gene names
ap1peaks    <- as.vector(AP1targets$Gene.Name[!duplicated(AP1targets$Gene.Name)])
fra1peaks   <- as.vector(fra1targets$Gene.Name[!duplicated(fra1targets$Gene.Name)])
fra2peaks   <- as.vector( fra2targets$Gene.Name[!duplicated(fra2targets$Gene.Name)])
fospeaks    <- as.vector( fostargets$Gene.Name[!duplicated(fostargets$Gene.Name)])
atf3peaks   <- as.vector(atf3targets$Gene.Name[!duplicated(atf3targets$Gene.Name)])
junbpeaks   <- as.vector( junbtargets$Gene.Name[!duplicated(junbtargets$Gene.Name)])
batfpeaks   <- as.vector(batftargets$Gene.Name[!duplicated(batftargets$Gene.Name)])
junap1peaks <- as.vector(junap1targets$Gene.Name[!duplicated(junap1targets$Gene.Name)])

FosJunAp1targetgenesNames       <- c(ap1peaks, fra1peaks, fra2peaks, fospeaks, atf3peaks, junbpeaks, batfpeaks, junap1peaks)
FosJunAp1targetgenes            <- as.data.frame(FosJunAp1targetgenesNames)
colnames(FosJunAp1targetgenes)  <- "Gene.Name"
FosJunAp1targetgenes            = merge(annotateddiffpeaks[positionsUP, c("Gene.Name", "logFC", "PValue")], 
                                        FosJunAp1targetgenes, by = "Gene.Name")
write.table(FosJunAp1targetgenes, file = "FosJunAp1targetgenesNames.txt",  row.names=F, quote=FALSE, sep ="\t")

##########GSEA#########
#select species and set
msigdbr_show_species()
m_df<- msigdbr(species = "Mus musculus", category = "H")
Hallmarks <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(Hallmarks)
m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
CGP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(CGP)
m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
KEGG <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(KEGG)
m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
REACTOME <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(REACTOME)
m_df<- msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT")
TFT <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(TFT)
m_df<- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
BP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(BP)
m_df<- msigdbr(species = "Mus musculus", category = "C6")
oncSig <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(oncSig)

#on all Diffpeaks
ranks <- Diffpeaks %>% 
  na.omit()%>%
  mutate(ranking=-log10(PValue)/sign(logFC))
ranks <- ranks$ranking
names(ranks) <- Diffpeaks$Gene.Name
head(ranks, 10)

BP_Diff <- fgsea(pathways = BP, 
                      stats = ranks,
                      minSize=10,
                      maxSize=500,
                      nperm=1000000)

Hallmarks_Diff %>% filter(abs(NES)>1 & pval<0.05)
View(CGP_Diff %>% filter(abs(NES)>1 & pval<0.05))
View(KEGG_Diff) 
REACTOME_Diff %>% filter(abs(NES)>1 & pval<0.05)
TFT_Diff %>% filter(abs(NES)>1 & pval<0.05) #no sign padj
BP_Diff  %>% filter(abs(NES)>1 & pval<0.05)


#########OVERLAP WITH TRANSCRIPTOMICS######
enrichmentscRNASeq <- rbind(Hallmarks_0vs4, CGP_0vs4, KEGG_0vs4,REACTOME_0vs4, TFT_0vs4, BP_0vs4)
View(enrichmentscRNASeq)
length(enrichmentscRNASeq$pathway)

enrichmentATAC <- rbind(Hallmarks_Diff, CGP_Diff, KEGG_Diff,REACTOME_Diff, TFT_Diff, BP_Diff)
length(enrichmentATAC$pathway)

overrepresentedRNA<- which(enrichmentscRNASeq$NES>1)
OVER_RNA <-enrichmentscRNASeq[overrepresentedRNA,]
underrepresentedRNA<- which(enrichmentscRNASeq$NES < -1)
UNDER_RNA <- enrichmentscRNASeq[underrepresentedRNA,]

overrepresentedATAC<- which(enrichmentATAC$NES>1 )
OVER_ATAC <-enrichmentATAC[overrepresentedATAC,]
underrepresentedATAC<- which(enrichmentATAC$NES < -1)
UNDER_ATAC <- enrichmentATAC[underrepresentedATAC,]

#hypergeometric test
set1 <-  OVER_ATAC$pathway
set2 <- OVER_RNA$pathway
overlap <- intersect(set1, set2)
allterms <- union(enrichmentscRNASeq$pathway, enrichmentATAC$pathway)
phyper(length(overlap), length(set1), length(allterms)-length(set1), length(set2), lower.tail=F)

#upset plot
Genesets <- list("enriched in ATAC 2d pi" = OVER_ATAC$pathway,
                #"depleted in ATAC 2d pi" =UNDER_ATAC$pathway,
                "enriched in scRNAseq 4d pi" =OVER_RNA$pathway)
                #"depleted in scRNASeq 4d pi" = UNDER_RNA$pathway)

upset(fromList(Genesets), nintersects = 8,  main.bar.color = "darkblue", 
      keep.order = T , matrix.color = "darkred", mainbar.y.label = "intersecting GSEA pathways", sets.x.label= "GSEA pathways",
      sets.bar.color = "lightgrey", point.size = 10 , text.scale = 2)

