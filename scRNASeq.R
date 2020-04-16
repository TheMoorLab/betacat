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

################PRE-PROCESSING OF INDIVIDUAL DATASETS################
#import rds files after zUMIs pipeline
control_0d  <- readRDS("/home/moorlab/mnt_rstudio/rstudio_projects/Costanza/H2_merged.dgecounts.rds") #heterozygous non injected
control_2d  <- readRDS("/home/moorlab/mnt_rstudio/rstudio_projects/Costanza/H1_merged.dgecounts.rds") #heterozygous 2d pi
control_4d  <- readRDS("/home/moorlab/mnt_rstudio/rstudio_projects/Costanza/H3_merged.dgecounts.rds") #heterozygous 4d pi
D164A_0d    <- readRDS("/home/moorlab/mnt_rstudio/rstudio_projects/Costanza/N2_merged.dgecounts.rds") #D164A non injected
D164A_2d    <- readRDS("/home/moorlab/mnt_rstudio/rstudio_projects/Costanza/N1_merged.dgecounts.rds") #D164A 2d pi
D164A_4d    <- readRDS("/home/moorlab/mnt_rstudio/rstudio_projects/Costanza/N3_merged.dgecounts.rds") #D164A 4d pi

source(Functions_KHandler.R)
control_0da       <- Annotation_mouse(control_0d)
control_0dS       <- CreateSeurat_200_3_mouse(control_0da, "control_0d", 3, 200, c(0.0125, 3), c(0.5, Inf), 30)
control_0d_r      <- Annotation_mouse_reads(control_0d, "control_0d")
control_0d_rS     <- CreateSeurat_200_3_mouse(control_0d_r, "control_0d", 3, 200, c(0.0125, 3), c(0.5, Inf), 30)
data_m            <- as.data.frame(as.matrix(control_0dS@meta.data))
data_reads_m      <- as.data.frame(as.matrix(control_0d_rS@meta.data))
a                 <- merge(data_m, data_reads_m, by = "row.names")
a$reads_per_umi   <- (as.numeric(a$nCount_RNA.y)) / (as.numeric(a$nCount_RNA.x))
a$reads_per_gene  <- (as.numeric(a$nCount_RNA.y)) / (as.numeric(a$nFeature_RNA.x))
median(a$reads_per_gene)
median(a$reads_per_umi)

control_2da       <- Annotation_mouse(control_2d) 
control_2dS       <- CreateSeurat_200_3_mouse(control_2da, "control_2d", 3, 200, c(0.0125, 3), c(0.5, Inf), 30)
control_2d_r      <- Annotation_mouse_reads(control_2d, "control_2d")
control_2d_rS     <- CreateSeurat_200_3_mouse(control_2d_r, "control_2d", 3, 200, c(0.0125, 3), c(0.5, Inf), 30)
data_m            <- as.data.frame(as.matrix(control_2dS@meta.data))
data_reads_m      <- as.data.frame(as.matrix(control_2d_rS@meta.data))
a                 <- merge(data_m, data_reads_m, by = "row.names")
a$reads_per_umi   <- (as.numeric(a$nCount_RNA.y)) / (as.numeric(a$nCount_RNA.x))
a$reads_per_gene  <- (as.numeric(a$nCount_RNA.y)) / (as.numeric(a$nFeature_RNA.x))
median(a$reads_per_gene)
median(a$reads_per_umi)

control_4da        <- Annotation_mouse(control_4d)
control_4dS        <- CreateSeurat_200_3_mouse(control_4da, "control_4d", 3, 200, c(0.0125, 3), c(0.5, Inf), 30)
control_4d_r       <- Annotation_mouse_reads(control_4d, "control_4d")
control_4d_rS      <- CreateSeurat_200_3_mouse(control_4d_r, "control_4d", 3, 200, c(0.0125, 3), c(0.5, Inf), 30)
data_m             <- as.data.frame(as.matrix(control_4dS@meta.data))
data_reads_m       <- as.data.frame(as.matrix(control_4d_rS@meta.data))
a                  <- merge(data_m, data_reads_m, by = "row.names")
a$reads_per_umi    <- (as.numeric(a$nCount_RNA.y)) / (as.numeric(a$nCount_RNA.x))
a$reads_per_gene   <- (as.numeric(a$nCount_RNA.y)) / (as.numeric(a$nFeature_RNA.x))
median(a$reads_per_gene)
median(a$reads_per_umi)

D164A_0da         <- Annotation_mouse(D164A_0d)
D164A_0dS         <- CreateSeurat_200_3_mouse(D164A_0da, "D164A_0d", 3, 200, c(0.0125, 3), c(0.5, Inf), 30)
D164A_0d_r        <- Annotation_mouse_reads(D164A_0d, "D164A_0d")
D164A_0d_rS       <- CreateSeurat_200_3_mouse(D164A_0d_r, "D164A_0d", 3, 200, c(0.0125, 3), c(0.5, Inf), 30)
data_m            <- as.data.frame(as.matrix(D164A_0dS@meta.data))
data_reads_m      <- as.data.frame(as.matrix(D164A_0d_rS@meta.data))
a                 <- merge(data_m, data_reads_m, by = "row.names")
a$reads_per_umi   <- (as.numeric(a$nCount_RNA.y)) / (as.numeric(a$nCount_RNA.x))
a$reads_per_gene  <- (as.numeric(a$nCount_RNA.y)) / (as.numeric(a$nFeature_RNA.x))
median(a$reads_per_gene)
median(a$reads_per_umi)

D164A_2da         <- Annotation_mouse(D164A_2d)
D164A_2dS         <- CreateSeurat_200_3_mouse(D164A_2da, "D164A_2d", 3, 200, c(0.0125, 3), c(0.5, Inf), 30)
D164A_2d_r        <- Annotation_mouse_reads(D164A_2d, "D164A_2d")
D164A_2d_rS       <- CreateSeurat_200_3_mouse(D164A_2d_r, "D164A_2d", 3, 200, c(0.0125, 3), c(0.5, Inf), 30)
data_m            <- as.data.frame(as.matrix(D164A_2dS@meta.data))
data_reads_m      <- as.data.frame(as.matrix(D164A_2d_rS@meta.data))
a                 <- merge(data_m, data_reads_m, by = "row.names")
a$reads_per_umi   <- (as.numeric(a$nCount_RNA.y)) / (as.numeric(a$nCount_RNA.x))
a$reads_per_gene  <- (as.numeric(a$nCount_RNA.y)) / (as.numeric(a$nFeature_RNA.x))
median(a$reads_per_gene)
median(a$reads_per_umi)

D164A_4da         <- Annotation_mouse(D164A_4d)
D164A_4dS         <- CreateSeurat_200_3_mouse(D164A_4da, "D164A_4d", 3, 200, c(0.0125, 3), c(0.5, Inf), 30)
D164A_4d_r        <- Annotation_mouse_reads(D164A_4d, "D164A_4d")
D164A_4d_rS       <- CreateSeurat_200_3_mouse(D164A_4d_r, "D164A_4d", 3, 200, c(0.0125, 3), c(0.5, Inf), 30)
data_m            <- as.data.frame(as.matrix(D164A_4dS@meta.data))
data_reads_m      <- as.data.frame(as.matrix(D164A_4d_rS@meta.data))
a                 <- merge(data_m, data_reads_m, by = "row.names")
a$reads_per_umi   <- (as.numeric(a$nCount_RNA.y)) / (as.numeric(a$nCount_RNA.x))
a$reads_per_gene  <- (as.numeric(a$nCount_RNA.y)) / (as.numeric(a$nFeature_RNA.x))
median(a$reads_per_gene)
median(a$reads_per_umi)

#################MERGE PRE-PROCESSED SAMPLES################
timecourse_merged <- merge(control_ni_rS, c(control_2d_rS, control_4d_rS, D164A_ni_rS, D164A_2d_rS, D164A_4d_rS), add.cell.ids = c("control 0d pi", "control 2d pi", "control 4d pi", "D164A 0d pi", "D164A 2d pi", "D164A 4d pi"))
length(timecourse_merged@active.ident)
FeatureScatter(timecourse_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
timecourse_merged <- subset(timecourse_merged, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito < 0.2)
length(timecourse_merged@active.ident)

#################CONOS###############
#standard Seurat preprocessing (as in Seurat wrapper vignette)
panel_timecourse <- SplitObject(timecourse_merged, split.by = "orig.ident")
for (i in 1:length(panel_timecourse)) {
  panel_timecourse[[i]] <- NormalizeData(panel_timecourse[[i]]) %>% FindVariableFeatures() %>% ScaleData(vars.to.regress = c("nCount_RNA", "percent.mito")) %>% 
    RunPCA(verbose = FALSE, set.seed(1))
}
panel_timecourse.con <- Conos$new(panel_timecourse)

#build joint graph
panel_timecourse.con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30, n.odgenes=2000, 
                     matching.method='mNN', metric='angular', 
                     score.component.variance=TRUE, verbose=TRUE)
plotComponentVariance(panel_timecourse.con, space='PCA')

#leiden clustering
panel_timecourse.con$findCommunities(method=leiden.community, resolution=1.3) 

#UMAP embedding
panel_timecourse.con$embedGraph(method="UMAP", min.dist=0.01, spread=15, n.cores=1, set.seed(1))
umap_dt  = data.table(cell_id=rownames(panel_timecourse.con$embedding), panel.con$embedding)
setnames(umap_dt, names(umap_dt), c("cell_id", "umap1", "umap2"))

#plot
panel_timecourse.con$plotGraph(clustering='leiden', size=1)
panel_timecourse.con$plotGraph(gene="Mki67")

#convert back to seurat object
timecourse <- as.Seurat(panel_timecourse.con, reduction = "UMAP") 
DimPlot(timecourse, group.by = "orig.ident",pt.size = .9)
length(timecourse@active.ident)

#################ANNOTATION OF CLUSTERS###############
markers_timecourse <- FindAllMarkers(object = timecourse, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers_timecourse %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC))

DimPlot(timecourse, group.by = "leiden", pt.size = .9, label = T)
markers.to.plot <-  c("Dclk1", "Rgs13", "Cd24a", #9 = tuft cell markers
                      "Cd3g", "Cd7", "Ccl5", #4 = immune cells 
                      "Gpx1", "Car4",  #6 = immature 
                      "Apoc3", "Sepp1", "Clec2h", #mature
                      "Creb3l3", "Nrli3", #proximal
                      "Olfm4", "Wdr43", "Mki67", "Ccnd1", #stem cells /TA$
                      "Mptx2", "Ang4","Lyz1", #5 = Paneth cells  
                      "Atoh1", "Muc2", "Tff3") #5=goblet cells)
#Dotplot
plot <- DotPlot(timecourse, features = markers.to.plot)
plot + 
  theme(axis.text.x = element_text(angle = 45, face="bold", hjust=1), axis.text.y = element_text(face="italic")) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ theme(legend.position="bottom")+
  coord_flip()

#annotation
current.cluster.ids <- c(1,2,3,4,5,6,7,8, 9, 10,11, 12, 13)
new.cluster.ids <- c( "Early progenitors",
                      "Goblet cells",
                      "Enterocytes 1",
                      "Paneth cells",
                      "Stem cells",
                      "Enterocytes 2",
                      "Enterocytes 2",
                      "Enterocytes 3",
                      "Enterocytes (immature)",
                      "T cells",
                      "Tuft cells",
                      "Goblet cells",
                      "Enterocytes (immature)" ) 
timecourse@active.ident <- plyr::mapvalues(x = timecourse@active.ident, from = current.cluster.ids, to = new.cluster.ids)
timecourse$leiden <- plyr::mapvalues(x = timecourse$leiden, from = current.cluster.ids, to = new.cluster.ids)

DimPlot(timecourse, reduction.used = "umap", 
        pt.size = .9, 
        label = TRUE, label.size = 5, repel =TRUE)

#nice ggplot
UMAP_data <- FetchData(timecourse, vars = c("UMAP_1", "UMAP_2")) 
UMAP_data$leiden <- as.factor(timecourse$leiden)
ggplot(data=UMAP_data, aes(UMAP_1, UMAP_2, fill=factor(leiden))) + 
  theme_classic()+ 
  geom_point(shape = 21,color="black", size =2.5)+ 
  guides(fill=guide_legend(title="CLUSTER"), colour=guide_legend(override.aes = list(size = 4)))+
  scale_fill_manual(values = c(wes_palette("Darjeeling1"), wes_palette("Cavalcanti1")))+ theme(legend.justification = "right") + theme(text = element_text(size=30))

#############Cluster breakdown per sample###########
numberofcells         <- table(timecourse$orig.ident, timecourse$leiden)
totalcellspersample   <- c(sum(numberofcells[1,]), sum(numberofcells[2,]), sum(numberofcells[3,]), sum(numberofcells[4,]), sum(numberofcells[5,]), sum(numberofcells[6,]))
a                     <- cbind(numberofcells,totalcellspersample)
totalcellspercluster  <- c(sum(a[,1]), sum(a[,2]), sum(a[,3]), sum(a[,4]), sum(a[,5]), sum(a[,6]), 
                          sum(a[,7]), sum(a[,8]),  sum(a[,9]), sum(a[,10]), sum(a[,11]))
b                     <- rbind(a, totalcellspercluster)

c0 <- (b[1:6,1]/totalcellspersample)*100
c1 <- (b[1:6,2]/totalcellspersample)*100
c2 <- (b[1:6,3]/totalcellspersample)*100
c3 <- (b[1:6,4]/totalcellspersample)*100
c4 <- (b[1:6,5]/totalcellspersample)*100
c5 <- (b[1:6,6]/totalcellspersample)*100
c6 <- (b[1:6,7]/totalcellspersample)*100
c7 <- (b[1:6,8]/totalcellspersample)*100
c8 <- (b[1:6,9]/totalcellspersample)*100
c9 <- (b[1:6,10]/totalcellspersample)*100

c <- rbind(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
colSums(c)
c
rownames(c) <-  c( "Early progenitors",
                   "Goblet cells",
                   "Enterocytes 1",
                   "Paneth cells",
                   "Stem cells",
                   "Enterocytes 2",
                   "Enterocytes 3",
                   "Enterocytes (immature)",
                   "T cells",
                   "Tuft cells")
c <-c[,c(5,4,6,2,1,3)]
colnames(c) <- c("D164A 4d pi", "D164A 2d pi", "D164A 0d pi", "control 4d pi", "control 2d pi", "control 0d pi")

par(mar=c(5,10,3,20))
par(xpd=TRUE)
barplot(c, horiz=TRUE,
        legend = TRUE,
        args.legend=list(bty = "n", x=225, cex=1.5),
        main = "Cluster breakdown per sample", 
        las = 1, 
        col=c(wes_palette("Darjeeling1"), wes_palette("Darjeeling2")),
        cex.axis=1.5, cex.names=1.5, cex.main=2)

#################SUBSETTING OF STEM CELLS AND EP TOGETHER###############
Idents(timecourse)  <- "leiden" #set leiden clustering as active identity 
crypt               <- subset(timecourse, idents = c("Stem cells", "Early progenitors"))
length(crypt@active.ident)
DimPlot(crypt, reduction.used = "umap", 
        pt.size = .9, 
        label = TRUE, label.size = 5, repel =TRUE, ncol=3)
FeaturePlot(crypt, features = "Lyz1")

crypt <- FindVariableFeatures(object = crypt, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
crypt <- ScaleData(object = crypt, features = rownames(x =crypt), vars.to.regress = c("nCount_RNA", "percent.mito", "orig.ident")) 
crypt <- RunPCA(object = crypt, features = VariableFeatures(object = crypt), npcs = 20, verbose = FALSE)
ElbowPlot(object = crypt)
DimPlot(crypt, reduction = "pca", group.by = "orig.ident")
crypt <- RunUMAP(object = crypt, dims = 1:10, set.seed = 5)
crypt <- FindNeighbors(crypt, dims = 1:10, set.seed = 5) 
DimPlot(crypt, reduction.used = "umap", label = FALSE, pt.size = .3, group.by = "orig.ident")

#################NORMALIZATION###############
Idents(crypt) <- "orig.ident"

#pull data matrix from samples (normalized, non-scaled data) and apply exp to remove the log
D164A_0d              <- subset(crypt, idents = "D164A_ni")
D164A_0d_norm.data    <- exp(D164A_0d[["RNA"]]@data)

D164A_2d              <- subset(crypt, idents = "D164A_2d")
D164A_2d_norm.data    <- exp(D164A_2d@assays[["RNA"]]@data)

D164A_4d              <- subset(crypt, idents = "D164A_4d")
D164A_4d_norm.data    <- exp(D164A_4d@assays[["RNA"]]@data)

control_0d            <- subset(crypt, idents = "control_ni")
control_0d_norm.data  <- exp(control_0d[["RNA"]]@data)

control_2d            <- subset(crypt, idents = "control_2d")
control_2d_norm.data  <- exp(control_2d@assays[["RNA"]]@data)

control_4d            <- subset(crypt, idents = "control_4d")
control_4d_norm.data  <- exp(control_4d@assays[["RNA"]]@data)

#compute means of controls
mean_control_0d         <- rowMeans(control_0d_norm.data, na.rm = TRUE)
names(mean_control_0d)  <- rownames(control_0d_norm.data)

mean_control_2d         <- rowMeans(control_2d_norm.data, na.rm = TRUE)
names(mean_control_2d)  <- rownames(control_2d_norm.data)

mean_control_4d         <- rowMeans(control_4d_norm.data, na.rm = TRUE)
names(mean_control_4d)  <- rownames(control_4d_norm.data)

#divide by D164A by mean of respective control, then apply log again (Seurat wants data in log scale)
normD164A_0d <- log(D164A_0d_norm.data/mean_control_0d)
normD164A_0d <- as(as.matrix(normD164A_0d), 'sparseMatrix')

normD164A_2d <- log(D164A_2d_norm.data/mean_control_2d)
normD164A_2d <- as(as.matrix(normD164A_2d), 'sparseMatrix')

normD164A_4d <- log(D164A_4d_norm.data/mean_control_4d)
normD164A_4d <- as(as.matrix(normD164A_4d), 'sparseMatrix')

#create Seurat object and input the data in the data slot (where data are after normalization), then scale
normD164A_0d_Seurat                       <- CreateSeuratObject(normD164A_0d)
normD164A_0d_Seurat@assays[["RNA"]]@data  <- normD164A_0d
normD164A_0d_Seurat                       <- ScaleData(normD164A_0d_Seurat)

normD164A_2d_Seurat                       <- CreateSeuratObject(normD164A_2d)
normD164A_2d_Seurat@assays[["RNA"]]@data  <- normD164A_2d
normD164A_2d_Seurat                       <- ScaleData(normD164A_2d_Seurat)

normD164A_4d_Seurat                       <- CreateSeuratObject(normD164A_4d)
normD164A_4d_Seurat@assays[["RNA"]]@data  <- normD164A_4d
normD164A_4d_Seurat                       <- ScaleData(normD164A_4d_Seurat)

#merge into normalized crypt Seurat object for downstream analysis
norm.crypt <- merge(normD164A_0d_Seurat, c(normD164A_2d_Seurat, normD164A_4d_Seurat))
length(norm.crypt@active.ident)
Idents(norm.crypt) <- "orig.ident"
RenameIdents(norm.crypt@active.ident, old.ident = "D164A ni", new.ident = "D164A 0d pi")
old.identities <- c("D164A ni", "D164A 2d pi", "D164A 4d pi")
new.identities <- c("D164A 0d pi", "D164A 2d pi", "D164A 4d pi")
norm.crypt@active.ident <- plyr::mapvalues(x = norm.crypt@active.ident, from = old.identities, to = new.identities)

norm.crypt <- FindVariableFeatures(norm.crypt)
norm.crypt <- ScaleData(norm.crypt)
norm.crypt <- RunPCA(norm.crypt)
ElbowPlot(object = norm.crypt) #Find significant PCAs; 15 are significant
norm.crypt <- RunUMAP(object = norm.crypt, dims = 1:10, set.seed = 5)
norm.crypt <- FindNeighbors(norm.crypt, dims = 1:15, set.seed = 5) 
DimPlot(norm.crypt, reduction.used = "umap", label = FALSE, pt.size =1,
        cols= wes_palette("Darjeeling1"))
FeaturePlot(norm.crypt, features="cc.score", pt.size=.6, sort.cell =T)+ scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))


######DIFFUSION MAP WITH DESTINY FOR DIFFERENTIATION TRAJECTORY########
# Get mean expression of genes of interest per cell: stem cell score
features <- c("Olfm4", "Lgr5", "Smoc2", "Soat1","Slc12a2","Ascl2", "Axin2", "Gkn3")
mean.exp <- colMeans(x = norm.crypt@assays$RNA[features, ], na.rm = TRUE)

# Add mean expression values in 'object@meta.data$gene.set.score'
if (all(names(x = mean.exp) == rownames(x = norm.crypt@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  norm.crypt@meta.data$SC.score <- mean.exp
}

# Get mean expression of genes of interest per cell: differentiation score
features <- c("Lyz1", "Muc2", "Defa24", "Defa17", "Alpi", "Krt20", "Ace2", "Apoa1", "Muc3")
mean.exp <- colMeans(x = norm.crypt@assays$RNA[features, ], na.rm = TRUE)

# Add mean expression values in 'object@meta.data$gene.set.score'
if (all(names(x = mean.exp) == rownames(x = norm.crypt@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  norm.crypt@meta.data$diff.score <- mean.exp
}

#convert to singlecellexperiment
D164A_sce <- as.SingleCellExperiment(norm.crypt)
dm        <- DiffusionMap(D164A_sce)

plot(dm2,1:2,
     pch = 16, # pch for prettier points
     col_by = "Olfm4",
     legend_main = 'norm TA and SC clusters')

#nice plot
ggplot(dm2) + theme_classic()+
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(data=dm, aes(DC1, DC2, fill=factor(orig.ident)), shape = 21,color="black", size =3) +
  scale_fill_manual(values = wes_palette("Darjeeling1"), name=" ") +
  theme(text = element_text(size=30)) +labs (x = "DC 1", y = "DC 2")

#violin plot
ggplot(dm2, aes(x=orig.ident, y=DC1, fill=factor(orig.ident))) + 
  geom_violin(color="black", width = 2)+ theme_classic()+ theme(legend.position = "none") + 
  scale_fill_manual(values = wes_palette("Darjeeling1"), name=" ") + 
  theme(text = element_text(size=30)) +
  stat_summary(fun.y=mean, geom="point", color="black") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),   axis.ticks.x=element_blank())

#write table to calculate Wilcoxon test (alternative to Student T test when you can't assume differences to be normally distributed)
View(dm2$DC1)
stat <- as.data.frame(dm2$DC1)
stat$time <- rownames(stat)
time0 <- stat[grep("D164A_ni", stat$time), ]
nrow(time0)
time2 <- stat[grep("D164A_2d", stat$time), ]
nrow(time2)
time4 <- stat[grep("D164A_4d", stat$time), ]
nrow(time4)

wilcox.test(time0$`dm2$DC1`, time2$`dm2$DC1`, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(time0$`dm2$DC1`, time4$`dm2$DC1`, alternative = "two.sided") #p-value < 2.2e-16

#plot showing expression of signatures
ggplot(dm) + theme_classic()+
  geom_point(data=dm, aes(DC1, DC2, fill=SC.score), shape = 21,color="gray47", size =2) +
  scale_fill_gradientn(colours= rev(brewer.pal(n = 11, name = "RdYlBu")), name=" ") +
  theme(text = element_text(size=15)) + labs (x = "DC 1", y = "DC 2")+
  labs(title = "Olfm4, Lgr5, Smoc2, Soat1, Slc12a2, \nAscl2, Axin2, Gkn3, Wdr12, Tnfrsf19") + 
  theme(plot.title = element_text(size=15, face="italic"))

ggplot(dm) + theme_classic()+
  geom_point(data=dm, aes(DC1, DC2, fill=diff.score), shape = 21, color="gray47", size =2) +
  scale_fill_gradientn(colours= rev(brewer.pal(n = 11, name = "RdYlBu")), name=" ") +
  theme(text = element_text(size=15)) +labs (x = "DC 1", y = "DC 2")+
  labs(title = "Lyz1, Muc2, Defa24, Defa17, Alpi, \nKrt20, Ace2, Apoa1, Muc3") + 
  theme(plot.title = element_text(size=15, face="italic"))

#################CELL CYCLE###############
cc.genes <- readLines(con = "/home/moorlab/mnt_rstudio/Coco/regev_lab_cell_cycle_genes.txt")
cc.genes <- str_to_title(cc.genes) #makes capital to small case
#use all cc.genes (no G2m or S phase difference)
tmp     <- CellCycleScoring(object = norm.crypt, s.features = cc.genes, g2m.features  = cc.genes, set.ident = TRUE)
head(tmp[[]], 50)

#ridgeplot
RidgePlot(tmp, features=(vars = c("G2M.Score")), group.by = "old.ident", wes_palette("Darjeeling1")) +
  theme_classic() +  
  geom_vline(xintercept = 0, linetype="dashed") +
  theme(text = element_text(size=25)) + labs (x = "Cell cycle score", y = "")

statCC      <- as.data.frame(tmp$G2M.Score)
statCC$time <- rownames(statCC)
time0CC     <- statCC[grep("D164A_ni", statCC$time), ]
nrow(time0CC)
time2CC     <- statCC[grep("D164A_2d", statCC$time), ]
time4CC     <- statCC[grep("D164A_4d", statCC$time), ]

wilcox.test(time0CC$`tmp$G2M.Score`, time2CC$`tmp$G2M.Score`, alternative = "two.sided") #p-value  5.92e-11
wilcox.test(time0CC$`tmp$G2M.Score`, time4CC$`tmp$G2M.Score`, alternative = "two.sided") #p-value  < 2.2e-16

#add cc.score to metadata
norm.crypt@meta.data$cc.score <- tmp$G2M.Score
FeaturePlot(object = norm.crypt, features = "cc.score") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

########PSUPERTIME##############
data_to_write_out  <- as.data.frame(as.matrix(norm.crypt@assays$RNA@scale.data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "/home/moorlab/mnt_rstudio/Coco/BetacatProject/scRNASeq/norm.crypt.txt")

# make sce
dt      = fread("/home/moorlab/mnt_rstudio/Coco/BetacatProject/scRNASeq/norm.crypt.txt")
sce     = SingleCellExperiment(assays=list(logcounts=as.matrix(dt[,-1, with=FALSE])))
rownames(sce)   = dt$V1
colnames(sce)   = names(dt)[-1]

# extract info for cells
col_dt          = data.table(raw_name = colnames(sce))
col_dt[, condition := factor(str_match(raw_name, '(\\w+) (\\w+) (\\w+) (\\w+)')[, 2])]
col_dt[, day := factor(str_match(raw_name, '(\\w+) (\\w+) (\\w+) (\\w+)')[, 3])]
col_dt[, day_cond := factor(str_match(raw_name, '(\\w+) (\\w+) (\\w+) (\\w+)')[, 5])]
col_dt[, cell_id := str_match(raw_name, '(\\w+) (\\w+) (\\w+) (\\w+)')[, 4]]

# check levels ok, add to sce
for (n in names(col_dt))
  print(levels(col_dt[[n]]))
colData(sce)    = as(col_dt, 'DataFrame')

# run psupertime
psuper_obj      = psupertime(sce, sce$day, sel_genes='all', scale=FALSE, min_expression=0)

# do plot
plot_identified_genes_over_psupertime(psuper_obj, label_name='time post injection', palette="YlOrRd")
plot_identified_gene_coefficients(psuper_obj)
plot_specified_genes_over_psupertime(psuper_obj, c("Olfm4", "Cd74", "Soat1","Smoc2", "Arglu1", "Slc12a2", "Mki67", "Ccnd2"), label_name='time post injection', palette="YlOrRd")

#cycling genes
plot_specified_genes_over_psupertime(psuper_obj, c("Olfm4", "Smoc2","Slc12a2", "Mki67", "Cenpe",  "Cdc20"), label_name='time post injection', palette="YlOrRd")

extra_genes <- c("Olfm4", "Smoc2","Slc12a2", "Mki67", "Cenpe",  "Cdc20")
proj_dt     = psuper_obj$proj_dt
beta_dt     = psuper_obj$beta_dt
x_data      = psuper_obj$x_data
params      = psuper_obj$params

# restrict to specified genes
extra_genes = intersect(extra_genes, colnames(x_data))
if (length(extra_genes)==0) {
  warning('genes not found; did not plot')
  return()
}

# set up data
plot_wide   = cbind(proj_dt, data.table(x_data[, extra_genes, drop=FALSE]))
plot_dt     = melt.data.table(
  plot_wide, 
  id 				= c("psuper", "label_input", "label_psuper"), 
  measure 		= extra_genes, 
  variable.name 	= "symbol"
)
plot_dt[, `:=`(symbol, factor(symbol, levels = extra_genes))]

# set up plot
n_genes 	= length(extra_genes)
ncol 		= ceiling(sqrt(n_genes*plot_ratio))
nrow 		= ceiling(n_genes/ncol)
plot_ratio=1

# plot
ggplot(plot_dt) + 
  aes( x=psuper, y=value ) +
  geom_point( size=3, aes(fill=factor(label_input) ), shape = 21,color="black", size =3)+
  geom_smooth(se=T, colour='black') + geom_hline(yintercept = 0, linetype="dashed")+
  scale_fill_manual( values=wes_palette("Darjeeling1"), name="days pi" ) +
  #scale_x_continuous( breaks=pretty_breaks() ) +
  #scale_y_continuous( breaks=pretty_breaks() ) +
  facet_wrap( ~ symbol, scales='fixed', nrow=nrow, ncol= 3 ) +
  ylim(-1,3)+
  theme_bw() + theme(text = element_text(size=20)) +
  theme(axis.text.x = element_blank()) +
  labs(x 		= 'psupertime',y 		= 'z-scored log2 expression')+
  labs(title = "Stem cell markers") + 
  theme(plot.title = element_text(size=20, face="bold"))+ theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

# run psupertime with mouse transcription factors
psuper_obj_TF      = psupertime(sce, sce$day, sel_genes='tf_mouse', scale=FALSE, min_expression=0)
plot_identified_genes_over_psupertime(psuper_obj_TF, label_name='time post injection', palette="YlOrRd")
plot_identified_gene_coefficients(psuper_obj_TF) +
  theme(text = element_text(size=15))+ labs(x = "Transcription factor")+
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) 

############DIFFERENTIAL EXPRESSION AND GSEA############
#set original identity as active identity 
Idents(timecourse) <- "orig.ident"

#markers and heatmap
markers         <- FindAllMarkers(norm.crypt)
View(markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC))
top100          <- markers %>% top_n(200, avg_logFC)
genes.to.plot   <- top100$gene
genes.to.label  <- c("Olfm4", "Slc12a2", "Mki67", "Ccnd1","Car2", "Ada", "Apoa4", "Apoa1", "Krt20", "Alpi", "Ace2" , "Muc3",  "Lyz1")
labels          <- rep(x = "transparent", times = length(x = genes.to.plot))
labels[c(31  , 43  ,50,139, 143, 151 ,163 ,181, 193 ,199) ] <- "black"

DoHeatmap(norm.crypt, slot="data", features = genes.to.plot, group.colors	= wes_palette("Darjeeling1"), angle=0 , draw.lines=T)+ 
  scale_fill_gradientn(colors = brewer.pal(11,"RdYlBu")[11:1])  +
  theme(axis.text.y = element_text(face = "bold", color = rev(x = labels))) 

geom_text_repel(aes(label=ifelse(timepoint>2,as.character(pathway),'')), direction="y", nudge_x=7, size = 5 , segment.size = 0.1, xjust  = 3) +
  
#GSEA
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

#on D164A 0d vs 2d
D164A_0vs2d <- FindMarkers(norm.crypt, ident.1 = "D164A 2d pi", ident.2 = "D164A 0d pi", verbose = FALSE)
ranks <- D164A_0vs2d %>% 
  na.omit()%>%
  mutate(ranking=-log10(p_val_adj)/sign(avg_logFC))
ranks <- ranks$ranking
names(ranks) <- rownames(D164A_0vs2d)
head(ranks, 10)

BP_0vs2 <- fgsea(pathways = BP, 
                  stats = ranks,
                  minSize=10,
                  maxSize=500,
                  nperm=1000000)

Hallmarks_0vs2
CGP_0vs2 %>% filter(abs(NES)>1.5 & pval<0.05)
KEGG_0vs2 %>% filter(abs(NES)>1.5 & pval<0.05)
REACTOME_0vs2 %>% filter(abs(NES)>1.5 & pval<0.05)
TFT_0vs2 %>% filter(abs(NES)>1.5 & pval<0.05) #no sign padj
BP_0vs2  %>% filter(abs(NES)>1.5 & pval<0.05)

#on D164A 0d vs 4d
D164A_0vs4d <- FindMarkers(norm.crypt, ident.1 = "D164A 4d pi", ident.2 = "D164A 0d pi", verbose = FALSE)
ranks <- D164A_0vs4d %>% 
  na.omit()%>%
  mutate(ranking=-log10(p_val_adj)/sign(avg_logFC))
ranks <- ranks$ranking
names(ranks) <- rownames(D164A_0vs4d)
head(ranks, 10)

BP_0vs4 <- fgsea(pathways = BP, 
                       stats = ranks,
                       minSize=10,
                       maxSize=500,
                       nperm=1000000)

Hallmarks_0vs4  %>% filter(abs(NES)>1.5 & padj<0.05)
CGP_0vs4        %>% filter(abs(NES)>1.5& padj<0.05)
KEGG_0vs4       %>% filter(abs(NES)>1.5 & padj<0.05)
REACTOME_0vs4   %>% filter(abs(NES)>1.5 & padj<0.05)
TFT_0vs4        %>% filter(abs(NES)>1.5 & pval<0.05)
BP_0vs4         %>% filter(abs(NES)>1.5 & padj<0.05)

#get data for plotting
Hallmarks_024   <- merge(Hallmarks_0vs2[,c("pathway","NES", "pval", "size")], Hallmarks_0vs4[,c("pathway","NES", "pval", "size")], by= "pathway")
KEGG_024        <- merge(KEGG_0vs2[,c("pathway","NES", "pval", "size")], KEGG_0vs4[,c("pathway","NES", "pval", "size")], by= "pathway")
REACTOME_024    <- merge(REACTOME_0vs2[,c("pathway","NES", "pval", "size")], REACTOME_0vs4[,c("pathway","NES", "pval", "size")], by= "pathway")
CGP_024         <- merge(CGP_0vs2[,c("pathway","NES", "pval", "size")], CGP_0vs4[,c("pathway","NES", "pval", "size")], by= "pathway")
TFT_024         <- merge(TFT_0vs2[,c("pathway","NES", "pval", "size")], TFT_0vs4[,c("pathway","NES", "pval", "size")], by= "pathway")
BP_024          <- merge(BP_0vs2[,c("pathway","NES", "pval", "size")], BP_0vs4[,c("pathway","NES", "pval", "size")], by= "pathway")

#prepare for plotting
enrichment        <- rbind(Hallmarks_024, KEGG_024, REACTOME_024,CGP_024, TFT_024, BP_024)
gsea              <- rbind(enrichment[,c("pathway", "NES.x", "pval.x", "size.x")], enrichment[,c("pathway", "NES.y", "pval.y", "size.y")], fill=TRUE)
gsea[2987:5972,2] <- gsea[2987:5972,5]
gsea[2987:5972,3] <- gsea[2987:5972,6]
gsea[2987:5972,4] <- gsea[2987:5972,7]
gsea[,c(5:7)]     <- NULL
time0             <- gsea[1:2986, c("pathway", "NES.x", "pval.x", "size.x")]
time0$NES.x       <- c(0)
time0$pval.x      <- c(0)
time0$size.x      <- c(1)
gseaplot          <- rbind(time0, gsea)
gseaplot$timepoint <- c(0)
gseaplot$timepoint[2987:5972] <- c(2)
gseaplot$timepoint[5973:8958] <- c(4)
colnames(gseaplot) <- c("pathway", "NES", "pval", "size", "timepoint")

ggplot(data = gseaplot, aes(x = timepoint, y = NES, group = pathway, label= pathway))+ 
  theme_classic()+ expand_limits(x = 6)+
  geom_smooth(aes(color=..x..), size=.25) +
  theme(axis.text.y = element_text(face="bold", size=15, colour = "black"))+
  theme(axis.text.x = element_text(face="bold", size=15, colour = "black"))+
  #theme(legend.position=c(0, 3), legend.direction = "vertical")+
  guides(col= FALSE)+
  scale_colour_gradientn(colours = colorRampPalette(wes_palette("Darjeeling1", n=3))(length(seq(0, 4.5, by = 0.04))))+
  #add darkred lines for selected pathways
  geom_smooth(data=gseaplot[c(71, 3057, 6043),], col="black", size=.5) + #GO_RIBOSOME
  annotate("label", x= 2.1, y=gseaplot$NES[3057], label = "ribosome", size=4.2, hjust = 0)+
  geom_smooth(data=gseaplot[c(160, 3146, 6132),], col="black", size=.5) + #REACTOME_EUKARYOTIC_TRANSLATION_INITIATION
  annotate("label", x= 2.1, y=gseaplot$NES[3146],label = "translation", size=4.2, hjust = 0)+
  geom_smooth(data=gseaplot[c(9, 2995, 5981),], col="black", size=.5) + #HALLMARK_E2F_TARGETS
  annotate("label", x= 4.1, y=gseaplot$NES[5981]-0.2,label = "E2F targets",col="black", size=4.2, hjust = 0)+
  geom_smooth(data=gseaplot[c(13, 2999, 5985),], col="black", size=.5) + #HALLMARK_G2M CHECKPOINT
  annotate("label", x= 4.1, y=gseaplot$NES[5985],label = "G2M checkpoint",col="black",  size=4.2, hjust = 0)+
  geom_smooth(data=gseaplot[c(43, 3029, 6015),], col="black", size=.5) + #REACTOME_CELL_CYCLE
  annotate("label", x= 4.1, y=gseaplot$NES[6015],label = "cell cycle",col="black",  size=4.2, hjust = 0)+
  geom_smooth(data=gseaplot[c(21, 3007, 5993),], col="black", size=.5) + #KRAS SIGNALLING
  annotate("label", x= 4.1, y=gseaplot$NES[5993],label = "Kras signalling", size=4.2, hjust = 0)+
  geom_smooth(data=gseaplot[c(2020, 5006, 7992),], col="black", size=.5) + #DIGESTION
  annotate("label", x= 4.1, y=gseaplot$NES[5992],label = "digestion", size=4.2, hjust = 0)+
  geom_smooth(data=gseaplot[c(1108, 4094, 7080),], col="black", size=.5) + 
  annotate("label", x= 4.1, y=gseaplot$NES[7080],label = "APC targets", col="black", size=4.2, hjust = 0)+
  geom_smooth(data=gseaplot[c(31, 3017, 6003),], col="black", size=.5) + 
  annotate("label", x= 4.1, y=gseaplot$NES[6003]-0.1,label = "TNFA signalling via NFKB", col="black", size=4.2, hjust = 0)+
  geom_smooth(data=gseaplot[c(2514, 5500, 8486),], col="black", size=.5) + 
  annotate("label", x= 4.1, y=gseaplot$NES[8486],label = "stress signalling", col="black", size=4.2, hjust = 0)+
  geom_smooth(data=gseaplot[c(609, 3595, 6581),], col="black", size=.5) + 
  annotate("label", x= 4.1, y=gseaplot$NES[6581],label = "downregulated upon Ctnnb1 KO", col="black", size=4.2, hjust = 0)+
  #geom_smooth(data=gseaplot[c(1574, 4560, 7546),], col="black", size=.5) + 
  #annotate("label", x= 4.1, y=gseaplot$NES[7546],label = "Lef1 targets", size=4.2, hjust = 0)+
  #geom_smooth(data=gseaplot[c(1546, 4532, 7518),], col="black", size=.5) + #LEF1 TARGETS
  #annotate("label", x= 4.1, y=gseaplot$NES[7518],label = "HNF1", size=4.2, hjust = 0)+
  geom_smooth(data=gseaplot[c(20, 3006, 5992),], col="black", size=.5) + #endogenous antigen presentation
  annotate("label", x= 4.1+0.7, y=gseaplot$NES[5992],label = "interferon g signalling",col="black", size=4.2, hjust = 0)+
  geom_smooth(data=gseaplot[c(1838, 4824, 7810),], col="black", size=.5) + #endogenous antigen presentation
  annotate("label", x= 4.1, y=gseaplot$NES[7810],label = "endog. antig. present.",col="black", size=4.2, hjust = 0)+
  geom_point(data=gseaplot, aes(timepoint, NES, fill=pval, size=size), shape = 21,color="black")+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_fill_gradientn(colours = colorRampPalette(c("gray0", "gray50", "gray100"))(length(seq(0, 4, by = 0.04))))+
  #scale_fill_manual( values=wes_palette("Darjeeling1"), name="days pi" ) 
  #geom_text(aes(label=ifelse(timepoint>2,as.character(pathway),'')),hjust=0,vjust=0) 
  #geom_text_repel(aes(label=ifelse(timepoint>2,as.character(pathway),'')), direction="y", nudge_x=7, size = 5 , segment.size = 0.1, xjust  = 3) +
  #geom_text_repel(aes(label=ifelse(timepoint==2&NES<0, as.character(pathway),'')), size = 5 ) +
  #geom_text_repel(aes(label=ifelse(timepoint>2&NES<0, as.character(pathway),'')), nudge_x=7, direction="y") +
  theme(text = element_text(size=15)) +
  labs(x 		= 'days post injection',y	= 'NORMALIZED ENRICHMENT SCORE')


####CLASSIFIER by A. Lafzi######
logistic.reg <- function(scale.data.ref, clus.ref, top.100.ref.markers, scale.data.test){
  
  train.data <- scale.data.ref
  test.data <- scale.data.test
  train.lables <- clus.ref[colnames(train.data)]
  
  var.train <- sapply(top.100.ref.markers, function(x) Matrix::colSums(train.data[x,]))
  p <- lapply(top.100.ref.markers, function(x) rownames(test.data)[which(rownames(test.data) %in% x)])
  var.test <- sapply(p, function(x) Matrix::colSums(test.data[x, ]))
  
  model.train <- data.frame(train.lables, var.train)
  model.test <- data.frame(var.test)
  
  mod <- glm(train.lables ~., data = model.train, family = binomial)
  fitted.results <- predict(mod, newdata = model.test, "response")
  
  return(fitted.results)
  
}

#training dataset: single cell survey of intestinal epithelium Haber et al 2017
atlas_UMIcounts     <- read.table(file = '~/mnt_rstudio/Coco/scRNASeq/RegevAtlas/GSE92332_atlas_UMIcounts.txt', as.is = TRUE)

#get cluster information from cell id 
metaData            <- data.frame(cellNames = colnames(atlas_UMIcounts))%>%
  mutate(cluster = factor(str_replace(cellNames,".+_","")))
rownames(metaData)  <- metaData$cellNames
print(metaData)

#standard Seurat processing
atlas1    <- CreateSeuratObject(counts = atlas_UMIcounts, meta.data = metaData)
colnames(x = atlas1[[]])
Idents(object = atlas1) <- 'cluster'
levels(x = atlas1)
length(atlas1$orig.ident)
FeatureScatter(atlas1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(atlas1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
atlas1    <- NormalizeData(atlas1, normalization.method = "LogNormalize", scale.factor = 10000)
atlas1    <- FindVariableFeatures(atlas1,  selection.method = "vst", nfeatures = 2000) 
VariableFeaturePlot(atlas1)
all.genes <- rownames(atlas1)
atlas1    <- ScaleData(atlas1, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mito", "orig.ident"))
atlas1    <- RunPCA(atlas1, verbose = FALSE, set.seed(1), features = VariableFeatures(object = atlas1), npcs = 20)
DimHeatmap(atlas1, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(atlas1)
atlas1    <- FindNeighbors(atlas1, dims = 1:20)
atlas1    <- FindClusters(atlas1, resolution = 0.5)
atlas1    <- RunUMAP(object = atlas1, dims = 1:10, set.seed = 5)
atlas1    <- FindNeighbors(atlas1, dims = 1:10, set.seed = 5) 
DimPlot(atlas1, reduction.used = "umap", label = FALSE, pt.size = .3, group.by = "cluster")

#subset SC and TAs
Idents(atlas1)  <- "cluster"
cryptREGEV      <- subset(atlas1, idents = c("Stem", "TA.G2", "TA.G1", "TA.Early"))
DimPlot(cryptREGEV)
                     
#make only one TA cluster                     
current.cluster.ids     <- c("Stem", "TA.G2",  "TA.G1", "TA.Early")
new.cluster.ids         <- c("Stem", "TA", "TA", "TA")
cryptREGEV@active.ident <- plyr::mapvalues(x = cryptREGEV@active.ident, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(cryptREGEV)
                     
#find markers of SC vs TAs
markers_regev <- FindAllMarkers(object = cryptREGEV, only.pos = T, min.pct = 0.25, logfc.threshold = 0.2)
gene_cl.ref   <- cut_markers(levels(markers_regev$cluster), markers_regev, ntop=100)
gene_cl.ref
                     
#test dataset                     
norm0 <- subset(norm.crypt, ident="D164A 0d pi")
norm2 <- subset(norm.crypt, ident="D164A 2d pi")

#input for logistic regression
scale.data.ref      <- cryptREGEV@assays$RNA@scale.data
clus.ref            <- cryptREGEV@active.ident
top.100.ref.markers <- gene_cl.ref

#classifier                     
TA.probs0         <- logistic.reg(scale.data.ref, clus.ref, top.100.ref.markers, norm0@assays$RNA@scale.data)
TA.probs2         <- logistic.reg(scale.data.ref, clus.ref, top.100.ref.markers, norm2@assays$RNA@scale.data)
sorted.TA.probs   <- sort(TA.probs)
TA.probs.df       <- data.frame(prob= sorted.TA.probs, cells = names(TA.probs))
sorted.TA.probs2  <- sort(TA.probs2)
TA.probs.df2      <- data.frame(prob= sorted.TA.probs2, cells = names(TA.probs2))

#plot
ggplot()+
  geom_density(data= TA.probs.df0, aes(x=prob, y= -..density..), col="black", size=.8, fill=wes_palette("Darjeeling1")[1])+
  geom_density(data= TA.probs.df2, aes(x=prob),  col="black", size=.8, fill=wes_palette("Darjeeling1")[2])+
  geom_vline(xintercept = 0.5, linetype="dashed") +
  annotate("text", x= 0.45, y=5,label = "IESC-like", size=6, hjust = 0)+
  annotate("text", x= 0.55, y=5,label = "TA-like", size=6, hjust = 0)+
  theme(text=element_text(size=17,  family="arial"))+ coord_flip()+
  labs(x 		= 'probability of being classified as TA vs IESC')+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

ks.test(TA.probs.df0$prob, TA.probs.df2$prob, alternative = c("two.sided"))
