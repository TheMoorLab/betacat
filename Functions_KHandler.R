CreateSeurat_200_3_mouse <- function(annotated_matrix, projectName, minCells, minFeatures, Variable_mean_cutoff, Variable_dispersion_cutoff, pc) {
  #this function creates a Seurat object with the initial cutoffs of min.cell=3 and min.features=200
  #then it also defines the mitochondrial gene content, adds it meta data
  #returns initial violin plots and median values of number of genes,...
  #makes the initial statistical tests for PCA filtering 
  #define file_name.rdata to save the Seruat object without any further filtering or cutoffs 
  #Variable_mean_cutoff is c(x, y; from tutorial c(0.0125, 3)
  #Dispersion cutoff is c(0.5, Inf) 
  tmp <- CreateSeuratObject(counts = annotated_matrix, project = projectName, min.cells = minCells, min.features = minFeatures)
  length(tmp$orig.ident)
  mito.features <- grep(pattern = "^mt-", x = rownames(x = tmp), value = TRUE)
  print(mito.features)
  percent.mito <- Matrix::colSums(x = GetAssayData(object = tmp, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = tmp, slot = 'counts'))
  tmp <- AddMetaData(object = tmp, metadata = percent.mito, col.name = "percent.mito")
  print (head(tmp@meta.data))
  print (VlnPlot(object = tmp, features = c("nCount_RNA", "nFeature_RNA", "percent.mito")))
  print (median(tmp@meta.data$nCount_RNA))
  print (median(tmp@meta.data$nFeature_RNA))
  print(median(tmp@meta.data$percent.mito))
  print(FeatureScatter(object = tmp, feature1 = "nCount_RNA", feature2 = "percent.mito"))
  print(FeatureScatter(object = tmp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  tmp <- NormalizeData(object = tmp, normalization.method = "LogNormalize", scale.factor = 1e4)
  tmp <- FindVariableFeatures(object = tmp, selection.method = 'mean.var.plot', mean.cutoff = Variable_mean_cutoff, dispersion.cutoff = Variable_dispersion_cutoff)
  print(length(x = VariableFeatures(object = tmp)))
  tmp <- ScaleData(object = tmp, features = rownames(x = tmp), vars.to.regress = c("nCount_RNA", "percent.mito"))
  tmp <- RunPCA(object = tmp, features = VariableFeatures(object = tmp), npcs = pc, verbose = FALSE)
  print(ElbowPlot(object = tmp))
  print(length(tmp$orig.ident))
  return(tmp)
}

Annotation_mouse <- function(zUMI_output) {
  #This function returns the annotated count matrix with gene names 
  ensembl<-useEnsembl(biomart="ensembl",GRCh=37)
  list<-listDatasets(ensembl)
  mart <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version=84)
  attributes<-listAttributes(mart)
  gene_ids<-getBM(attributes = c("ensembl_gene_id","external_gene_name"), mart = mart)
  tmp<-as.matrix(zUMI_output$umicount$exon$all)
  tmp<-mutate(as.data.frame(tmp),ensembl_gene_id=rownames(tmp))
  tmp<-full_join(tmp,gene_ids)
  print(length(unique(tmp$external_gene_name)))
  tmp<-tmp[!duplicated(tmp$external_gene_name),]
  tmp[is.na(tmp)]<-0 #make all empty value to zero
  rownames(tmp)<-tmp$external_gene_name
  tmp<-dplyr::select(tmp,-external_gene_name,-ensembl_gene_id)
  return(tmp)
}

Annotation_mouse_reads <- function(zUMI_output, prefix) {
  #This function returns the annotated count matrix with gene names 
  ensembl<-useEnsembl(biomart="ensembl",GRCh=37)
  list<-listDatasets(ensembl)
  mart <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version=84)
  attributes<-listAttributes(mart)
  gene_ids<-getBM(attributes = c("ensembl_gene_id","external_gene_name"), mart = mart)
  tmp<-as.matrix(zUMI_output$readcount$exon$all)
  colnames(tmp)<-paste(colnames(tmp), prefix) 
  tmp<-mutate(as.data.frame(tmp),ensembl_gene_id=rownames(tmp))
  tmp<-full_join(tmp,gene_ids)
  print(length(unique(tmp$external_gene_name)))
  tmp<-tmp[!duplicated(tmp$external_gene_name),]
  tmp[is.na(tmp)]<-0 #make all empty value to zero
  rownames(tmp)<-tmp$external_gene_name
  tmp<-dplyr::select(tmp,-external_gene_name,-ensembl_gene_id)
  return(tmp)
}
