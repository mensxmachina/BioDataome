#' Run all steps to download, preprocess and annotate an RNASeq dataset from Recount
#'
#' @param x a recount dataset ID
#' @param y the path to write the output
#' @return writes in the given path two data frames, the preprocessed data and the metadata file with phenotype information
#' @examples
#' curateRecountRNASeq("SRP032775",getwd())
#' @export
#' @importFrom recount scale_counts
#' @importFrom DESeq2 DESeqDataSetFromMatrix varianceStabilizingTransformation estimateSizeFactors
#' @importFrom SummarizedExperiment assay

curateRecountRNASeq<-function(x,y){
  #define geometric means function
  gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
  #download and load the data matrix
  downloadRecount(x)
  load(file.path(x, 'rse_gene.Rdata'))
  # Scale counts by taking into account the total coverage per sample
  rse <- recount::scale_counts(rse_gene)
  #access count matrix
  data<-SummarizedExperiment::assay(rse)
  #access phenotype information
  pheno<-SummarizedExperiment::colData(rse)
  #construct a DESeqDataSet object for further analysis
  dds<-DESeq2::DESeqDataSetFromMatrix(data,pheno, ~ 1)
  #find genometric mean for the count data
  geoMeans = apply(SummarizedExperiment::assay(dds), 1, gm_mean)
  #Estimate the size factors for the count data
  dds <- DESeq2::estimateSizeFactors(dds, geoMeans=geoMeans)
  #account for heteroscedasticity
  vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)
  #normalized data matrix
  dataNorm<-SummarizedExperiment::assay(vsd)
  save(dataNorm, file=paste0(y,"/",x,".Rda"))

  samples<-colnames(dataNorm)
  gse<-recountIDtoGSE(x)
  if(!is.na(gse)){
    metadata<-GSEmetadata(gse,"gpl11154")
    #match samples in data matrix
    sids<-match(pheno$geo_accession,metadata$samples)
    metadata<-metadata[sids,]
    rownames(metadata)<-samples
  } else {
    metadata<-as.data.frame(pheno)
    #just for consistency write columns samples and class
    metadata<-cbind(samples=samples,class="unknown",metadata)
  }

  save(metadata, file=paste0(y,"/",x,"_Annot.Rda"))


}
