#' Run all steps to download, preprocess and annotate an RNASeq dataset from Recount
#'
#' @param x a recount dataset ID
#' @param y the path to write the output
#' @return writes in the given path two data frames, the preprocessed data and the metadata file with phenotype information
#' @examples
#' curatedRecount<-curateRecountRNASeq("SRP032775",getwd())
#' @export
#' @importFrom recount scale_counts
#' @importFrom DESeq2 DESeqDataSetFromMatrix varianceStabilizingTransformation estimateSizeFactors
#' @importFrom SummarizedExperiment assay

curateRecountRNASeq<-function(x,y){
  if (missing(x))
    stop("Need to specify a recount accession id, i.e 'SRP032775'")
  if ( (missing(y)) | (!file.exists(y)) )
    stop("Need to specify a valid path, i.e getwd()")

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
  #Estimate the size factors for the count data
  dds <- DESeq2::estimateSizeFactors(dds)
  #account for heteroscedasticity
  vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)
  #normalized data matrix
  dataNorm<-SummarizedExperiment::assay(vsd)
  samples<-colnames(dataNorm)
  gse<-recountIDtoGSE(x)
  if(!is.na(gse)){
    metadata<-GSEmetadata(gse,"GPL11154")
    #match samples in data matrix
    sids<-match(pheno$geo_accession,metadata$samples)
    metadata<-metadata[sids,]
    rownames(metadata)<-samples
  } else {
    metadata<-as.data.frame(pheno)
    #just for consistency write columns samples and class
    metadata<-cbind(samples=samples,class="unknown",metadata)
  }

  curatedRecount<-list(metadata=metadata,dataNorm=dataNorm)

  return(curatedRecount)

}
