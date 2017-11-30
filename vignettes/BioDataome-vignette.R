## ----style, echo = FALSE, results = 'asis'---------------------------------
BiocStyle::markdown()

## ----setup, include = FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## --------------------------------------------------------------------------
#install.packages("devtools")

## ---- results= "hide"------------------------------------------------------
library(devtools)

## --------------------------------------------------------------------------
#install_github("mensxmachina/DBIR", subdir="BioDataome")

## ---- warning=FALSE,message=FALSE------------------------------------------
library(BioDataome)

## --------------------------------------------------------------------------
#Download the filtered list matching specific criteria

metadata<-read.csv(file="E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/BioDataome/metadata.csv", sep=",", header=TRUE, stringsAsFactors = F)

#Let's filter out further the results by selecting datasets with sample size greater or equal to 25 and those with no common samples to any other dataset.
subsetDsets<-metadata[which( (metadata$samples>=25) & (metadata$duplicates==" ") ),]

dataList<-lapply(subsetDsets$file,read.csv)

#first five rows and columns of the data
dataList[[1]][1:5,1:5]

## ----eval=FALSE, warning=FALSE---------------------------------------------
#  #curate dataset GSE10026 measured with technology GPL570.
#  
#  GSE10026<-curateGSE("GSE10026","GPL570",getwd())
#  
#  #first five rows and columns of the metadata
#  GSE10026[[1]][1:5,1:5]
#  
#  #first five rows and columns of the preprocessed data
#  GSE10026[[2]][1:5,1:5]
#  

## ----eval=FALSE, warning=FALSE---------------------------------------------
#  #curate dataset SRP032775.
#  
#  SRP032775<-curateRecountRNASeq("SRP032775",getwd())
#  
#  #first five rows and columns of the metadata
#  SRP032775[[1]][1:5,1:5]
#  
#  #first five rows and columns of the preprocessed data
#  SRP032775[[2]][1:5,1:5]

## --------------------------------------------------------------------------
#download dataset SRP050971

downloadRecount("SRP050971")
load(file.path("SRP050971", 'rse_gene.Rdata'))
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

## ---- fig.height=5,fig.width=7---------------------------------------------
#source("https://bioconductor.org/biocLite.R")
#biocLite("geneplotter")
library('geneplotter')

multidensity(SummarizedExperiment::assay(dds),
              xlab="mean counts", xlim=c(0, 1000),
              ylab="Density", legend=F, main="")

#account for heteroscedasticity
vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)
dataNorm<-SummarizedExperiment::assay(vsd)
#show the effect of the variance stabilizing transformation, by plotting the first sample against the second
lims <- c(-2, 20)
plot(SummarizedExperiment::assay(dds)[,1:2],
     pch=16, main="before VST")

lims <- c(-2, 20)
plot(dataNorm[,1:2],
       pch=16, main="after VST")

## ---- fig.height=5,fig.width=7---------------------------------------------
#install.packages("RColorBrewer")
library(RColorBrewer)

distsRL <- dist(t(dataNorm))
mat <- as.matrix(distsRL)
condition<-SummarizedExperiment::colData(vsd)$title
condition<-strsplit(condition,"_")
condition<-sapply(condition,"[[",2)

rownames(mat) <-  condition
colnames(mat) <-  SummarizedExperiment::colData(vsd)$sample

hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
#install.packages("pheatmap")
library(pheatmap)
pheatmap(mat,fontsize_row=8)

## --------------------------------------------------------------------------
#download two preprocessed datasets from BioDataome in .rda
d1<-get(load(url("http://dataome.mensxmachina.org/data/Homo%20sapiens/GPL570/GSE86013.Rda")))
d2<-get(load(url("http://dataome.mensxmachina.org/data/Homo%20sapiens/GPL570/GSE86015.Rda")))
# compare two preprocessed datasets
commons<-compareDsets(d1,d2)
#Number of samples that d1 and d2 share
commons

## --------------------------------------------------------------------------
#the path to GSE8601 dataset in BioDataome
x<-"http://dataome.mensxmachina.org/data/Homo%20sapiens/GPL570/GSE86013.csv"

## --------------------------------------------------------------------------
#create a character vector of the paths in BioDataome
y<-c("GSE86015.csv","GSE9008.csv","GSE9119.csv")
y<-paste0("http://dataome.mensxmachina.org/data/Homo%20sapiens/GPL570/",y)
#find with which of the three datasets our dataset of interest shares samples
commonGSEs<-compareDsetList(x,y)
commonGSEs

## --------------------------------------------------------------------------
#Assign Disease-Ontology terms to study "GSE10245"
diseases<-GSEtoDisease("GSE10006")
diseases

## --------------------------------------------------------------------------
#Assign Disease-Ontology terms to study "SRP032775"
gse<-recountIDtoGSE("SRP032775")
diseases<-GSEtoDisease(gse)
diseases

## --------------------------------------------------------------------------
#download GSE8671 preprocessed dataset from BioDataome in .rda
d<-get(load(url("http://dataome.mensxmachina.org/data/Homo%20sapiens/GPL570/GSE8671.Rda")))

#download GSE8671 preprocessed dataset from BioDataome in .rda
pheno<-get(load(url("http://dataome.mensxmachina.org/data/Homo%20sapiens/GPL570/GSE8671_Annot.Rda")))
#dimensions of dataset (rows: probes, columns:samples)
dim(d)

#install and load limma package

#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
library(limma)
#Fit a linear model for each gene in the expression data given the design matrix specified in column class of the phenotype data
fit <- lmFit(d, design=model.matrix(~ pheno$class))
#calculate the differential expression by empirical Bayes shrinkage
fit <- eBayes(fit)
#A table with the top 10 most statistically significant differentially expressed genes between the groups sorted by adjusted p-value:
tt <- topTable(fit, coef=2)

#create a heatmap of those top 10 highly significant genes. 
heatmap(d[match(rownames(tt),row.names(d)),], labCol = FALSE)

#Find gene symbols of those top 10 highly significant genes.

#install and load necessary packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("hgu133plus2.db")
#biocLite("annotate")

library("annotate")
library('hgu133plus2.db')

Symbol <- getSYMBOL(row.names(d), "hgu133plus2.db")
#Find gene symbols for the top 10 differentially expressed genes 
top10genes<-Symbol[match(rownames(tt),row.names(d))]
top10genes

## --------------------------------------------------------------------------

#Perform enrichment analysis for the top 10 differentially expressed genes 

#install and download enrichR package
#install.packages("enrichR")
library('enrichR')

#select annotated gene set libraries
dbs <- c("GO_Molecular_Function_2017", "GO_Cellular_Component_2017", "GO_Biological_Process_2017",
         "Reactome_2016","KEGG_2016")
#Perform enrichment analysis
enriched <- enrichr(top10genes, dbs)

#sort Gene Ontology molecular function terms by adjusted p-values
GOMF<-enriched[[1]]$Term[order(enriched[[1]]$Adjusted.P.value)]
#show top-5 
head(GOMF, n=5)
#sort genes by the number of Gene Ontology molecular function terms
GeneMembershipGOMF<-sort(table(enriched[[1]]$Genes),decreasing = T)
GeneMembershipGOMF

#sort Gene Ontology cellular component terms by adjusted p-values
GOCC<-enriched[[2]]$Term[order(enriched[[2]]$Adjusted.P.value)]
#show top-5 
head(GOCC, n=5)
#sort genes by the number of Gene Ontology cellular component terms
GeneMembershipCC<-sort(table(enriched[[2]]$Genes),decreasing = T)
GeneMembershipCC

#sort Gene Ontology biological process terms by adjusted p-values
GOBC<-enriched[[3]]$Term[order(enriched[[3]]$Adjusted.P.value)]
#show top-5 
head(GOBC, n=5)
#sort genes by the number of Gene Ontology biological process terms
GeneMembershipBC<-sort(table(enriched[[2]]$Genes),decreasing = T)
GeneMembershipBC

#sort Reactome pathway terms by adjusted p-values
Reactome<-enriched[[4]]$Term[order(enriched[[4]]$Adjusted.P.value)]
#show top-5 
head(Reactome, n=5)
#sort genes by the number of Reactome pathway terms
GeneMembershipReactome<-sort(table(enriched[[4]]$Genes),decreasing = T)
GeneMembershipReactome

#sort KEGG pathway terms by adjusted p-values
KEGG<-enriched[[5]]$Term[order(enriched[[5]]$Adjusted.P.value)]
#show top-5 
head(KEGG,n=5)
#sort genes by the number of GKEGG pathway terms
GeneMembershipKEGG<-sort(table(enriched[[5]]$Genes),decreasing = T)
GeneMembershipKEGG

