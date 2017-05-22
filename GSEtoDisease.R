#this function takes as input a vector of GSEs (i.e GSE5343)
#and returns a data frame with 3 columns
#1st:GSE code
#2nd: diseaseCategory
#3rd: diseaseSubCategory

#THIS CODE NEEDS MUCH CLEANING

GSEtoDisease<-function(GSE){
  #load neccessary files
  source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/GSEtoDiseaseGEO.R')
  path2<-"E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/"
  diseaseSubCategoryALLU<-get(load("E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/diseaseSubCategoryALLU.Rda"))
  diseaseCategoryALLU<-get(load("E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/diseaseCategoryALLU.Rda"))
  leafs<-get(load("E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/leafs.Rda"))
  annot<-get(load(paste0(path2,"diseasesLevelALLUnique2.Rda")))
  #unlist all gses per disease
  tt<-strsplit(annot$GSEs,";")
  
  pubmed<-c()
  diseasesA<-c()
  diseasesALL<-c()
  
  for (i in 1:length(GSE)){
    diseases<-c()
    #ANNOTATE by PubTator
    #first find gse's pubmedid
    a<-paste0(GSE[i]," [ACCN] AND gse[ETYP]")
    
    r_search<- entrez_search(db="gds", term=a,retmax =10000, use_history=TRUE)
    if (r_search$count!=0){
      b<-entrez_summary(db="gds", id=r_search$ids[1])
      pubmed<-b$pubmedids
      if (length(pubmed)!=0){
        #use PubTator's API to annotate gse
        for (j in 1:length(pubmed)){
          uri<-paste("https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/Disease/",pubmed[j],"/BioC",sep="") 
          t<-getURL(uri)
          ll<-xmlToList(t)
          list1<-ll$document[[3]]
          if (length(list1)>3){
            d<-c()
            for (k in 1:(length(list1)-3)){
              d[k]<-list1[[k+3]]$text   
            }
            diseases<-c(diseases,tolower(d))
          } 
          
        }  
        diseases<-unique(diseases)
        #if PubTator returns an empty set
        #try with GEO
        if (length(diseases)==0){
          diseasesALL[i]<-sapply(GSE[i],GSEtoDiseaseGEO)
        } else {
          #keep first 3 diseases as most relevant
          if (length(diseases)>=3){
            diseases<-diseases[1:3] 
          } 
          #find all leafs
          leaf<-c()
          for (j in 1:length(diseases)){
            leaf[j]<-match(tolower(diseases[j]),tolower(leafs[,1]))
          }
          #if we have only one leaf keep that
          if (length(which(!is.na(leaf)))==1){
            diseasesALL[i]<-diseases[which(!is.na(leaf))]
          } else if (length(which(!is.na(leaf)))>1){
            #if we have more than one leafs find parent categories of all entiites
            #and keep the leaf that belongs to the most common category
            aa<-diseaseSubCategoryALLU[match(tolower(diseases),tolower(diseaseSubCategoryALLU[,1])),4]
            commonCategory<-names(table(aa)[which(table(aa)==max(table(aa)))]) 
            leafCats<-aa[which(!is.na(leaf))]
            #which leafs belong to the most common category
            leafsToKeep<-which(leafCats %in% commonCategory)
            #if no leaf belong to the common category keep them all
            if (length(leafsToKeep)==0){
              diseasesALL[i]<-paste(diseases[which(!is.na(leaf))], sep=";", collapse=";")
            } else {
              diseasesALL[i]<-paste(diseases[leafsToKeep], sep=";", collapse=";")
            }
          } else {
            #if we only have nodes
            #find deepest nodes
            depth<-c()
            for (j in 1:length(diseases)){
              depth[j]<-as.numeric(diseaseSubCategoryALLU[match(tolower(diseases[j]),tolower(diseaseSubCategoryALLU[,1])),3 ] )
            } 
            #if we have more than one category with the same depth, keep the most representative as in leafs
            
            bb<-which(depth==max(depth,na.rm = T))
            if (length(bb)==0){
              diseasesALL[i]<-sapply(GSE[i],GSEtoDiseaseGEO) 
            } else if (length(bb)==1){
              diseasesALL[i]<-diseases[which(depth==max(depth,na.rm = T))]
            } else {
              aa<-diseaseSubCategoryALLU[match(tolower(diseases),tolower(diseaseSubCategoryALLU[,1])),4]
              commonCategory<-names(table(aa)[which(table(aa)==max(table(aa)))])
              nodesToKeep<-which(aa[bb] %in% commonCategory)
              if (length(nodesToKeep)==0){
                diseasesALL[i]<-paste(diseases[bb], sep=";", collapse=";") 
              } else {
                diseasesALL[i]<-paste(diseases[bb[nodesToKeep]], sep=";", collapse=";")
              }
            }
            
          }
        }
        
        if (is.na(diseasesALL[i])){
          diseasesALL[i]<-sapply(GSE[i],GSEtoDiseaseGEO)
        }
        
        
      } else {
        diseasesALL[i]<-sapply(GSE[i],GSEtoDiseaseGEO)
      }

    } else {
      diseasesALL[i]<-NA
    }  
    
  }
  
  towrite<-cbind(GSE,diseasesALL)
  
  towrite2<-strsplit(diseasesALL,";")
  towrite2<-sapply(towrite2,"[[",1)
  towrite2<-cbind(GSE,towrite2)
  
  #map combined to DO3
  DOlevel3<-c()
  for (i in 1:nrow(towrite2)){
    dis<-towrite2[i,2]
    if (!is.na(towrite2[i,2])){
      aa<-diseaseSubCategoryALLU[match(tolower(dis),tolower(diseaseSubCategoryALLU[,1])),4]
      #    if (!is.na(aa)){
      #      commonCategory<-table(aa)
      #      DOlevel3[i]<-names(commonCategory[which(commonCategory==max(commonCategory))])
      DOlevel3[i]<-aa
      #    }
    } else {
      DOlevel3[i]<-NA
    }
    
  }
  
  DOlevel2<-c()
  for (i in 1:nrow(towrite2)){
    dis<-towrite2[i,2]
    if (!is.na(towrite2[i,2])){
      aa<-diseaseCategoryALLU[match(tolower(dis),tolower(diseaseCategoryALLU[,1])),4]
      #    if (!is.na(aa)){
      #      commonCategory<-table(aa)
      #      DOlevel3[i]<-names(commonCategory[which(commonCategory==max(commonCategory))])
      DOlevel2[i]<-aa
      #    }
    } else {
      DOlevel2[i]<-NA
    }
    
  }
  
  toreturn<-data.frame(cbind(towrite2,DOlevel2,DOlevel3))
  names(toreturn) <-c("GSE","disease","DOlevel2","DOlevel3")
  return(toreturn)
}







