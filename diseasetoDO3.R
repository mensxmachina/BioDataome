#map a disease to DO-level2

diseasetoDO3<-function(x){
  diseaseSubCategoryALLU<-get(load("E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/diseaseSubCategoryALLU.Rda"))
  do3<-diseaseSubCategoryALLU[match(tolower(x),tolower(diseaseSubCategoryALLU[,1])),4]
return(do3)
}
