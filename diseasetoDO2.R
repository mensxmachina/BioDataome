#map a disease to DO-level2

diseasetoDO2<-function(x){
diseaseCategoryALLU<-get(load("E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/diseaseCategoryALLU.Rda"))
do2<-diseaseCategoryALLU[match(tolower(x),tolower(diseaseCategoryALLU[,1])),4]
return(do2)
}
