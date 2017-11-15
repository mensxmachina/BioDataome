#map a disease to DO-level2

diseasetoDO2<-function(x){
#diseaseCategoryALLU<-get(load("E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/diseaseCategoryALLU.Rda"))
do2<-BioDataome::diseaseCategoryALLU[match(tolower(x),tolower(BioDataome::diseaseCategoryALLU[,1])),4]
return(do2)
}
BioDataome::diseaseCategoryALLU

t<-diseasetoDO2(x)
