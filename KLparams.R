#a function that calculates all parameters needed in the KL computation
#input:the output of microarray normalization, namely
#a dataset with #probes x #samples
#and t as the threshold of variance explained (in % i.e 0.5)

KLparams<-function(dataset,th){
  #make rows:samples, columns:probe sets
  dataset<-t(dataset)
  #scale it
  dataset = scale(dataset)
  
  #if zero columns exist
  #dataset[is.nan(dataset)] <- 0
  #compute PC
  pc<-prcomp(dataset[,1:ncol(dataset)])
  #number of PC that explain variance greater than t
  npc<-min(which(stats:::summary.prcomp(pc)$importance[3, ]>=th))
  N<-ncol(dataset) #number of variables
  #number of total PCs
  npcT<-length(pc$sdev)
  #sum of the reamining variance
  sigma<-sum( (stats:::summary.prcomp(pc)$importance[1, (npc+1):npcT])^2 )/N
  #loadings of the selected PCs
  pcL = pc$rotation[,1:npc]
  #variance of the selected PCs
  pcvars<-stats:::summary.prcomp(pc)$importance[1, 1:npc]^2
  #output
  params<-list(npc,npcT,sigma,pcL,pcvars)
  names(params)<-c("npc","npcT","sigma","pcL","pcvars")
  return(params)

}