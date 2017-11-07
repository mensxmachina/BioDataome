#this function downloads raw CEL files from GEO from a specified study
#It receives as input a GEO Series (GSE) id, i.e."GSE11761" and
#the path to save the compressed downloaded file.
downloadRaw <- function(x,y) {
  raw <- NULL
  a <- TRUE
  while(a) {
    raw <- try(getGEOSuppFiles(x, makeDirectory = F, baseDir = y))
    if(!inherits(raw,"try-error")) a <- FALSE
    else cat("NCBI was busy, trying again\n")
  }
  return(raw)
}
