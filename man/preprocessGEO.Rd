% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessGEO.R
\name{preprocessGEO}
\alias{preprocessGEO}
\title{Preprocess CEL files with SCAN}
\usage{
preprocessGEO(x, y)
}
\arguments{
\item{x}{the path where the CEL files are stored}

\item{y}{the number of cores to run in parallel}
}
\value{
a matrix of dimensions: probes x samples with the normalized expression values
}
\description{
This function calls the SCAN method as described in Piccolo SR, Sun Y, Campbell JD, Lenburg ME,
Bild AH and Johnson WE (2012). A single-sample microarray normalization method to facilitate personalized-medicine workflows.
Genomics, 100(6), pp. 337-344.
Example: Assuming that CEL files are located in working directory
}
\examples{
\dontrun{
preprocessGEO(getwd(),3) }
normCELs<-BioDataome:::normalized
}
