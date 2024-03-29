% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/champ.normv2.R
\name{champ.normv2}
\alias{champ.normv2}
\title{champ.normv2}
\usage{
champ.normv2(
  beta = myLoad$beta,
  rgSet = myLoad$rgSet,
  mset = myLoad$mset,
  resultsDir = "./CHAMP_Normalization/",
  method = "BMIQ",
  plotBMIQ = FALSE,
  arraytype = "mouse",
  cores = 4
)
}
\arguments{
\item{beta}{Beta Matrix}

\item{rgSet}{RGset (optional)}

\item{mset}{Mset (optional)}

\item{resultsDir}{directory with results}

\item{method}{type of normalisation, defaults to BMIQ}

\item{plotBMIQ}{should output be plotted? F by default}

\item{arraytype}{type of array.}

\item{cores}{used cores, default = 4}
}
\value{
Normalised beta matrix
}
\description{
ChAMP normalisation for the Illumina Human Methylation EPIC v2 array. Like the ChAMP champ.norm package, it can use SWAN, funNorm, PBC, or BMIQ normalisation. Depending on normalisation algorithm, different inputs may be needed (beta matrix, rgSet, mset). Please refer to the ChAMP package manual.
}
