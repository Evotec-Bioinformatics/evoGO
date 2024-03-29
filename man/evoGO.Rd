% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evogo-class.R
\name{evoGO}
\alias{evoGO}
\title{Create an input evoGO object for the gene set analysis}
\usage{
evoGO(graphTable, geneSets, annotation = NULL, minGenes = 3, nCores = 1)
}
\arguments{
\item{graphTable}{a data.frame defining relationship between the GO terms.
The first column must be called "child" and the second is "parent".}

\item{geneSets}{a named list with character vectors of gene groups per GO term,
e.g. list(`GO:01` = c("ENS01", "ENS02"), `GO:02` = c("ENS03")), where
the names are the GO term identifiers.
The list items may contain the genes belonging to all children terms or
include only the genes with exclusion of the children terms.}

\item{annotation}{a data.frame with first column "id" with the term identifiers
and columns with the annotation information.}

\item{minGenes}{minimal size of a term to keep in the object.}

\item{nCores}{number of cores to use.}
}
\description{
Create an input evoGO object for the gene set analysis
}
