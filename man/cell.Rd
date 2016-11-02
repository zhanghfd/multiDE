\name{cell}
\alias{cell}
\docType{data}
\title{cellular phenotype data}
\description{
a study of Homo sapiens hormone embryonic stem cells was downloaded from the NCBI GEO database (accession ID: GSE36552). To find the casual relationship between gene expression network and cellular phenotype, Yan et al. derived embryonic stem cells from donated human pre-implantation embryos, prepared cDNA and sequenced them by Illumina HiSeq2000.}
\usage{data(cell)}
\format{
  A data frame with 6,526 observations on the following 10 variables.
  \describe{
    \item{\code{nameOfGene}}{Genes name}
    \item{\code{2cell_e1}}{Read count for sample 1 of 2-cell stage}
    \item{\code{2cell_e2}}{Read count for sample 2 of 2-cell stage}
    \item{\code{2cell_e3}}{Read count for sample 3 of 2-cell stage}
    \item{\code{4cell_e1}}{Read count for sample 1 of 4-cell stage}
    \item{\code{4cell_e2}}{Read count for sample 2 of 4-cell stage}
    \item{\code{4cell_e3}}{Read count for sample 3 of 4-cell stage}
    \item{\code{8cell_e1}}{Read count for sample 1 of 8-cell stage}
    \item{\code{8cell_e2}}{Read count for sample 2 of 8-cell stage}
    \item{\code{8cell_e3}}{Read count for sample 3 of 8-cell stage}
  }
}
\references{
Yan L., Yang M., Guo H., Yang L., Wu J., Li R., Liu P., Lian Y., Zheng X., Yan J., et al. (2013). Single-cell RNA-Seq profiling of human preimplantation embryos and embryonic stem cells. Nature Structural & Molecular Biology 20(9), 1131-1139.
}
\keyword{datasets}
