\name{psoriatic}
\alias{psoriatic}
\docType{data}
\title{psoriatic study data}
\description{
In this study, the major interest was to detect the influence of aryl hydrocarbon receptor (AhR) on RNA expression profiles of psoriatic lesion cells. Each of eight patients were treated with culture treatment of DMSO (vehicle control), AhR agonist FICZ, and AhR antagonist CH-2233191. RNA-seq data were obtained using Illumina Genome Analyzer II platform for each of three treated lesion tissue samples. Therefore, this was a matched sample design. The RNA-seq read counts data were derived from the GEO database  (accession ID: GSE47944). There were altogether 13,416 genes with maximal read counts greater than 50 in each treatment condition.}
\usage{data(psoriatic)}
\format{
  A data frame with 13,416 observations on the following 25 variables.
  \describe{
    \item{\code{Gene}}{Genes name}
    \item{\code{Sample_PDM_K1_2}}{Read count for patient 1 in vehicle control group}
    \item{\code{Sample_PDM_K2_2}}{Read count for patient 2 in vehicle control group}
    \item{\code{Sample_PDM_K3_2}}{Read count for patient 3 in vehicle control group}
    \item{\code{Sample_K4.2}}{Read count for patient 4 in vehicle control group}
    \item{\code{Sample_K5.2}}{Read count for patient 5 in vehicle control group}
    \item{\code{Sample_K6_2}}{Read count for patient 6 in vehicle control group}
    \item{\code{Sample_K8_2}}{Read count for patient 7 in vehicle control group}
    \item{\code{Sample_K9_2}}{Read count for patient 8 in vehicle control group}   
    \item{\code{Sample_PDM_K1_3}}{Read count for patient 1 in AhR agonist group}
    \item{\code{Sample_PDM_K2_3}}{Read count for patient 2 in AhR agonist group}
    \item{\code{Sample_PDM_K3_3}}{Read count for patient 3 in AhR agonist group}
    \item{\code{Sample_K4.3}}{Read count for patient 4 in AhR agonist group}
    \item{\code{Sample_K5.3}}{Read count for patient 5 in AhR agonist group}
    \item{\code{Sample_K6_3}}{Read count for patient 6 in AhR agonist group}
    \item{\code{Sample_K8_3}}{Read count for patient 7 in AhR agonist group}
    \item{\code{Sample_K9_3}}{Read count for patient 8 in AhR agonist group}    
    \item{\code{Sample_PDM_K1_4}}{Read count for patient 1 in AhR antagonist group}
    \item{\code{Sample_PDM_K2_4}}{Read count for patient 2 in AhR antagonist group}
    \item{\code{Sample_PDM_K3_4}}{Read count for patient 3 in AhR antagonist group}
    \item{\code{Sample_K4.4}}{Read count for patient 4 in AhR antagonist group}
    \item{\code{Sample_K5.4}}{Read count for patient 5 in AhR antagonist group}
    \item{\code{Sample_K6_4}}{Read count for patient 6 in AhR antagonist group}
    \item{\code{Sample_K8_4}}{Read count for patient 7 in AhR antagonist group}
    \item{\code{Sample_K9_4}}{Read count for patient 8 in AhR antagonist group}
  }
}
\references{
Di Meglio P, Duarte J.H., Ahlfors H., Owens N.D., Li Y.,
Villanova F., Tosi I., Hirota K., Nestle F.O., Mrowietz U., et al. (2014).
Activation of the aryl hydrocarbon receptor dampens the severity of
in ammatory skin conditions. Immunity 40(6), 989-1001.}

\keyword{datasets}
