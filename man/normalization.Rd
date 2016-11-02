

\name{normalization}
\alias{normalization}

\title{
Normalize the raw RNA-seq read count data.
}
\description{
This function normalizes raw count data. The imputs include a raw count data matrix and the normalization method. There are four options of normalizations methods, namely meandian normalization (median), total normalization (total), quantile normalization (quantile), and trimmed mean of M-values (TMM). The outputs is a list containing the 
estimated size factor and normalized count matrix.
}
\usage{
normalization(count,method="median")
}

\arguments{
  \item{count}{Raw RNA-seq read count data matrix, rows for genes and columns for samples.
}
  \item{method}{
The normalization method to be used. This should be one of "median", "total", "quantile", or "TMM". The default normalization method is "median". 
}
}
\value{

 \item{count}{Normalized read counts matrix, rows for genes and columns for samples.}

 \item{sizeFactor}{Estimated size factor.}
 
}
\references{
Kang G, Du L, Zhang H (2016). multiDE: A dimension reduced model based statistical method for differential expression analysis using RNA-sequencing data with multiple treatment conditions.
}
\author{
Guangliang Kang and Hong Zhang
}
\examples{

data(psoriatic);

count = as.matrix(psoriatic[,-1]);

M = normalization(count,method="median");

}


\keyword{Normalization}



