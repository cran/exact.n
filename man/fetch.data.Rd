% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch.data.r
\name{fetch.data}
\alias{fetch.data}
\title{Download a target power library}
\usage{
fetch.data(alpha, delta, prin = FALSE)
}
\arguments{
\item{alpha}{value of nominal size of test}

\item{delta}{value of clinically relevant difference}

\item{prin}{If TRUE, error and/or warning messages are printed.}
}
\value{
No return value. A library object will appear in the global environment.
}
\description{
Function downloads one of 85 power libraries from chrislloyd.com.au.
The libraries are all four column matrices with roughly 20 million
rows and will be around 1.4Gb within R. Download should take roughly
a minute. The object will have a name of the form LIB.alpha.delta.Rdata.
}
\examples{
#'
\dontrun{
fetch.data(alpha=0.05,delta=0.10)
# Object LIB.a050.d10 will appear with 83,289,960 elements
}

}
\references{
C.J. Lloyd (2022) Exact samples sizes for clinical trials subject to size and power constraints. \doi{10.13140/RG.2.2.11828.94085}
}
\author{
Chris J. Lloyd
}
