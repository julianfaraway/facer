% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/facer.R
\name{arraymean}
\alias{arraymean}
\title{Array mean}
\usage{
arraymean(A)
}
\arguments{
\item{A}{a 3D array}
}
\value{
a matrix with dimensions equal to the first two dimensions of A
}
\description{
computes array mean across the third dimension of an array
}
\details{
surprisingly faster than using apply because colMeans is fast
}
