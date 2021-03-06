% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/facer.R
\name{makesym}
\alias{makesym}
\title{Make a face symmetrical}
\usage{
makesym(sm, mfconfig)
}
\arguments{
\item{sm}{a face matrix}

\item{mfconfig}{face config info including left right pairs}
}
\value{
scores for each marker
}
\description{
Make a face symmetrical by averaging with reflection
}
\details{
Face is reflected and averaged with reflection
}
