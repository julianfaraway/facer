% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/facer.R
\name{bridgefix}
\alias{bridgefix}
\title{Location shift onto template}
\usage{
bridgefix(templface, face, marker = 4, meanfix = FALSE)
}
\arguments{
\item{templface}{the template face}

\item{face}{the face to be adjusted}

\item{marker}{the marker on the template to be used as the anchor}

\item{meanfix}{logical indicating whether to use the mean rather than a specific marker.}
}
\value{
the adjusted face
}
\description{
Face is location-adjusted to make selected marker coincide with the template
}
\details{
Usually do this with some (relatively) fixed landmark like the bridge of the nose
}
