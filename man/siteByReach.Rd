% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topology.r
\name{siteByReach}
\alias{siteByReach}
\title{Produce a site by reach connectivity matrix}
\usage{
siteByReach(ws, points, names, self = TRUE)
}
\arguments{
\item{ws}{A watershed object}

\item{points}{A vector of pixel id numbers for sites}

\item{names}{An optional vector of site names}

\item{self}{If TRUE, a reach is considered connected to itself}
}
\value{
A [Matrix::sparseMatrix()]
}
\description{
A value of 1 at indices `[i,j]` indicates that reach `i` is connected to (i.e., downstream of)
reach `j`.
}
