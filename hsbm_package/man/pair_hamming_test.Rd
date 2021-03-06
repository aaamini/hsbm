% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation.R
\name{pair_hamming_test}
\alias{pair_hamming_test}
\title{Testing in-group versus random pairs}
\usage{
pair_hamming_test(A, grps, npairs = 5000, type = 1)
}
\arguments{
\item{A}{a list of (sparse) adj. matrices -- should be a list; assuming all
networks have the same number of nodes}

\item{grps}{grps a list of indices within each group}
}
\description{
Testing in-group versus random pairs
}
