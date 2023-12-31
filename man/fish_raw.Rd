% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fish_raw.R
\docType{data}
\name{fish_raw}
\alias{fish_raw}
\title{Fish eDNA metabarcoding dataset}
\format{
A list containing an array of sequence read count \code{y},
a vector of the total number of mismatched bases in the priming region of
the forward and reverse primers for each species \code{mismatch},
and a factor indicating whether the riverbank of each site lacks aquatic
and riparian vegetation \code{riverbank}.
The sequence reads were obtained from three replicates (collected from the
center of the river and near the left and right riverbanks, respectively)
from 50 sites across the watershed, of which read counts from 6 samples are
missing.
The resulting sequence counts of the 50 freshwater fish taxa detected are
recorded.
}
\source{
K. Fukaya, N. I. Kondo, S. S. Matsuzaki, T. Kadoya (2021) Data from: Multispecies site occupancy modeling and study design for spatially replicated environmental DNA metabarcoding. Dryad Digital Repository. \doi{10.5061/dryad.3bk3j9kkm}
}
\usage{
fish_raw
}
\description{
A dataset of fish eDNA metabarcoding collected in the Kasumigaura watershed, Japan.
}
\keyword{datasets}
