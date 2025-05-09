% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PROTREC.R
\name{pairwise_recovery}
\alias{pairwise_recovery}
\title{This works out the recovery significance for FCS,HE and GSEA}
\usage{
pairwise_recovery(
  predict_list,
  original_prot_list,
  check_prot_list,
  complex_vec
)
}
\arguments{
\item{predict_list}{A vector of significant proteins for a given sample, default set 0.05 p-value as significant}

\item{original_prot_list}{The original set of proteins observed for that particular sample}

\item{check_prot_list}{A set of proteins observed in a second replicate}

\item{complex_vec}{A list of complex objects}
}
\value{
output A vector of five pieces of information: the observed proporetion of overlap,
the significance of this overlap (p-value), the number of verified proteins, the total number of predicted missing proteins,
and the list of validated proteins separated by 'a'.
}
\description{
This works out the recovery significance for FCS,HE and GSEA
}
