% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PROTREC.R
\name{fcs}
\alias{fcs}
\title{This assigns probabilities to complexes based on FCS}
\usage{
fcs(data, complex_vector, sim_size, threshold)
}
\arguments{
\item{data}{Proteomics expression data}

\item{complex_vector}{protein complex}

\item{sim_size}{how many iterations for simulation, recommend 1000}

\item{threshold}{the minimum size of the complex required}
}
\value{
output_mat A matrix of FCS probabilities assigned to complexes for the dataset
}
\description{
This assigns probabilities to complexes based on FCS
}
