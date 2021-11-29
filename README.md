PROTREC
===
`PROTREC` is an R package containing several functions for predicting and validating missing proteins in proteomics data based on Kong W, Wong B J H, Gao H, et al. PROTREC: A probability-based approach for recovering missing proteins based on biological networks. Journal of Proteomics, 2022, 250: 104392.

`PROTREC` methods are currently divisible into five functional classes
- Protein Recovery (PROTREC)
- Functional Class Scoring (FCS)
- Hypergeometic Enrichment (HE)
- Gene Set Enrichment Analysis (GSEA)
- Performance metrics (Recovery Rate)

## Getting Started

##Example data included in the package
First we need to find some network information to act as the feature vector. PROTREC works well with real complexes and this data can be obtained from CORUM (http://mips.helmholtz-muenchen.de/genre/proj/corum/).

An example complex dataset (complex_vector) is available, which is already processed CORUM complex 2018 release. 

    library(PROTREC)
    complexes <- data(complex_vector)


Alternatively, PROTREC can also use as its feature vector a list of predicted network clusters or pathways. 

One proteomics expression datasets are provided. The renal cancer dataset (RC) comprises 12 normal (RC_N) and 12 cancer (RC_C) samples. Both datasets may be called by their names (in brackets).

A peptide based validation file for RC (RC_peptides_uniq) is available. This is akin to a list of unique and ambiguous PSMs that can be used for checking if there is at least one peptide that points to the presence of a predicted missing protein. The protein and peptide file can be processed as follows:

    RC_cancer <- data(RC_C)
    RC_normal <- data(RC_N)
    peptide_support <- data(RC_peptides_uniq)
    

## Protein Reovery Methods

### Protein Recovery (PROTREC)
PROTREC is a probability-based scoring schema for assigning probabilities to observed and predicted missing proteins based on the false discovery rate of the proteomics screen, and the probability its constituent complex (if any) is present in the sample.

### PROTREC complex probability
Assigns a PROTREC-based probability based on a complex vector and a set of observed proteins

      PROTREC_cplx_prob <- function(data, complex, fdr, threshold)

Where data is a data matrix, complex is a list of complexes, fdr is the false discovery rate of the screen, set as 0.01 or 0.05 normally. Threshold is the minimum size of the complex to consider, normally set as 5. 

For example, 

      rc_nprotrec <- PROTREC_cplx_prob(RC_N, complex_vector, 0.01, 5)

#### PROTREC individual protein probability
Assigns a PROTREC-based probability to an individual protein based on a series of complex-based PROTREC probabilities

      PROTREC_protprob <- function(cplx, p, prot, fdr)

Where cplx is the list of complex components, p is the probability complex exists, and prot is the list of proteins in the screen

For example, 

      PROTREC_prot_rc_n_1 <- data.frame(PROTREC_protprob(complex_vector, 1-rc_nprotrec[1,], rownames(RC_N),0.01))

### FCS

#### FCS complex probability

FCS generates a matrix of p-values based on significant enrichment of observed proteins against a vector of complexes.

      fcs <- function(data, complex_vector, sim_size, threshold)

It takes a data matrix (For example, RC_N) and a vector of complex features (For example, complex_vector) as its primary inputs. Sim_size is the number of simulations and should be set to 1000 typically. Threshold is the minimal complex size to consider (default is usually size 5). 

      rcnfcs <- fcs(RC_N,complex_vector,1000,5)

Note that FCS can take a while to run, especially if there are many samples, and a large feature vector to consider.

#### FCS individual protein probability

fcs_prot_prob_assign assigns probabilities to individual proteins based on the FCS probability.

      fcs_prot_prob_assign <- function(cplx, p)
      
Where cplx is the complex vector and p is a vector of complex-based probabilities derived from FCS. Since FCS provides p-values, then p is simply (1 - FCS p-values). For example:

      fcs_prot_rc_n_1 <- prot_prob_fcs(complex_vector, 1 - rcnfcs[1,])
