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

### FCS
FCS generates a matrix of p-values based on significant enrichment of observed proteins against a vector of complexes. It takes a data matrix (RCC_original) and a vector of complex features as its primary inputs. Sim_size is the number of simulations and should  be set to 1000 typically. Threshold is the minimal complex size to consider (default is usually size 5). 

      fcs(RCC_original, complex_vector,sim_size=1000, threshold=5)

Note that FCS can take a while to run, especially if there are many samples, and a large feature vector to consider. A example output is provided in fcs_rc_n.
