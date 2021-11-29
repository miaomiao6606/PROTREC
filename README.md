PROTREC
===
`PROTREC` is an R package containing several functions for predicting and validating missing proteins in proteomics data based on Kong W, Wong B J H, Gao H, et al. PROTREC: A probability-based approach for recovering missing proteins based on biological networks. Journal of Proteomics, 2022, 250: 104392.

`PROTREC` methods are currently divisible into five functional classes
- Protein Recovery (PROTREC)
- Functional Class Scoring (FCS)
- Hypergeometic Enrichment (HE)
- Gene Set Enrichment Analysis (GSEA)
- Performance metrics

## Getting Started

## Example data included in the package
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

      rc_nfcs <- fcs(RC_N,complex_vector,1000,5)

Note that FCS can take a while to run, especially if there are many samples, and a large feature vector to consider.

#### FCS individual protein probability

fcs_prot_prob_assign assigns probabilities to individual proteins based on the FCS probability.

      fcs_prot_prob_assign <- function(cplx, p)
      
Where cplx is the complex vector and p is a vector of complex-based probabilities derived from FCS. Since FCS provides p-values, then p is simply (1 - FCS p-values). For example:

      fcs_prot_rc_n_1 <- data.frame(prot_prob_fcs(complex_vector, 1 - rcnfcs[1,]))

### HE

#### HE complex probability

HE generates a matrix of p-values based on significant enrichment of observed proteins against a vector of complexes.

      hgtest <- function(data, complex_vector, threshold)

It takes a data matrix (For example, RC_N) and a vector of complex features (For example, complex_vector) as its primary inputs. Threshold is the minimal complex size to consider (default is usually size 5). 

      rc_nhg<- data.frame(hgtest(RC_N, complex_vector,5))

#### HE individual protein probability

hgtest_prot_prob_assign assigns probabilities to individual proteins based on the HE probability.

      hgtest_prot_prob_assign <- function(cplx, p)
      
Where cplx is the complex vector and p is a vector of complex-based probabilities derived from HE. Since HE provides p-values, then p is simply (1 - HE p-values). For example:

      rc_nhg <- data.frame(hgtest_prot_prob_assign(complex_vector,1-rc_nhg[1,] ))

### GSEA

#### GSEA complex probability

GSEA generates a matrix of p-values based on significant enrichment of observed proteins against a vector of complexes.

      repgsea <- function(data, complex_vector)

It takes a data matrix (For example, RC_N) and a vector of complex features (For example, complex_vector) as its primary inputs. 
      rc_ngsea_p <- repgsea(RC_N,complex_vector)

#### GSEA individual protein probability

hgtest_prot_prob_assign assigns probabilities to individual proteins based on the GSEA probability.

      gsea_prot_prob_assign <- function(cplx, p)
      
Where cplx is the complex vector and p is a vector of complex-based probabilities derived from GSEA. Since GSEA provides p-values, then p is simply (1 - GSEA p-values). For example:

      rc_ngsea<- data.frame(gsea_prot_prob_assign(complex_vector,1-rc_ngsea_p[1,] ))

## Performance Metrics

### Recovery rate for FCS, HE and GSEA
Evaluates the significance of the proportion of verified proteins given a set of predicted proteins.

      pairwise_recovery <- function(predict_list, original_prot_list, check_prot_list, complex_vec)

Where predict_list is a vector of predicted(by FCS, HE or GSEA) significant values for a given sample, original_prot_list is the original set of proteins observed for that particular sample, check_prot_list is the set of proteins in the cross-replicate for verification, and complex_vector is the feature vector set. The function take p-value <0.05 as default setting for judging significance. For example:

      pairwise_recovery(t(rc_nfcs)[,1],rownames(RC_N)[(RC_N[,1])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,1])!=0],complex_vector)
      
The output will contain a vector of five pieces of information: the observed proporetion of overlap, the significance of this overlap (p-value), the number of verified proteins, the total number of predicted missing proteins, and the list of validated proteins separated by the character 'a'.

### Recovery rate for PROTREC
Evaluates the significance of the proportion of verified proteins given a set of predicted proteins.

      pairwise_recovery_protrec <- function(prot_predict_list, original_prot_list, check_prot_list, complex_vec,protrecscoreset=0.95)

Where prot_predict_list is a vector of predicted significant values by PROTREC for a given sample, original_prot_list is the original set of proteins observed for that particular sample, check_prot_list is the set of proteins in the cross-replicate for verification, complex_vector is the feature vector set, and the protrecscoreset is the PROTREC score significance cutoff, default is 0.95. For example:

      PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-rc_cprotrec[1,], rownames(RC_C),0.01))
      rownames(PROTREC_prot)<- PROTREC_prot[,1]
      PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
      PROTREC_prot<- data.frame(PROTREC_prot)
      pairwise_recovery_protrec(PROTREC_prot,rownames(RC_C)[(RC_C[,1])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,1])!=0],complex_vector)
      
The output will contain a vector of five pieces of information: the observed proporetion of overlap, the significance of this overlap (p-value), the number of verified proteins, the total number of predicted missing proteins, and the list of validated proteins separated by the character 'a'.

