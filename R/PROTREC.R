###DATA######
#' @title This is the RC protein data included in my package
#' @name RC
#' @docType data
#' @references \url{http://www.nature.com/nm/journal/v21/n4/full/nm.3807.html}
#'
#'
#'
NULL

#' @title This is the RC cancer protein data included in my package
#' @name RC_C
#' @docType data
#' @references \url{http://www.nature.com/nm/journal/v21/n4/full/nm.3807.html}
#'
#'
#'
NULL

#' @title This is the RC normal protein data included in my package
#' @name RC_N
#' @docType data
#' @references \url{http://www.nature.com/nm/journal/v21/n4/full/nm.3807.html}
#'
#'
#'
NULL

#' @title This is the RC peptide data included in my package
#' @name RC_peptides_uniq
#' @docType data
#' @references \url{http://www.nature.com/nm/journal/v21/n4/full/nm.3807.html}
#'
#'
#'
NULL

#' @title This is the complex_vector included in my package
#' @name complex_vector
#' @docType data
#' @references \url{http://mips.helmholtz-muenchen.de/corum/}
#'
#'
#'
NULL

###PROTREC functions########
#' @title This assigns probabilities to complexes based on PROTREC
#' @name PROTREC_cplx_prob
#' @param data Proteomics expression data
#' @param complex protein complex
#' @param fdr false discovery rate of proteomics screen
#' @param threshold the minimum size of the complex required
#' @return output_mat A matrix of PROTREC (1-probabilities) assigned to complexes for the dataset
#' @export
PROTREC_cplx_prob <- function(data, complex, fdr, threshold)
{
  output_mat <- c()
  for (x in 1:ncol(data))
  {
    p_val <- c()
    for (i in 1:length(complex))
    {
      size_complex <- length(complex[[i]])

      if(size_complex < threshold)
      {
        p_val <- append(p_val, 1)
      }
      else
      {
        detected_prot <- rownames(data)[as.numeric(data[,x])!=0]
        prob_complex_exists <- (length(intersect(detected_prot, complex[[i]]))/size_complex) * (1 - fdr)
        prob_complex_exists_pval <- (1 - prob_complex_exists)

        p_val <- append(p_val, prob_complex_exists_pval)
      }
    }
    output_mat <- rbind(output_mat, p_val)
  }
  colnames(output_mat) <- names(complex)
  rownames(output_mat) <- colnames(data)
  return(output_mat)
}

#' @title This assigns probabilities to individual proteins based on PROTREC
#' @name PROTREC_protprob
#' @param cplx complex vector
#' @param p is the complex-based PROTREC probability, typically should be the 1 minus the result from certain sample of PROTREC_cplx_prob
#' @param prot is the list of observed proteins in the screen
#' @param fdr is the false discovery rate of the screen
#' @return The vector of PROTREC probabilities for the proteins of a sample
#' @export

PROTREC_protprob <- function(cplx, p, prot, fdr)

{

  new_mat <- c()
  for (i in 1:length(cplx))
  {
    protein_probabilities <- c()
    size_complex <- length(cplx[[i]])
    complex_names_vec <- rep(names(cplx)[i], size_complex)
    complex_prot_names_vec <- cplx[[i]]

    for(j in 1:length(cplx[[i]]))
    {
      if (cplx[[i]][j] %in% prot)
      {
        protein_probabilities <- append(protein_probabilities, (p[names(p) %in% names(cplx)[i]]) + (1 - fdr)*(1 - p[names(p) %in% names(cplx)[i]]))
      }
      else
      {
        protein_probabilities <- append(protein_probabilities, p[names(p) %in% names(cplx)[i]])
      }
    }

    temp_mat <- cbind(complex_prot_names_vec, complex_names_vec, round(protein_probabilities, 4))
    new_mat <- rbind(new_mat, temp_mat)
  }

  prot_list <- unique(sort(unlist(cplx)))
  output <- c()
  for (i in 1:length(prot_list))
  {
    output<- rbind(output,cbind(prot_list[i], max(unlist(new_mat[which(new_mat[,1]%in%prot_list[i]),3]))))
  }
  return(output)
}

###FCS functions########
#' @title This assigns probabilities to complexes based on FCS
#' @name fcs
#' @param data Proteomics expression data
#' @param complex_vector protein complex
#' @param sim_size how many iterations for simulation, recommend 1000
#' @param threshold the minimum size of the complex required
#' @return output_mat A matrix of FCS probabilities assigned to complexes for the dataset
#' @export

fcs <- function(data, complex_vector, sim_size, threshold)
{
  len=c()
  for(i in 1:length(complex_vector))
  {
    tp=length(complex_vector[[i]])
    if(!(tp %in% len) && tp>=5) len=append(len, tp)
  }
  testlis=list()
  for(o in 1:length(len))
  {
    test_mat=c()
    for(m in 1:sim_size)
    {
      set.seed(m+as.numeric(format(Sys.time(), "%S")))
      tmp=sample(unique(unlist(complex_vector)), len[o],replace=T)
      tmp=matrix(tmp,nrow=1)
      test_mat=rbind(test_mat,tmp)
    }
    testlis[[len[o]]]=test_mat
  }
  output_mat <- c()
  for (x in 1:ncol(data))
  {
    print(x)
    p_val <- c()
    for (i in 1:length(complex_vector))
    {
      size_complex <- length(complex_vector[[i]])

      if(size_complex < threshold)
      {
        p_val <- append(p_val, 1)
      }
      else
      {
        detected_prot <- rownames(data)[as.numeric(data[,x])!=0]
        obs_overlap <- length(intersect(detected_prot, complex_vector[[i]]))/size_complex
        test_intersects <- c()
        test_mt=testlis[[size_complex]]
        for (j in 1:nrow(test_mat))
        {
          test_intersects <- append(test_intersects, length(intersect(test_mt[j,], detected_prot))/ncol(test_mt))
        }

        p_val <- append(p_val,sum(as.numeric(test_intersects) >= as.numeric(obs_overlap))/length(test_intersects))
      }
    }
    output_mat <- rbind(output_mat, p_val)
  }
  colnames(output_mat) <- names(complex_vector)
  rownames(output_mat) <- colnames(data)
  return(output_mat)
}


#' @title This assigns probabilities to individual proteins based on FCS
#' @name fcs_prot_prob_assign
#' @param cplx complex vector
#' @param p is the complex-based FCS probability, typically should be the 1 minus the result from certain sample of fcs
#' @return The vector of FCS probabilities for the proteins of a sample
#' @export

fcs_prot_prob_assign <- function(cplx, p)
{
  complex_names_vec <- rep(names(cplx), lapply(cplx, length))
  data<- data.frame(table(complex_names_vec))
  data[,1]=as.numeric(as.character(data[,1]))
  data=data[order(data[,1]),]
  prot_vec <- unlist(cplx, use.names=F, recursive=F)
  new_mat <- as.matrix(cbind(prot_vec, complex_names_vec, rep(p, data[,2])))
  prot_list <- unique(prot_vec)
  output <- c()
  for (i in 1:length(prot_list))
  {
    output = rbind(output,cbind(prot_list[i], max(unlist(new_mat[which(new_mat[,1]%in%prot_list[i]), 3]))))
  }
  return(output)
}

###Hypergeometric Enrichment(HE) functions########
#' @title This assigns probabilities to complexes based on HE
#' @name hgtest
#' @param data Proteomics expression data
#' @param complex_vector protein complex
#' @param threshold the minimum size of the complex required
#' @return output_mat A matrix of HE probabilities assigned to complexes for the dataset
#' @export

hgtest <- function(data, complex_vector, threshold)
{
  output_mat <- c()
  for (x in 1:ncol(data))
  {
    print(x)
    p_val <- c()
    for (i in 1:length(complex_vector))
    {
      size_complex <- as.numeric(length(complex_vector[[i]]))
      if(size_complex < threshold)
      {
        p_val <- append(p_val, 1)
      }
      else
      {
        detected_prot <- rownames(data)[as.numeric(data[,x])!=0]
        obs_cplx_prot <- length(intersect(detected_prot, complex_vector[[i]]))
        total_prot <- as.numeric(length(union(detected_prot, unique(unlist(complex_vector, use.names=F, recursive=F)))))
        sample_size <- length(detected_prot) #k
        p_val <- append(p_val, phyper(q=obs_cplx_prot -1, m=size_complex, n=(total_prot- size_complex), k=sample_size, lower.tail=FALSE))
      }
    }
    output_mat <- rbind(output_mat, p_val)
  }

  colnames(output_mat) <- names(complex_vector)
  rownames(output_mat) <- colnames(data)
  return(output_mat)
}


#' @title This assigns probabilities to individual proteins based on HE
#' @name hgtest_prot_prob_assign
#' @param cplx complex vector
#' @param p is the complex-based HE probability, typically should be the 1 minus the result from certain sample of hgtest
#' @return output The vector of HE probabilities for the proteins of a sample
#' @export

hgtest_prot_prob_assign <- function(cplx, p)
{

  complex_names_vec <- rep(names(cplx), lapply(cplx, length))
  data<- data.frame(table(complex_names_vec))
  data[,1]=as.numeric(as.character(data[,1]))
  data=data[order(data[,1]),]
  prot_vec <- unlist(cplx, use.names=F, recursive=F)
  new_mat <- as.matrix(cbind(prot_vec, complex_names_vec, rep(p, data[,2])))
  prot_list <- unique(prot_vec)
  output <- c()
  for (i in 1:length(prot_list))
  {
    output<- rbind(output,cbind(prot_list[i], max(unlist(new_mat[which(new_mat[,1]%in%prot_list[i]),3]))))
  }
  return(output)
}


###Gene Set Enrichment Analysis(GSEA) functions########
#' @title This performs t-test inside GSEA function
#' @name mat.ttest
#' @param data Proteomics expression data
#' @return output_mat A matrix containing all p-values after t-test
#' @export
mat.ttest <- function(data)
{ output_mat<- c()
datat=as.data.frame(data)
for (i in 1:ncol(data))
{rep_pval<-c()
for (j in 1:nrow(data))
{ ttmp=as.matrix(datat[j,])
pval <- rowTtest(as.matrix(ttmp),y=NULL, mu=datat[j,i])$p.value
rep_pval<- append(rep_pval,pval)
}
output_mat<- cbind(output_mat,rep_pval)
rep_pval=c()

}

rownames(output_mat)<- rownames(data)
colnames(output_mat)<- colnames(data)
return(output_mat) }

#' @title This assigns probabilities to complexes based on GSEA
#' @name repgsea
#' @param data Proteomics expression data
#' @param complex_vector protein complex
#' @return gsea_pvals A matrix of GSEA probabilities assigned to complexes for the dataset
#' @export

repgsea <- function(data, complex_vector)
{
  ttest_out <- mat.ttest(data)
  gsea_pvals <- c()
  for (k in 1:ncol(data))
  {
    ranks<- rank(ttest_out[,k])
    rep_gseapval<-c()
    for (x in 1:length(complex_vector))
    {if (length(ranks[which(names(ranks) %in% complex_vector[[x]])]) >= 1)
    {ks_pval <- ks.test(jitter(ranks[which(names(ranks) %in% complex_vector[[x]])]), jitter(ranks[which(!names(ranks) %in% complex_vector[[x]])]))$p.value
    rep_gseapval <- append(rep_gseapval, ks_pval)
    }
      else
      {rep_gseapval <- append(rep_gseapval, 1)
      }}
    gsea_pvals<- rbind(gsea_pvals,rep_gseapval) }
  colnames(gsea_pvals)<-names(complex_vector)
  rownames(gsea_pvals)<-colnames(data)

  return(gsea_pvals) }


#' @title This assigns probabilities to individual proteins based on GSEA
#' @name gsea_prot_prob_assign
#' @param cplx complex vector
#' @param p is the complex-based GSEA probability, typically should be the 1 minus the result from certain sample of repgsea
#' @return output The vector of GSEA probabilities for the proteins of a sample
#' @export

gsea_prot_prob_assign <- function(cplx, p)
{

  complex_names_vec <- rep(names(cplx), lapply(cplx, length))
  data<- data.frame(table(complex_names_vec))
  data[,1]=as.numeric(as.character(data[,1]))
  data=data[order(data[,1]),]
  prot_vec <- unlist(cplx, use.names=F, recursive=F)
  new_mat <- as.matrix(cbind(prot_vec, complex_names_vec, rep(p, data[,2])))
  prot_list <- unique(prot_vec)
  output <- c()
  for (i in 1:length(prot_list))
  {
    output<- rbind(output,cbind(prot_list[i], max(unlist(new_mat[which(new_mat[,1]%in%prot_list[i]),3]))))
  }
  return(output)
}

###Recovery Rate########

#for checking recovery and providing a significance value for the recovery for FCS,HE and GSEA
#' @title This works out the recovery significance for FCS,HE and GSEA
#' @name pairwise_recovery
#' @param predict_list A vector of significant proteins for a given sample, default set 0.05 p-value as significant
#' @param original_prot_list The original set of proteins observed for that particular sample
#' @param check_prot_list A set of proteins observed in a second replicate
#' @param complex_vec A list of complex objects
#' @return output A vector of five pieces of information: the observed proporetion of overlap,
#' the significance of this overlap (p-value), the number of verified proteins, the total number of predicted missing proteins,
#' and the list of validated proteins separated by 'a'.
#' @export

pairwise_recovery <- function(predict_list, original_prot_list, check_prot_list, complex_vec)
{

  predict_list <- predict_list[as.numeric(predict_list) <= 0.05]

  add_proteins <- setdiff(unlist(complex_vec[names(predict_list)]), original_prot_list)
  obs_intersect <- round(length(intersect(add_proteins, check_prot_list))/length( add_proteins), 3)
  if(length(add_proteins)==0) obs_intersect=0
  theoretical_vec <- c()

  for (i in 1:1000)
  {
    S_prime <- sample(unlist(complex_vec), length(add_proteins))
    add_proteins_prime <- setdiff(S_prime, original_prot_list)
    theoretical_intersect <- round(length(intersect(add_proteins_prime, check_prot_list))/length(add_proteins_prime), 3)
    if(!is.na(theoretical_intersect))
      theoretical_vec <- append(theoretical_vec, theoretical_intersect)
  }

  p_value <- sum(theoretical_vec >= obs_intersect)/length(theoretical_vec)
  a=length(add_proteins)
  b=length(intersect(add_proteins, check_prot_list))
  if(obs_intersect==0)
  {
    p_value=0
    a=0
    b=0
  }
  return(c(obs_intersect, p_value, a, b, add_proteins,"a",intersect(add_proteins, check_prot_list)))
}


#for checking recovery and providing a significance value for the recovery for PROTREC
#' @title This works out the recovery significance for PROTREC
#' @name pairwise_recovery_protrec
#' @param prot_predict_list A vector of significant proteins for a given sample
#' @param original_prot_list The original set of proteins observed for that particular sample
#' @param check_prot_list A set of proteins observed in a second replicate
#' @param complex_vec A list of complex objects
#' @param protrecscoreset The threshold of PROTREC score, default 0.95
#' @return output A vector of five pieces of information: the observed proporetion of overlap,
#' the significance of this overlap (p-value), the number of verified proteins, the total number of predicted missing proteins,
#' and the list of validated proteins separated by 'a'.
#' @export

pairwise_recovery_protrec <- function(prot_predict_list, original_prot_list, check_prot_list, complex_vec,protrecscoreset=0.95)
{

  prot_predict_list <- prot_predict_list[which(as.numeric(prot_predict_list[,2]) >= protrecscoreset),]
  add_proteins <- setdiff(prot_predict_list[,1], original_prot_list)
  obs_intersect <- round(length(intersect(add_proteins, check_prot_list))/length( add_proteins), 3)

  theoretical_vec <- c()

  for (i in 1:1000)
  {
    S_prime <- sample(unlist(complex_vec), length(add_proteins))
    add_proteins_prime <- setdiff(S_prime, original_prot_list)
    theoretical_intersect <- round(length(intersect(add_proteins_prime, check_prot_list))/length(add_proteins_prime), 3)
    theoretical_vec <- append(theoretical_vec, theoretical_intersect)
  }

  p_value <- sum(as.numeric(theoretical_vec >= obs_intersect))/length(theoretical_vec)

  return(c(obs_intersect, p_value, length(add_proteins), length(intersect(add_proteins, check_prot_list)), add_proteins,"a",intersect(add_proteins, check_prot_list)))
}







