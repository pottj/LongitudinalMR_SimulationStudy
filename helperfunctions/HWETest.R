#' @title Test for Hardy-Weinberg Equilibrium
#' @description This function tests the genotype distribution of one bi-allelic SNP in a population for deviation from Hardy-Weinberg equilibrium (HWE).  
#' @param AA number of individuals with genotype AA 
#' @param AB number of individuals with genotype AB 
#' @param BB number of individuals with genotype BB
#' @return A data.table with 
#' * observed allele and genotype frequencies
#' * expected genotype frequencies under HWE
#' * $\chi^2$ and p-value for 1 degree of freedom
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  HWETest(625, 3750, 5625)
#'  }
#' }
#' @rdname HWETest
#' @export 
HWETest = function(AA,AB,BB){
  # AA = dim(myDistTab[genotype==0,])[1]
  # AB = dim(myDistTab[genotype==1,])[1]
  # BB = dim(myDistTab[genotype==2,])[1]
  
  # step 1: get allele frequency
  N = AA + AB+ BB
  p = (2*AA + AB)/(2*N)
  q = (2*BB + AB)/(2*N)
  
  # step 2: get expected genotype frequency
  Exp<-c()
  Exp[1]<-p^2
  Exp[2]<-2*p*q
  Exp[3]<-q^2
  
  # step 3: get observed genotzpe frequency
  Obs<-c()
  Obs[1]<-AA/N
  Obs[2]<-AB/N
  Obs[3]<-BB/N
  
  # step 4: calculate chi^2
  x<-N*(  ((Obs[1]-Exp[1])^2/Exp[1]) + 
            ((Obs[2]-Exp[2])^2/Exp[2]) + 
            ((Obs[3]-Exp[3])^2/Exp[3]))
  pval = c(1-pchisq(x,df=1))
  
  # step 5: return results
  res = data.table::data.table(AF_A = p,
                               AF_B = q,
                               AA_obs = Obs[1],
                               AB_obs = Obs[2],
                               BB_obs = Obs[3],
                               AA_exp = Exp[1],
                               AB_exp = Exp[2],
                               BB_exp = Exp[3],
                               CHISQ = x,
                               Pvalue = pval)
  return(res)
  
}
