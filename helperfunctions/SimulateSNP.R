#' @title Simulate genotype data
#' @description This functions simulates genotypes in Hardy-Weinberg equilibrium given a certain minor allele frequency
#' @param n_samples Number of individuals for which the genotypes should be simulated
#' @param p Desired allele frequency of the SNP
#' @param SNP_centered A Boolean, default is F. Should the genotypes be centered to 0 (TRUE) or not?
#' @return A numeric vector containing the genotypes coded as 0, 1, 2 for AA, AB, and BB. p(A) will roughly correspond to input p, not necessarily the minor allele. 
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  table(SimulateSNP(100,0.25,F))
#'  }
#' }
#' @rdname SimulateSNP
#' @export 
SimulateSNP = function(n_samples, p, SNP_centered=F){
  # debug
  # n_samples = 10000
  # p = 0.25
  # SNP_centered=F
  
  # check input
  #stopifnot(p<=0.5)
  
  # HWE distribution
  q=1-p
  p_AA = p^2
  p_AB = 2*p*q
  p_BB = q^2
  
  G = sample(x = c(0,1,2), 
             size = n_samples, 
             replace = T, 
             prob = c(p_AA, p_AB, p_BB))
  
  if(SNP_centered==T){
    table2 = table(G)
    p_obs = (2*table2[1] + table2[2])/n_samples
    q_obs = (2*table2[3] + table2[2])/n_samples
    G =  G - (2*q_obs)
    # G =  G/sqrt(2*p_obs*q_obs)
  }
  
  return(G)
  
}
