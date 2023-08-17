#' @title MVMR function to get all causal estimates over several outcomes
#' @description In this function, I test for each outcome an MVMR model (classical = as defined in the MendelianRandomization package, correction = using mrest_me function to correct for possible measurement error)
#' @param data_GX_long data.table with SNP associations with any exposure levels 
#' @param data_GY_long data.table with SNP associations with any outcome types 
#' @param filterBadSNPs A Boolean, default is T. Should bad SNPs (weak instruments) be filtered for the analysis?  
#' @param filterBadSNPs_threshold threshold value to filter bad SNPs in case of filterBadSNPs == T, default: 1e-06
#' @return data.table with columns for exposure, outcome, beta_IVW, SE_IVW, pval_IVW (1 = classic, 2 = correction)
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname MVMR_jp
#' @export
MVMR_jp = function(data_GX_long, data_GY_long, filterBadSNPs = T, filterBadSNPs_threshold = 1e-6){
  #debug
  # data_GX_long = copy(myTab_SNPAssocs_X)
  # data_GY_long = copy(myTab_SNPAssocs_Y)
  # filterBadSNPs = MR_filterBadSNPs
  # filterBadSNPs_threshold = 1e-6
  
  data_GX_wide = dcast(data_GX_long, SNP ~ exposure, value.var = c("beta","SE","tval","pval"))
  outcomes = unique(data_GY_long$outcome)
  exposures = unique(data_GX_long$exposure)
  NR_exposures = length(exposures)
  
  if(filterBadSNPs == T){
    filt = data_GX_long[, pval<filterBadSNPs_threshold]
    goodSNPs = data_GX_long[filt,SNP]
    goodSNPs = unique(goodSNPs)
    data_GX_wide = data_GX_wide[SNP %in% goodSNPs,]
    data_GY_long = data_GY_long[SNP %in% goodSNPs,]
    
  }else{
    goodSNPs = data_GX_long[,SNP]
    goodSNPs = unique(goodSNPs)
  }
  
  stopifnot("There are not enough SNPs left to do an MVMR!" = length(goodSNPs)>NR_exposures)
  
  dumTab_MR = foreach(i = 1:length(outcomes))%do%{
    #i=1
    myOutcome = outcomes[i]
    myTab_GY = copy(data_GY_long)
    myTab_GY = myTab_GY[outcome == myOutcome,]
    
    filt1 = grepl("beta",names(data_GX_wide))
    data_beta = copy(data_GX_wide)
    data_beta = data_beta[,filt1,with=F]
    filt2 = grepl("SE",names(data_GX_wide))
    data_SE = copy(data_GX_wide)
    data_SE = data_SE[,filt2,with=F]
    
    mvmr_obj = mr_mvinput(bx = as.matrix(data_beta),
                          bxse = as.matrix(data_SE),
                          by = myTab_GY$beta, 
                          byse = myTab_GY$SE,
                          exposure = exposures,
                          outcome = myOutcome)
    
    res3 = mrest_me(mvmr_obj)
    res2 = mr_mvivw(mvmr_obj)    
    
    res = data.table(exposure = c(res2@Exposure),
                     outcome = rep(res2@Outcome,NR_exposures),
                     beta_IVW1 = c(res2@Estimate),
                     SE_IVW1 = c(res2@StdError),
                     pval_IVW1 = c(res2@Pvalue),
                     HeteroStat = rep(res2@Heter.Stat[1],NR_exposures),
                     HeteroStat_pval = rep(res2@Heter.Stat[2],NR_exposures), 
                     beta_IVW2 = c(res3$thest),
                     SE_IVW2 = c(sqrt(res3$Var[1,1]),sqrt(res3$Var[2,2])))
    res[, pval_IVW2 := pnorm(-abs(beta_IVW2/SE_IVW2))*2]
    res
  }
  myTab_MR = rbindlist(dumTab_MR)
  
  return(myTab_MR)
}
