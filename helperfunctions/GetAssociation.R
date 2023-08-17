#' @title Get SNP associations with exposure or outcome
#' @description This function estimates the beta coefficients and standard errors for all SNP-phenotype association necessary in my simulation setup
#' @param data data.table with phenotype column
#' @param method parameter indicating what type of phenotype is given in the phenotype column and which regression model should be used. Must be either "linReg", "linMixed", "gamlss".
#' @param genotypes Genotype matrix
#' @param dep_var_name name of phenotype column 
#' @return data.table containing one row per SNP from the genotype matrix (in case of "linReg") or two rows per SNP (in case of "linMixed": intercept and slope; in case of "gamlss": mu and sigma). Columns are SNP number, exposure (as given in dep_var_name), beta, SE, tval, and pval. 
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetAssociation
#' @export
GetAssociation = function(data,method,genotypes,dep_var_name){
  # debug
  # genotypes = G
  # method = "linMixed"
  # data = copy(TestData)
  # dep_var_name = "myX"
  
  stopifnot(method %in% c("linReg","linMixed","gamlss"))
  
  # step 1: get number of SNPs to be tested
  SNPs_NR = dim(G)[2]
  
  # step 2: loop per SNP
  if(method == "linReg"){
    
    modTab = foreach(j = 1:SNPs_NR)%do%{
      #j=1
      Gj = G[,j]
      
      data2 = copy(data)
      data2[, mySNP := Gj]
      data2[, myVar := get(dep_var_name)]
      
      mod1 = lm(myVar ~ mySNP, data = data2)
      #summary(mod1)
      data2[, mySNP := NULL]
      data2[, myVar := NULL]
      
      res = data.table(SNP = j,
                       exposure = dep_var_name,
                       beta = summary(mod1)$coef[2,1],
                       SE = summary(mod1)$coef[2,2],
                       tval = summary(mod1)$coef[2,3],
                       pval = summary(mod1)$coef[2,4])
      res
    }
    
  }else if(method == "linMixed"){
    time_NR = data[,.N,by=ID]
    stopifnot(length(unique(time_NR$N))==1)
    modTab = foreach(j = 1:SNPs_NR)%do%{
      #j=1
      Gj = G[,j]
      data2 = copy(data)
      data2[,mySNP := rep(Gj, each=time_NR$N[1])]
      data2[,myVar := get(dep_var_name)]
      
      modX_G = lmer(myVar ~ mySNP + time + (1|ID) + mySNP:time, data = data2)
      summary(modX_G)
      
      res = data.table(SNP = c(j,j),
                       exposure = c("intercept","slope"),
                       beta = c(summary(modX_G)$coef[2,1],summary(modX_G)$coef[4,1]),
                       SE = c(summary(modX_G)$coef[2,2],summary(modX_G)$coef[4,2]),
                       tval = c(summary(modX_G)$coef[2,3],summary(modX_G)$coef[4,3]))
      res[,pval := pnorm(-abs(beta/SE))*2]
      res
    }
    
  }else if(method == "gamlss"){
    time_NR = unique(data$time)
    modTab = foreach(j = 1:SNPs_NR)%do%{
      #j=1
      Gj = G[,j]
      data2 = copy(data)
      data2[,mySNP := rep(Gj, each=length(time_NR))]
      data2[,myVar := get(dep_var_name)]
      
      modX_G = gamlss(myVar ~ mySNP + scale(time), random=~1|ID, sigma.formula = ~mySNP, 
                      data = data2, family = "NO")
      dummy = summary(modX_G)
      
      res = data.table(SNP = c(j,j),
                       exposure = c("mu","sigma"),
                       beta = c(dummy[2,1],dummy[5,1]),
                       SE = c(dummy[2,2],dummy[5,2]),
                       tval = c(dummy[2,3],dummy[5,3]),
                       pval = c(dummy[2,4],dummy[5,4]))
      res
    }
  }
  myAssocs = rbindlist(modTab)
  
  # step 3: return results
  return(myAssocs)
  
}
