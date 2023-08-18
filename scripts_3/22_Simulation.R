#' ---
#' title: "Simulation Study - Scenario 22"
#' subtitle: "Shared SNPs - all time points"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output: github_document
#' ---
#'
#' # Initialize ####
#' ***

rm(list = ls())
time0<-Sys.time()

source("../settings_3/22_ParameterSettings.R")

set.seed(2023)

if(dir.exists(paste0("../",outfiles_dir,"/",outfiles_prefix))==F){
  dir.create(paste0("../",outfiles_dir,"/",outfiles_prefix))
  message("Created results folder ",paste0("../",outfiles_dir,"/",outfiles_prefix))
}else{
  message("Using pre-existing results folder ",paste0("../",outfiles_dir,"/",outfiles_prefix))
}

message("Fixed parameters: 
        \n - number of samples: ",n_samples,
        "\n - number of time points: ",n_times,
        "\n - number of simulations: ",n_sim,
        "\n - number of SNPs: ",SNPs_NR,
        "\n - number of cores: ",n_cores)

message("Global setting: 
        \n - number of SNP sets: ",length(unique(SNPs_classes))," (1 = shared set; 2=distinct sets for exposure 1 and 2)         \n - timepoints per individual: ",n_times_random, " (FALSE = no missing timepoints; TRUE = random missing timepoints per individual)")

message("Scenario settings: 
        \n - second function: ",X_useAs2ndFunction," (1 = fluctuation using sinus (amplitude); 2 = trend-like using sinus (wavelength); 3 = trend using linear slope) \n - exposure model: ",AssocModel,
        "\n - correction ov MVMR estiamtes: ",MR_doCorrection," (FALSE = no correction; TRUE = correct with beta coefficient of linear regression of allele scores on exposures; only relevant for model meanSD and eigenfunc)")

#' # Simulation ####
#' ***

registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
#registerDoParallel(n_cores)

counter = seq(1,n_sim,n_sim/10)
SimTab = foreach(s = 1:n_sim)%dorng%{
  #s=1
  #message("Working on simulation ",s)
  source("../settings_3/22_ParameterSettings.R")
  
  #' ## Step 0: create directory to store data
  #' ***
  {
    outdir_sim = paste0("../",outfiles_dir,"/",outfiles_prefix,"/Simulation_",s,"/")
    if(save_data_perSim == T | s %in% counter){
      if(dir.exists(outdir_sim)==F){
        dir.create(outdir_sim)
        message("Created results folder ",outdir_sim, " for simulation ",s)
      }else{
        message("Using pre-existing results folder ",outdir_sim, " for simulation ",s)
      }
    }
  }
  
  #' ## Step 1: Simulate genotypes G ####
  #' ***
  #' Here, parameter *SNPs_classes* defines shared or distinct SNP sets for mean and variability
  {
    #' Get genotypes 
    G = matrix(data = NA, nrow = n_samples, ncol=SNPs_NR)
    
    HWETab = foreach(j = 1:SNPs_NR)%do%{
      #j=1
      
      # simulate SNP
      SNP = SimulateSNP(n_samples = n_samples,
                        p = SNPs_EAF,
                        SNP_centered = SNPs_centering)
      G[,j] = SNP
      
      table1 = table(SNP)
      names(table1)[names(table1)==0] = "AA"
      names(table1)[names(table1)==1] = "AB"
      names(table1)[names(table1)==2] = "BB"
      
      AA = table1[names(table1)=="AA"]
      AB = table1[names(table1)=="AB"]
      BB = table1[names(table1)=="BB"]
      
      if(length(AA)==0) AA = 0
      if(length(AB)==0) AB = 0
      if(length(BB)==0) BB = 0
      
      # test HWE
      test = HWETest(AA = AA, 
                     AB = AB,
                     BB = BB)
      test[,NR := j]
      test
    }
    HWETab = rbindlist(HWETab)
    
    #' Check for HWE violation
    HWETab[Pvalue <= 0.05/SNPs_NR,]
    HWETab[Pvalue <= 0.05,]
    if(save_data_perSim == T | s %in% counter) save(HWETab, file = paste0(outdir_sim, "/01_HWE.RData"))
    
    #' Check for LD violation
    CorTab = cor(G)^2
    filt = CorTab>= 0.1
    stopifnot(sum(filt) == SNPs_NR)
    stopifnot(sum(diag(CorTab)) == SNPs_NR)
    if(save_data_perSim == T | s %in% counter) save(G, file = paste0(outdir_sim, "/01_GenotypeMatrix.RData"))
    
    #' Get SNP effect 
    if(length(unique(SNPs_classes))==2){
      X_beta_func1 = rnorm(sum(SNPs_classes=="A"), X_beta_mean[1], X_beta_sd[1])
      X_beta_func2 = rnorm(sum(SNPs_classes=="B"), X_beta_mean[2], X_beta_sd[2])
      filtA = SNPs_classes == "A"
      filtB = SNPs_classes == "B"
    }else{
      X_beta_func1 = rnorm(SNPs_NR, X_beta_mean[1], X_beta_sd[1])
      X_beta_func2 = rnorm(SNPs_NR, X_beta_mean[2], X_beta_sd[2])
      filtA = SNPs_classes == "A"
      filtB = SNPs_classes == "A"
    }
    AS1 = G[,filtA] %*% X_beta_func1
    if(SNPs_Correct_ASs==T) AS1 = (AS1 - mean(AS1))/sd(AS1)
    AS2 = exp(-0.5*G[,filtB]%*%X_beta_func2)
    if(SNPs_Correct_ASs==T) AS2 = AS2/0.0375
    
  }
  
  #' ## Step 2: Simulate exposure X and GX association ####
  #' ***
  #' Here, parameter *X_useAs2ndFunction* defines the second function used to simulate the time-dependent exposure
  #' 
  #' Here, parameter *AssocModel* defines which association model and exposure type should be used
  {
    #' Get random effects
    CoVarMatrix = diag(x=X_var_random,nrow = 2)
    REff = MASS::mvrnorm(n_samples, 
                         mu=X_mean_random, 
                         Sigma=(CoVarMatrix)^2)
    
    #' Get exposure
    myTabX_long= data.table(ID = rep(1:n_samples, each=n_times),
                            time = rep(1:n_times, times=n_samples),
                            random_intercept = rep(REff[,1], each=n_times),
                            random_trend = rep(REff[,2], each=n_times),
                            genetic_const=rep(AS1, each=n_times),
                            genetic_trend=rep(AS2, each=n_times))
    
    myTabX_long[,time := (time-1)*3]
    func1 = get(X_useAs1stFunction)
    func2 = get(X_useAs2ndFunction)
    myTabX_long[,genetic_func1 := func1(t = time, g = genetic_const, u = random_intercept)]
    myTabX_long[,genetic_func2 := func2(t = time, g = genetic_trend, u = random_trend)]
    
    X = with(myTabX_long, (genetic_func1 + genetic_func2) 
             + rnorm(n=n_samples*n_times, mean=X_error_mean, sd=sqrt(X_error_sd)))
    myTabX_long[,X := X]
    
    myTabX_wide = dcast(myTabX_long, ID + genetic_const + genetic_trend + random_intercept + random_trend ~ time, value.var = "X")
    names(myTabX_wide)[6:(5+n_times)] = paste0("X", names(myTabX_wide)[6:(5+n_times)])
    myTabX_reduced = cbind(myTabX_wide[,1:5])
    
    #' Get some plots
    if((do_plotting == T & save_data_perSim==T) | s %in% counter){
      plot1 = ggplot()+
        geom_line(aes(y=X, 
                      x=time, 
                      group=ID, 
                      colour=genetic_const), 
                  data=myTabX_long, 
                  show.legend = TRUE) + 
        labs(color="PGS for X") +
        theme(legend.position = "none") + theme_classic()
      
      tiff(filename = paste0(outdir_sim, "02_TimeVsExposure_PGS1.tiff"), 
           width = 900, height = 600, res=125, compression = 'lzw')
      print(plot1)
      dev.off()
      
      plot2 = ggplot()+
        geom_line(aes(y=X, 
                      x=time, 
                      group=ID, 
                      colour=genetic_trend), 
                  data=myTabX_long, 
                  show.legend = TRUE) + 
        labs(color="PGS for X") +
        theme(legend.position = "none") + theme_classic()
      
      tiff(filename = paste0(outdir_sim, "02_TimeVsExposure_PGS2.tiff"), 
           width = 900, height = 600, res=125, compression = 'lzw')
      print(plot2)
      dev.off()
    }
    
    if(n_times_random == T){
      #' Filter for random time points
      myTabX_long[,time_alt := time/3]
      uniqueTimepoints = myTabX_long[,unique(time_alt)]
      myTabX_reduced[,NR_timePoints := runif(n=n_samples,min = n_times_random_min,max = n_times_random_max)]
      myTabX_reduced[,NR_timePoints := round(NR_timePoints,0)]
      dumTab = foreach(i = 1:n_samples)%do%{
        #i=2
        start = (i-1)*30 + 1
        end = i*30
        test = myTabX_long[start:end,]
        filt_timepoints = sample(uniqueTimepoints,myTabX_reduced[i,NR_timePoints],replace = F)
        test = test[time_alt %in% filt_timepoints,]
        test
      }
      myTabX_long = rbindlist(dumTab)
      
      #' Get some plots
      if((do_plotting == T & save_data_perSim==T) | s %in% counter){
        plot1 = ggplot()+
          geom_line(aes(y=X,
                        x=time,
                        group=ID,
                        colour=genetic_const),
                    data=myTabX_long,
                    show.legend = TRUE) +
          labs(color="PGS for X") +
          theme(legend.position = "none") + theme_classic()
        
        tiff(filename = paste0(outdir_sim, "02_TimeVsExposure_PGS1_filtered.tiff"),
             width = 900, height = 600, res=125, compression = 'lzw')
        print(plot1)
        dev.off()
        
        plot2 = ggplot()+
          geom_line(aes(y=X,
                        x=time,
                        group=ID,
                        colour=genetic_trend),
                    data=myTabX_long,
                    show.legend = TRUE) +
          labs(color="PGS for X") +
          theme(legend.position = "none") + theme_classic()
        
        tiff(filename = paste0(outdir_sim, "02_TimeVsExposure_PGS2_filtered.tiff"),
             width = 900, height = 600, res=125, compression = 'lzw')
        print(plot2)
        dev.off()
      }
      
    }
    
    #' Do dimension reduction
    if(AssocModel == "eigenfunc"){
      fPCA_input = MakeFPCAInputs(IDs = myTabX_long[,ID], 
                                  tVec=myTabX_long[,time], 
                                  t(myTabX_long[,X]))
      fPCA_dense = FPCA(fPCA_input$Ly, fPCA_input$Lt)
      if(do_plotting == T){
        plot(fPCA_dense)
      }
      fPCA_dense$lambda/sum(fPCA_dense$lambda)
      
      # select relevant eigenfunctions only e.g. lambda/sum(lambdas)>5%
      filt1 = fPCA_dense$lambda/sum(fPCA_dense$lambda)
      filt2 = filt1>eigenfunc_cutoff
      myEigenfunctions = fPCA_dense$xiEst
      myEigenfunctions = myEigenfunctions[,filt2]
      NR_cols_reduced = dim(myTabX_reduced)[2]
      myTabX_reduced = cbind(myTabX_reduced,myEigenfunctions)
      NR_eigenfunctions = dim(myTabX_reduced)[2]-NR_cols_reduced
      names(myTabX_reduced)[(NR_cols_reduced+1):(NR_cols_reduced+NR_eigenfunctions)] = paste0("eigfunc",c(1:(NR_eigenfunctions)))
      
      factor_1 = summary(lm(eigfunc1 ~ genetic_const + genetic_trend, data=myTabX_reduced))$coeff[2,1]
      if(NR_eigenfunctions>=2) factor_2 = summary(lm(eigfunc2 ~ genetic_const + genetic_trend, data=myTabX_reduced))$coeff[3,1]
      if(NR_eigenfunctions>=3) factor_3 = summary(lm(eigfunc3 ~ genetic_const + genetic_trend, data=myTabX_reduced))$coeff[3,1]
      
    }else if(AssocModel == "meanSD"){
      mean = myTabX_long[,mean(X),by = ID]
      sd = myTabX_long[,sd(X),by = ID]
      myTabX_reduced[,mean := mean$V1]
      myTabX_reduced[,sd := sd$V1]
      myTabX_reduced[,sd := log(sd)]
      
      factor_1 = summary(lm(mean ~ genetic_const + genetic_trend, data=myTabX_reduced))$coeff[2,1]
      factor_2 = summary(lm(sd ~ genetic_const + genetic_trend, data=myTabX_reduced))$coeff[3,1]
    }
    
    #' Save data
    if(save_data_perSim == T | s %in% counter){
      save(myTabX_wide, myTabX_long, myTabX_reduced,
           file = paste0(outdir_sim, "02_ExposureData.RData"))
    }
    
    #' Get SNP association with exposure
    if(AssocModel %in% c("eigenfunc","meanSD")){
      if(n_times_random == T){
        myVars = names(myTabX_reduced)[-c(1:6)]
      }else{
        myVars = names(myTabX_reduced)[-c(1:5)]
      }
      myTab_SNPAssocs = foreach(i = 1:length(myVars))%do%{
        #i=1
        myVar = myVars[i]
        TestData = copy(myTabX_reduced)
        TestData[,myX := get(myVar)]
        #TestData[,time := 0]
        
        myTab_GX1 = GetAssociation(data = TestData, 
                                   method = "linReg", 
                                   genotypes = G, 
                                   dep_var_name = "myX")
        
        myTab_GX1[,exposure := myVar]
        myTab_GX1
      }
      myTab_SNPAssocs = rbindlist(myTab_SNPAssocs)
    }else if(AssocModel == "linMixed"){
      if(linMixed_random == T){
        dumTab = foreach(i = 1:n_samples)%do%{
          #i=2
          test = myTabX_long[ID == i,]
          uniqueTimepoints_i = test$time
          filt_timepoints = sample(uniqueTimepoints_i,linMixed_NRtimepoints,replace = F)
          test = test[time %in% filt_timepoints,]
          test
        }
        TestData = rbindlist(dumTab)
      }else{
        dumTab = foreach(i = 1:n_samples)%do%{
          #i=2
          test = myTabX_long[ID == i,]
          test = test[c(1:linMixed_NRtimepoints),]
          test
        }
        TestData = rbindlist(dumTab)
      }
      
      myTab_SNPAssocs = GetAssociation(data = TestData, 
                                       method = "linMixed", 
                                       genotypes = G, 
                                       dep_var_name = "X")
      
    }else{
      myTab_SNPAssocs = GetAssociation(data = myTabX_long, 
                                       method = "gamlss", 
                                       genotypes = G, 
                                       dep_var_name = "X")
    }
    myTab_SNPAssocs_X = copy(myTab_SNPAssocs)
    
    #' Save data
    if(save_data_perSim == T | s %in% counter){
      save(myTab_SNPAssocs_X,
           file = paste0(outdir_sim, "02_AssociationData.RData"))
    }
  }
  
  #' ## Step 3: Simulate outcome Y and GY association ####
  #' ***
  #' Here, parameter *Y_AS_new* defines if new allele scores should be calculated using the betas from the GX association table
  {
    if(Y_AS_new == T){
      X_beta_func1_new = myTab_SNPAssocs_X[exposure %in% c("mean","eigfunc1","intercept","mu"),beta]
      AS1_new = G[,filtA] %*% X_beta_func1_new[filtA]
      X_beta_func2_new = myTab_SNPAssocs_X[exposure %in% c("sd","eigfunc2","slope","sigma"),beta]
      if(length(X_beta_func2_new) == 0){
        AS2_new = AS2
      }else{
        AS2_new = G[,filtB] %*% X_beta_func2_new[filtB]  
      }
    }else{
      AS1_new = AS1
      AS2_new = AS2
    }
    
    myTabX_reduced[,Y1 := rnorm(n=n_samples,mean = Y_mean_random,sd = Y_var_random)]
    myTabX_reduced[,Y2 := Y_alpha[1] * AS1_new + rnorm(n=n_samples,mean = Y_mean_random,sd = Y_var_random)]
    myTabX_reduced[,Y3 := Y_alpha[2] * AS2_new + rnorm(n=n_samples,mean = Y_mean_random,sd = Y_var_random)]
    myTabX_reduced[,Y4 := Y_alpha[1] * AS1_new + Y_alpha[2] * AS2_new + 
                     rnorm(n=n_samples,mean = Y_mean_random,sd = Y_var_random)]
    
    #' Save data
    if(save_data_perSim == T | s %in% counter){
      save(myTabX_reduced,
           file = paste0(outdir_sim, "03_OutcomeData.RData"))
    }
    
    #' Get associations
    myVars = c("Y1","Y2","Y3","Y4")
    
    myTab_SNPAssocs = foreach(i = 1:length(myVars))%do%{
      #i=1
      myVar = myVars[i]
      TestData = copy(myTabX_reduced)
      TestData[,myX := get(myVar)]
      
      myTab_GY1 = GetAssociation(data = TestData, 
                                 method = "linReg", 
                                 genotypes = G, 
                                 dep_var_name = "myX")
      
      setnames(myTab_GY1,"exposure","outcome")
      myTab_GY1[,outcome := myVar]
      myTab_GY1
    }
    myTab_SNPAssocs_Y = rbindlist(myTab_SNPAssocs)
    
    #' Save data
    if(save_data_perSim == T | s %in% counter){
      save(myTab_SNPAssocs_Y,
           file = paste0(outdir_sim, "03_AssociationData.RData"))
    }
  }
  
  #' ## Step 4: Run MVMR ####
  #' ***
  #' Here, parameter *MR_doCorrection* defines if the MVMR estimates for mean, sd, eigfunc1 and eigfunc2 should be corrected or not
  {
    #' Run analysis
    if(!exists("NR_eigenfunctions")) NR_eigenfunctions = 2
    if(AssocModel != "eigenfunc" | (AssocModel == "eigenfunc" & NR_eigenfunctions>=2)){
      MR_analysis = MVMR_jp(data_GX_long = myTab_SNPAssocs_X,
                            data_GY_long = myTab_SNPAssocs_Y, 
                            filterBadSNPs = MR_filterBadSNPs,
                            filterBadSNPs_threshold = MR_filterBadSNPs_treshold)
      
      MRTab = copy(MR_analysis)
      
      if(AssocModel %in% c("meanSD","eigenfunc") & MR_doCorrection == T){
        MRTab[exposure %in% c("sd","eigfunc2"),beta_IVW2 := beta_IVW2*factor_2]
        MRTab[exposure %in% c("sd","eigfunc2"),beta_IVW1 := beta_IVW1*factor_2]
        MRTab[exposure %in% c("mean","eigfunc1"),beta_IVW2 := beta_IVW2*factor_1]
        MRTab[exposure %in% c("mean","eigfunc1"),beta_IVW1 := beta_IVW1*factor_1]
        if(NR_eigenfunctions>2){
          MRTab[exposure %in% c("eigfunc3"),beta_IVW1 := beta_IVW1*factor_3]
          MRTab[exposure %in% c("eigfunc3"),beta_IVW2 := beta_IVW2*factor_3]
        } 
      }
      
      #' Get some plots
      if((do_plotting == T & save_data_perSim==T) | s %in% counter){
        plot3 = ggplot(MRTab, aes(x=outcome, y=beta_IVW1 )) + 
          facet_wrap(~exposure,scales = "free")+
          #geom_line() +
          geom_point()+
          geom_errorbar(aes(ymin=beta_IVW1 -1.96*SE_IVW1, ymax=beta_IVW1 +1.96*SE_IVW1), width=.2)
        tiff(filename = paste0(outdir_sim, "04_MREstimates.tiff"), 
             width = 900, height = 600, res=125, compression = 'lzw')
        print(plot3)
        dev.off()
      }
      
    }else{
      MR_analysis = MR_jp(data_GX_long = myTab_SNPAssocs_X,
                          data_GY_long = myTab_SNPAssocs_Y, 
                          filterBadSNPs = MR_filterBadSNPs,
                          filterBadSNPs_threshold = MR_filterBadSNPs_treshold)
      MRTab = copy(MR_analysis)
      
      if(AssocModel %in% c("meanSD","eigenfunc") & MR_doCorrection == T){
        MRTab[exposure %in% c("mean","eigfunc1"),beta_IVW2 := beta_IVW2*factor_1]
        MRTab[exposure %in% c("mean","eigfunc1"),beta_IVW1 := beta_IVW1*factor_1]
      }
    }
    
    
    #' Save data
    if(save_data_perSim == T | s %in% counter){
      save(MRTab,
           file = paste0(outdir_sim, "04_MVMRResults.RData"))
    }
    
  }
  
  MRTab[,n_sim := s]
  MRTab
}
SimTab = rbindlist(SimTab,fill=T)
save(SimTab, file = paste0("../",outfiles_dir,"/",outfiles_prefix, "_SimAll.RData"))

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
