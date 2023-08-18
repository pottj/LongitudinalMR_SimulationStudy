#############################
# this is a template source file
# please change all parameters accordingly
#############################

#############################
# R library and R packages
#############################
.libPaths()

suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(doRNG))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(gamlss))
suppressPackageStartupMessages(library(fdapace))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(MendelianRandomization))

#############################
# Other helper functions 
#############################
# Functions to simulate the SNP data
source("../helperfunctions/SimulateSNP.R")
source("../helperfunctions/HWETest.R")

# Functions to get SNP associations
source("../helperfunctions/GetAssociation.R")

# MVMR function (helperfunctions in addition to MendelianRandomization)
source("../helperfunctions/MVMR_jp.R")
source("../helperfunctions/mrest_me.R")
source("../helperfunctions/MR_jp.R")

#############################
# General parameters simulation
#############################
n_samples = 10000
n_times = 30
n_times_random = F
n_times_random_min = 20 
n_times_random_max = n_times
n_sim = 1000
outfiles_prefix = "22_Sim_trend_meanSD_MVMRcor"
outfiles_dir = "results_2"
save_data_perSim = F
do_plotting = F
n_cores = 10

#############################
# Parameters for SNP simulation
#############################
SNPs_NR = 20
SNPs_EAF = 0.25
SNPs_centering = F
SNPs_classes = c(rep("A",SNPs_NR/2),rep("B",SNPs_NR/2))
SNPs_Correct_ASs = T

#############################
# Parameters for exposure simulation
#############################
X_beta_mean = c(0.2,0.2)
X_beta_sd = c(0.05,0.05)
X_mean_random = c(0,1)
X_var_random = c(0.1,0.1)
X_covar_random = 0
X_error_mean = 0
X_error_sd = 0.1

X_func0 = function(t,g,u) g+0*t+u
X_func1 = function(t,g,u) g*sin(u+t)
X_func2 = function(t,g,u) sin(t*pi*g*u/180)
X_func3 = function(t,g,u) t*g*u/60                   

X_useAs1stFunction = "X_func0"
X_useAs2ndFunction = "X_func2" # possible functions X_func1, X_func2, X_func3

#############################
# Parameters for outcome simulation
#############################
Y_alpha = c(0.3,0.3)
Y_mean_random = 0
Y_var_random = 1
Y_AS_new = F # get new AS for variability based on GX association with SD, eigenfunction 2, interaction term or variance (depending on AssocModel)

#############################
# Parameters for SNP association 
#############################
AssocModel = "meanSD"      # possible models: meanSD, eigenfunc, linMixed, gamlss
linMixed_random = F        # only relevant in linMixed assoc model, should these time points be selected randomly or not
linMixed_NRtimepoints = 5  # only relevant in linMixed assoc model and random = T, how many time points should be used for lineare mixed regression
eigenfunc_cutoff = 0.01

#############################
# Parameters for MVMR
#############################
MR_filterBadSNPs = T
MR_filterBadSNPs_treshold = 1e-6
MR_doCorrection = T

