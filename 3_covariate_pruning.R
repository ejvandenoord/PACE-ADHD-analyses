#Covariate pruning
#Code written by Mandy Meijer, Karolina Aberg and Edwin van den Oord, 2020. Center for Biomarker Research and Precision Medicine (BPM), Virginia Commonwealth University
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#This script needs to be iterated
#You will decide which additional covariates to use and how many principal components  for the final MWAS model
#First, we will run the MWAS without any covariates, after which we will include more and more covariates
#Based on the correlation between ( EvdO corrected TYPO)  covariates with your outcome, the optimal amount of covariates will be selected for all three phenotypes. Eventually PCs with remaining variance will be added.

dirproject = "/nfs/home2/mmeijer/ramwas/adhd" #provide here your project directory
collect_files   = c("PC_vs_covs_pvalue.txt","PC_plot_covariates_removed.pdf", "UsedSettings.txt", "screeplot_techPC_values.pdf", "QQ_plot.pdf", "DegreesOfFreedom.txt") # these are the files we want to collect now

# read covariates. 
covariates <- read.csv(paste0(dirproject,"/covariates/covariates_ramwas_PACE_ADHD.csv") ) 

modeloutcome <- "ADHD_z-score" #or Inatt_z-score or Hyp_z-score

#This needs to be set for every iteration -- until we start the RAMWAS itself
techPCs <- 2 #put here the amount of technical PCs you want to include

#Here we will select the covariates we will put in the model. First run the model with no covariates (NULL) and run the script to get a feeling of what covariates could be important for your cohort. 
#Below we will provide some examples on which covariates could be important, but decide yourself based on the PC_vs_covs_pvalue.txt file (see below)
modelcovariates <- NULL 
#modelcovariates <- c("Age", "Sex", "Smoking", "Ethnicity",  
#                     "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")

#modelcovariates = c(modelcovariates, "Array", "xSlide")

#modelcovariates = c(modelcovariates,"mMed", "uMed",paste0("techPC",1:techPCs))


analysis_label = "GSMS_no_covariates_ADHD_bulk" #name for the  CURRENT analysis. The new subdirectory in your project directory will have this name. Please make sure that it includes your cohortname

modelPCs <- 0 #insert here the number of PCs you decided to use (start with 0, increase in new iteration only after you decided on the covariates)


### end settings

### start RAMWAS
# In the subdirectory paste0(param$dirproject,analysis_label,"\\"), the following files can be found:
# File with your used settings; QQ-plot; p-values for correlation between covariates and PCs; model PC plot; Manhattan plot; file with degrees of freedom

library(zip) 
library(ramwas)
setwd(dirproject) 
source("0_functions_PACE_ADHD.R")


#Run the model and inspect PC_vs_covs_pvalues.txt, to include significant extra covariates that correlate with your model PCs. Do not include them all at once, but one at a time -> check again, until you found an optimal amount of covariates
#Note here that the p-value threshold doesn't have to be 0.05. Use common sense, and determine how the model improves. 
#When you found the right amount of covariates, look at the scree plot of the model PCs. Are there any PCs that still explain variation in the data? Also add these to the model
#We will use the same covariates for all phenotypes (ADHD_z-score; Inatt_z-score; Hyp_z-score)


#MWAS with covariates
#This needs to be iterated (see description above)

param = ramwasParameters(
  dircoveragenorm = 'rw',
  covariates = covariates,
  modelcovariates = modelcovariates,
  modeloutcome = modeloutcome,
  modelPCs = modelPCs)


message("Scree plots for the PCA on technical probes and methylation data are in: ",paste0(dirproject,"/",analysis_label,"/"))


ramwas4PCA(param) #The covariates in the model will be regressed out before every PCA, so this line of code needs to be run every time you change something to the model
ramwas5MWAS(param)
qqPlotFast(getMWAS(param)$`p-value`)
title(paste("QQ_plot_", modeloutcome, "_", modelcovariates, "_", modelPCs))
#If you see in the QQ-plot that the lambda is not close to 1, go back to the model to see if there are not other covariates or model PCs that need to be included 

# Collecting the files


zip_file = paste0(param$dirproject,analysis_label,"/analysis_label","_results.zip") #This is the name of the file that you need to send to us - we only want to receive the folder of the model for the three phenotypes that you chose to be the best one
collect_results(param,analysis_label,collect_files) # This will collect the results and archive them


