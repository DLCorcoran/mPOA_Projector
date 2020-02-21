# mPOA_Projector
Pace of Age calculator for Illumina methyl-array data

mPoAProjector.R -- 20200220
 
Given a set of methylation beta values and directory containing the probe
model weights, these functions compute the model's predictive score.
 
#### Usage:
  source("mPoAProjector.R")  
  load("betas")  
  project(betas) -> mPoA_results.list  
  
## Input:
####  betas:
    Matrix or data.frame of beta values where rownames are probe ids and  
    column names should correspond to sample names.  
    N.B. ensure beta values are numeric  
    Missing values should be coded as 'NA'  

####  outputDirectory:  
    This should be a directory name to save the output (a .csv file per model analyzed).  
    By default, it is set to create a 'results' subdirectory in the current folder.  
    Changing this parameter to be blank "", or NA, will result in no files being generated  
  
####  proportionOfProbesRequired:  
    This is the proportion of probes to have a non-missing value for both the sample to have a  
    mPoA calculated, as well as to determine if we can impute the mean from the current cohort  
    By default, this is set to 0.8  
  
####  modelEnvironmentLocation:   
     This is the remote location of the environment containing the data for the models.  Do not adjust!  
  
## Output:  
  [1] A list containing the mPoAs for each model  
  [2] If a <outputDirectory> is specified, it will generate, for each model,   
      a file of the following form: <outputdir>/<model>_results.csv  

