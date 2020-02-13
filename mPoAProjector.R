# projector.R -- 20180522
# 
# Given a set of methylation beta values and directory containing the probe
# model weights, these functions compute the model's predictive score.
# 
# Usage:
#   source("projector.R")
#   load("betas")
#   project(betas)
#
# Input:
#   betas:
#     Matrix or data.frame of beta values where rownames are probe ids and
#     column names should correspond to sample names.
#     N.B. ensure beta values are numeric
#   
#   weights:
#     Two column tab separated file, containing probes and weights of the model.
#     a value for "Intercept" is expected.
#     In this instance, this is the provided file "poa-all.tsv"
#  
#   modeldir:
#     Directory containing weights for all models.
#     In this instance, this is the directory "mpoa_models", containing one file,
#     and should be called explicitly under projector function.
#  
#   outputdir: 
#     Directory to save output, created if it does not exist.
#  
# Output:
#   For each model a file of the following form will be generated:
#   <outputdir>/<model>.csv
#
####################################################################################################

# First run the following functions:


projector = function( betas, modelEnvironmentLocation="", pctProbesRequired=0.9 ) {
  # Needed to download model environment
  require(RCurl)
  
  # Download environment
  
  
  
  
}

projector = function(betas, modeldir="models", outdir="results") {
  models = read_models(modeldir)
  betas = t(check_betas(betas, models))

  # Iterate through models and write scores
  for (modelname in names(models)) {
    model = models[[modelname]]
    weights = model$weights
    found = intersect(weights$probe, colnames(betas))
    score = model$intercept + (betas[, found] %*% weights[found, "weight"])
    write_score(score, found, model, modelname, outdir)
  }
}

read_models = function(path) {
  # Reads model weights
  files = list.files(path, pattern=".tsv", full.names=TRUE)
  names(files) = sub(".tsv", "", basename(files))
  ret = lapply(files, function(fn) {
    dat = read.csv(fn, sep="\t", header=FALSE, col.names=c("probe", "weight"), as.is=TRUE)
    rownames(dat) = dat$probe
    list(
      "intercept"=dat["Intercept", "weight"],
      "weights"=dat[setdiff(dat$probe, "Intercept"),])
  })
  return(ret)
}

check_betas = function(betas, models) {
  # Remove probes not in any model and deals with missing values
  probes = unique(do.call(c, lapply(models, function(x) x$weights$probe)))
  ret = betas[which(rownames(betas) %in% probes),]
  ret = check_na(ret)
  return(ret)
}

check_na = function(betas) {
  # Removes probes with missing values in greater than 5% of samples 
  # Imputes missing values for probes with at least 95% of samples 
  maxmis = ncol(betas) * 0.05
  hasna = apply(betas, 1, function(x) length(which(is.na(x))))

  remove = names(which(hasna > maxmis))
  if (length(remove) > 0) {
    betas = betas[setdiff(rownames(betas), remove), ]
  }

  impute = names(which(hasna < maxmis & hasna > 0))
  if (length(impute) > 0) {
    betas[impute, ] = apply(betas[impute, ], 1, function(x) {
      x[is.na(x)] = mean(x, na.rm=TRUE)
    })
  }
  return(betas)
}

write_score = function(score, found, model, modelname, outdir) {
  res = data.frame("id"=rownames(score), "score"=score[,1])
  header = paste0(
    "# Using ", length(found), " of ", nrow(model$weights), 
    " selected probes\n# Intercept: ", model$intercept, "\n")

  filename = file.path(outdir, paste0(modelname, ".csv"))
  dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
  cat(header, file=filename, append=FALSE)
  write.table(
    res, file=filename, sep=',', append=TRUE, quote=FALSE, 
    row.names=FALSE, col.names=FALSE)
}

####################################################################################################

# To generate the mPoA values, call the function "projector" with appropriate arguments.
# For example, if DNA methylation beta values object is named "betas", 
# call projector(betas, "mpoa_models")

###################################################################################################
