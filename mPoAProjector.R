# projector.R -- 20180522
# 
# Given a set of methylation beta values and directory containing the probe
# model weights, these functions compute the model's predictive score.
# 
# Usage:
#   source("mPoAProjector.R")
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

projector = function( betas, proportionOfProbesRequired=0.8, outputDirectory="./", modelEnvironmentLocation="" ) {
  # Remotely load model environment
#  load_url(modelEnvironmentLocation)
  
  # loop through models
  model_results <- lapply(mPOA_Models$model_names, function(model_name) {
    # make sure it has been converted to a matrix
    if( !is.numeric(as.matrix(betas)) ) { stop("betas matrix/data.frame is not numeric!") }
    probeOverlap <- length(which(rownames(betas) %in% mPOA_Models$model_probes[[model_name]])) / length(mPOA_Models$model_probes[[model_name]])
    # make sure enough of the probes are present in the data file
    if( probeOverlap < proportionOfProbesRequired ) { 
      result <- rep(NA, ncol(betas))
      names(result) <- colnames(betas)
      result
    } else {
      # Work with a numeric matrix of betas
      betas.mat <- as.matrix(betas[which(rownames(betas) %in% mPOA_Models$model_probes[[model_name]]),])
      # If probes don't exist, we'll add them as rows of 'NA's
      probesNotInMatrix <- mPOA_Models$model_probes[[model_name]][which(mPOA_Models$model_probes[[model_name]] %in% rownames(betas.mat) == F)]
      if( length(probesNotInMatrix) > 0 ) {
        for( probe in probesNotInMatrix ) {
          tmp.mat <- matrix(NA, nrow=1, ncol=ncol(betas.mat))
          rownames(tmp.mat) <- probe
          colnames(tmp.mat) <- colnames(betas.mat)
          betas.mat <- rbind(betas.mat, tmp.mat)
        }
      }
      
      # Identify samples with too many missing probes and remove them from the matrix
      samplesToRemove <- colnames(betas.mat)[which(apply(betas.mat, 2, function(x) { length(which(is.na(x))) / length(x) > proportionOfProbesRequired}))]
      if( length(samplesToRemove) > 0 ) {
        betas.mat <- betas.mat[,-which(colnames(betas.mat) %in% samplesToRemove)]
      }
      if(ncol(betas.mat) > 0) { 
        # Identify missingness on a probe level
        pctValuesPresent <- apply( betas.mat, 1, function(x) { 1 - (length(which(is.na(x))) / length(x)) } )
        # If they're missing values, but less than the proportion required, we impute to the cohort mean
        if( length(which(pctValuesPresent < 1 & pctValuesPresent >= proportionOfProbesRequired)) > 0 ) {
          betas.mat[which(pctValuesPresent < 1 & pctValuesPresent >= proportionOfProbesRequired),] <- t(apply( betas.mat[which(pctValuesPresent < 1 & pctValuesPresent >= proportionOfProbesRequired),], 1 , function(x) { 
            x[is.na(x)] = mean( x, na.rm = TRUE )
            x
          }))
        }
        # If they're missing too many values, everyones value gets replaced with the mean from the Dunedin cohort
        if( length(which(pctValuesPresent < proportionOfProbesRequired)) > 0 ) {
          probesToReplaceWithMean <- rownames(betas.mat)[which(pctValuesPresent < proportionOfProbesRequired)]
          for( probe in probesToReplaceWithMean ) {
            betas.mat[probe,] <- rep(mPOA_Models$model_means[[model_name]][probe], ncol(betas.mat))
          }
        }
        # Calculate score:
        score = mPOA_Models$model_intercept[[model_name]] + rowSums(t(betas.mat[mPOA_Models$model_probes[[model_name]],]) %*% diag(mPOA_Models$model_weights[[model_name]]))
        names(score) <- colnames(betas.mat)
        if( length(samplesToRemove) > 0 ) {
          score.tmp <- rep(NA, length(samplesToRemove))
          names(score.tmp) <- samplesToRemove
          score <- c(score, score.tmp)
        }
        score <- score[colnames(betas)]        
        score
      } else {
        result <- rep(NA, ncol(betas.mat))
        names(result) <- colnames(betas.mat)
        result
      }
    }
  })
  names(model_results) <- mPOA_Models$model_names
  if( length(outputDirectory) > 0 & is.na(outputDirectory) == FALSE ) {
    dir.create(outputDirectory, recursive=TRUE, showWarnings=FALSE)
    sapply(mPOA_Models$model_names, function(model_name) {
      results.df <- data.frame(SampleID=names(model_results[[model_name]]), mPoA=model_results[[model_name]])
      write.table(results.df, file=paste0(outputDirectory, "/", model_name, "_results.csv"), sep=',', quote=FALSE, row.names=FALSE)
    })
  }
  model_results
}

### Pulled from https://stackoverflow.com/questions/24846120/importing-data-into-r-rdata-from-github
### Code to be able to load the environment file from github
load_url <- function (url, ..., sha1 = NULL) {
  # based very closely on code for devtools::source_url
  stopifnot(is.character(url), length(url) == 1)
  temp_file <- tempfile()
  on.exit(unlink(temp_file))
  request <- httr::GET(url)
  httr::stop_for_status(request)
  writeBin(httr::content(request, type = "raw"), temp_file)
  file_sha1 <- digest::digest(file = temp_file, algo = "sha1")
  if (is.null(sha1)) {
    message("SHA-1 hash of file is ", file_sha1)
  }
  else {
    if (nchar(sha1) < 6) {
      stop("Supplied SHA-1 hash is too short (must be at least 6 characters)")
    }
    file_sha1 <- substr(file_sha1, 1, nchar(sha1))
    if (!identical(file_sha1, sha1)) {
      stop("SHA-1 hash of downloaded file (", file_sha1, 
           ")\n  does not match expected value (", sha1, 
           ")", call. = FALSE)
    }
  }
  load(temp_file, envir = .GlobalEnv)
}



