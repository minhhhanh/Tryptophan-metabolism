# The random slope is the subject-specific deviation from the population slope,
# not the subject-specific slope itself. Therefore, we reconstruct each
# subject's posterior slope as:
#   population slope + treatment-specific slope modification + random slope.
# Because the population slope differs between treatment arms, the
# treatment-by-time interaction must also be included.

get_subject_slopes <- function(jm, data_long = dat_long, dat_surv = dat_surv, metabolite_idx, 
                               time_var = "day", trt_var = "armDexamethasone", trt_vector = NULL) {
  slope_idx_met <- 2 * metabolite_idx
  
  met_name_model <- names(jm$model_data$extra_parms)[metabolite_idx]
  met_name_data <- names(data_long)[make.names(names(data_long)) == met_name_model]
  
  if (length(met_name_data) != 1) {
    stop("Could not uniquely map metabolite name.")
  } else {
    cat(
      "\n==================",
      "\nComputing:", met_name_data,
      "\n==================\n"
    )
  }
  
  b_slopes <- do.call(rbind, lapply(jm$mcmc$b, function(b_chain) {
    t(b_chain[, slope_idx_met, ])
  }
  ))
  
  n_iter <- nrow(b_slopes);
  n_subj <- ncol(b_slopes)
  
  beta_mat <- do.call(rbind, lapply(jm$mcmc[[paste0("betas", metabolite_idx)]], as.matrix)) # 20k rows, 4 columns
  cn <- colnames(beta_mat)
  beta_day <- beta_mat[, cn == time_var] # a 20k long vector
  
  # treatment modified slope
  interaction_term <- paste0(time_var, ":", trt_var)
  idx <- which(cn == interaction_term) # index of interaction_term in the beta_mat
  
  if (length(idx) == 0) {
    beta_int <- rep(0, nrow(beta_mat))
  } else if (length(idx) == 1) {
    beta_int <- beta_mat[, idx] # 20k long vector
  } else {
    stop("Multiple interaction terms found")
  }
  
  ids <- levels(jm$model_data$idT)
  stopifnot(length(ids) == ncol(b_slopes))
  
  if (is.null(trt_vector)) {
    trt <- ifelse(dat_surv$arm == "Dexamethasone", 1, 0)
    trt_vector <- trt[match(ids, dat_surv$Subject_id)]
  }
  if (anyNA(trt_vector))
    stop("Some subjects could not be matched.")
  
  beta_day_mat <- matrix(beta_day,
                         nrow = n_iter,
                         ncol = n_subj)
  beta_int_mat <- beta_int %o% trt_vector
  # beta_int_mat <- beta_int %*% t(trt_vector)
  subject_slopes <- b_slopes + beta_day_mat + beta_int_mat
  colnames(subject_slopes) <- ids
  
  return(subject_slopes)
}