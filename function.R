#> Function 1: CV QC
poolCV <- function(df, cv_cutoff = 0.3) {
  if (!"Method" %in% colnames(df))
    stop("`Method` column not found in input dataframe")
  method <- unique(df$Method)
  if (length(method) != 1) stop("Multiple methods detected — please subset by method before using this function.")
  if (method %in% c("HILIC-pos", "C8-pos")) {
    pooled_sample <- "PREFB"
    drift_sample  <- "PREFA"
  } else if (method %in% c("C18-neg", "HILIC-neg")) {
    pooled_sample <- "PREFA"
    drift_sample  <- "PREFB"
  } else {
    stop("Check `Method` column — unrecognized method: ", method)
  }
  PREP_cols <- grep(pooled_sample, colnames(df))
  PREP <- df[, PREP_cols, drop = FALSE]
  PREP <- as.data.frame(lapply(PREP, function(col) as.numeric(as.character(col))))
  cv_bio <- apply(PREP, 1, function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  })
  n_mets <- names(which(cv_bio >= cv_cutoff))
  return(n_mets)
}
#> Function 2: Half-min imputation
imputeHalfMin <- function(x){
  lod <- min(x, na.rm=T)
  x[is.na(x)] <- 0.5 * lod
  return(x)
}
# Function 3: Harmonizing levels in factored data
harmonize <- function(v1, v2) {
  all_levels <- union(levels(v1), levels(v2))
  v1 <- factor(v1, levels = all_levels)
  v2 <- factor(v2, levels = all_levels)
  list(v1 = v1, v2 = v2)
}
harmonize_dfs <- function(df_long, df_surv, vars) {
  for (v in vars) {
    if (!v %in% names(df_long)) {
      stop(sprintf("Variable '%s' not found in df_long", v))
    }
    if (!v %in% names(df_surv)) {
      stop(sprintf("Variable '%s' not found in df_surv", v))
    }
    
    tmp <- harmonize(df_long[[v]], df_surv[[v]])
    
    if (!is.list(tmp) || !all(c("v1", "v2") %in% names(tmp))) {
      stop(sprintf("harmonize() did not return list(v1, v2) for '%s'", v))
    }
    
    df_long[[v]] <- tmp$v1
    df_surv[[v]] <- tmp$v2
  }
  list(long = df_long, surv = df_surv)
}
#> Function 4: Standardizing the name in clinnical data
clean_names <- function(x) {
  x %>%
    str_to_upper() %>%
    str_replace_all("PROTEIN", "PRO") %>%
    str_replace_all("NEUTLE", "NEU") %>%
    str_replace_all("LYMLE", "LYM") %>%
    str_replace_all("GLUC", "GLU")
}
#> Function 5: Extracting slopes posteriors from a univariate JM object
extract_slope_posteriors <- function(jm_obj) {
  
  # stack chains × iterations into one posterior matrix
  beta_draws <- do.call(
    rbind,
    lapply(jm_obj$mcmc$betas1, as.matrix)
  )
  
  stopifnot(all(c("day", "day:armDexamethasone") %in% colnames(beta_draws)))
  
  beta_placebo <- beta_draws[, "day"]
  beta_delta   <- beta_draws[, "day:armDexamethasone"]
  beta_dex     <- beta_placebo + beta_delta
  
  data.frame(
    beta_placebo = beta_placebo,
    beta_dex     = beta_dex,
    delta_beta   = beta_delta
  )
}
#> Function 6: Summarizing beta into table
summarise_beta <- function(x) {
  c(
    median = median(x),
    lwr  = quantile(x, 0.025),
    upr  = quantile(x, 0.975),
    prob_pos = mean(x > 0)
  )
}
#> Function 7: Computing Wald hypothesis testing
Wald_P_cal <- function(jm_object,post_beta="day:armDexamethasone") {
  df <- summary(jm_object)$Outcome1
  ## calculate pvalue by wald test
  b=df[post_beta,"Mean"]
  sd=df[post_beta,"StDev"]
  z=b/sd
  options(scipen = 0, digits = 10)
  2 * (1 - pnorm(abs(z)))
}
#> Function 8: Extracting slopes posteriors from a multi-variate JM object
extract_slope_posteriors2 <- function(jm_obj) {
  
  # stack chains × iterations
  beta_draws <- do.call(
    rbind,
    lapply(jm_obj$mcmc$betas1, as.matrix)
  )
  
  stopifnot(all(c(
    "day",
    "day:armDexamethasone",
    "day:hiv_sttpositive"
  ) %in% colnames(beta_draws)))
  
  beta_day   <- beta_draws[, "day"]
  beta_arm   <- beta_draws[, "day:armDexamethasone"]
  beta_hiv   <- beta_draws[, "day:hiv_sttpositive"]
  
  data.frame(
    slope_dex_hivneg  = beta_day,
    slope_placebo_hivneg = beta_day + beta_arm,
    slope_dex_hivpos  = beta_day + beta_hiv,
    slope_placebo_hivpos = beta_day + beta_arm + beta_hiv,
    delta_arm = beta_arm,
    delta_hiv = beta_hiv
  )
}
#> Function 7: Computing Wald hypothesis testing for multivariate JM
Wald_P_cal2 <- function(jm_object, post_beta) {
  
  df <- summary(jm_object)$Outcome1
  
  b  <- df[post_beta, "Mean"]
  sd <- df[post_beta, "StDev"]
  
  z  <- b / sd
  
  2 * (1 - pnorm(abs(z)))
}
#> Function 8: Getting all covariance matrix from a MJM
get_D <- function(D_row, i, j) {
  if (i >= j) {
    return(D_row[paste0("D[", i, ", ", j, "]")])
  } else {
    return(D_row[paste0("D[", j, ", ", i, "]")])
  }
}
#> Function 9: covariance to correlation at time t
get_cov_ij_t <- function(D_row, i, j, t) {
  i0 <- intercept_idx[i]
  i1 <- slope_idx[i]
  j0 <- intercept_idx[j]
  j1 <- slope_idx[j]
  
  get_D(D_row, i0, j0) +
    t * get_D(D_row, i1, j0) +
    t * get_D(D_row, i0, j1) +
    t^2 * get_D(D_row, i1, j1)
}
get_cor_matrix_t <- function(D_row, t, met_names) {
  q <- length(met_names)
  cov_mat <- matrix(0, q, q)
  
  for (i in 1:q) {
    for (j in 1:q) {
      cov_mat[i, j] <- get_cov_ij_t(D_row, i, j, t)
    }
  }
  cor_mat <- cov2cor(cov_mat)
  rownames(cor_mat) <- met_names
  colnames(cor_mat) <- met_names
  return(cor_mat)
}
#> Function 10: Getting slope by group
get_slope_by_group <- function(b_iter, beta_day, beta_int, group_vec) {
  library(coda)
  q <- ncol(b_iter) / 2
  slope_idx <- seq(2, 2*q, by = 2)
  
  slopes_random <- b_iter[, slope_idx, drop = FALSE]
  
  idx_ctrl <- which(group_vec == "Placebo")
  idx_trt  <- which(group_vec == "Dexamethasone")
  
  slopes_ctrl <- sweep(slopes_random[idx_ctrl, , drop = FALSE], 2, beta_day, "+")
  slopes_trt  <- sweep(slopes_random[idx_trt,  , drop = FALSE], 2, beta_day + beta_int, "+")
  
  list(
    control = slopes_ctrl,
    treated = slopes_trt
  )
}
#> Function 11: Calculating significance level of correlations 
get_pmat <- function(mat) {
  q <- ncol(mat)
  pmat <- matrix(NA, q, q)
  
  for(i in 1:q) {
    for(j in 1:q) {
      pmat[i,j] <- cor.test(mat[,i], mat[,j])$p.value
    }
  }
  
  rownames(pmat) <- colnames(mat)
  colnames(pmat) <- colnames(mat)
  
  pmat
}
#> Function 12: Calculating slope-slope correlations
slope_cor <- function(b_array, beta_list, group_vec, met_names) {
  
  n_iter <- dim(b_array)[1]
  q <- length(met_names)
  
  cor_ctrl <- vector("list", n_iter)
  cor_trt  <- vector("list", n_iter)
  
  p_ctrl <- vector("list", n_iter)
  p_trt  <- vector("list", n_iter)
  
  for (m in 1:n_iter) {
    
    b_iter <- b_array[m, , , drop = FALSE]
    b_iter <- b_iter[1,,]
    
    slopes <- get_slope_by_group(
      b_iter,
      beta_list$day[m, ],
      beta_list$interaction[m, ],
      group_vec
    )
    cor_ctrl[[m]] <- cor(
      slopes$control,
      use = "pairwise.complete.obs"
    )
    
    cor_trt[[m]] <- cor(
      slopes$treated,
      use = "pairwise.complete.obs"
    )
    
    p_ctrl[[m]] <- get_pmat(slopes$control)
    
    p_trt[[m]] <- get_pmat(slopes$treated)
    
    rownames(cor_ctrl[[m]]) <- met_names
    colnames(cor_ctrl[[m]]) <- met_names
    
    rownames(cor_trt[[m]]) <- met_names
    colnames(cor_trt[[m]]) <- met_names
    
    rownames(p_ctrl[[m]]) <- met_names
    colnames(p_ctrl[[m]]) <- met_names
    
    rownames(p_trt[[m]]) <- met_names
    colnames(p_trt[[m]]) <- met_names
  }
  
  return(list(
    control = list(
      cor = cor_ctrl,
      p = p_ctrl
    ),
    treated = list(
      cor = cor_trt,
      p = p_trt
    )
  ))
}
#> Function 14: Calculating value for a metabolite at time t
get_Y_t <- function(b_mat, t, intercept_idx, slope_idx) {
  intercepts <- b_mat[, intercept_idx, drop = FALSE]
  slopes <- b_mat[, slope_idx, drop = FALSE]
  intercepts + t*slopes
}


