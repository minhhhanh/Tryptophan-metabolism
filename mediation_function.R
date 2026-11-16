prepare_mediator <- function(x, clip_quantile = NULL, scale = FALSE) {
  
  if (!is.null(clip_quantile)) {
    q <- quantile(x, c(clip_quantile, 1 - clip_quantile), na.rm = TRUE)
    x <- pmin(pmax(x, q[1]), q[2])
  }
  
  if (scale)
    x <- as.numeric(scale(x))
  
  x
}
default_cmest <- list(
  exposure = "arm",
  mediator = "mediator",
  outcome = "tt_death_12m",
  event = "ev_death_12m",
  EMint = TRUE,
  estimation = "imputation",
  inference = "bootstrap"
)

fit_mediation <- function(dat, cmest_args) {
  
  old_opt <- options(warn = 2)
  on.exit(options(old_opt), add = TRUE)
  do.call(
    CMAverse::cmest,
    c(list(data = dat),
      cmest_args
    )
  )
}
safe_extract <- function(res, row, col) {
  if (row %in% rownames(res) && col %in% colnames(res)) {
    return(res[row, col])
  } else {
    return(NA)
  }
}
empty_result <- function(draw, warning = NA) {
  data.frame(
    draw = draw,
    TE_est = NA_real_,
    TNDE_est = NA_real_,
    PNDE_est = NA_real_,
    TNIE_est = NA_real_,
    PNIE_est = NA_real_,
    PM_est = NA_real_,
    warning = warning
  )
}
extract_estimates <- function(res) {
  effects <- c(
    TE = "Rte",
    TNDE = "Rtnde",
    PNDE = "Rpnde",
    TNIE = "Rtnie",
    PNIE = "Rpnie",
    PM = "pm"
  )
  out <- lapply(effects, \(x)
                safe_extract(res, x, "Estimate"))
  names(out) <- paste0(names(effects), "_est")
  as.data.frame(out)
}
summarize_effect <- function(x) {
  tibble::tibble(
    est = median(x, na.rm = T),
    lwr = quantile(x, 0.025, na.rm = T),
    upw = quantile(x, 0.975, na.rm = T)
  )
}
run_mediation_simplify <- function(slope_mat,
                          analysis_df,
                          cmest_args,
                          n_draws = 2000,
                          nboot = 500,
                          basec = c("AGE", "SEX", "LTA4H", "TBMGRADE", "CSFWBC", "CSFPROTEIN"),
                          scale_mediator = FALSE,
                          clip_quantile = NULL,
                          seed = 1234) {
  if (!requireNamespace("CMAverse", quietly = TRUE))
    stop("Package 'CMAverse' is required.")
  set.seed(seed)
  # Prepare data
  subjects <- colnames(slope_mat)
  analysis_df <- analysis_df |>
    dplyr::distinct(Subject_id, .keep_all = TRUE) |>
    dplyr::mutate(arm = factor(arm, levels = c("Placebo", "Dexamethasone")))
  idx <- match(subjects, analysis_df$Subject_id)
  
  if (anyNA(idx)) {
    stop("Missing Subject_id(s): ", paste(subjects[is.na(idx)], collapse = ", "))
  }
  
  analysis_df <- analysis_df[idx, ]
  
  stopifnot(identical(analysis_df$Subject_id, subjects))
  
  n_draws <- min(n_draws, nrow(slope_mat))
  draw_idx <- sample.int(nrow(slope_mat), size = n_draws, replace = FALSE)
  
  res_list <- vector("list", length(draw_idx))
  pb <- txtProgressBar(min = 0, max = n_draws, style = 3)
  on.exit(close(pb), add = TRUE)
  
  cmest_args <- modifyList(
    default_cmest,
    list(
      model = "gformula",
      basec = basec,
      mreg = list("linear"),
      yreg = "coxph",
      nboot = nboot,
      a = "Dexamethasone",
      astar = "Placebo",
      mval = list(median(dat_tmp$mediator, na.rm = TRUE))
    )
  )
  for (k in seq_along(draw_idx)) {
    setTxtProgressBar(pb, k)
    mediator <- prepare_mediator(slope_mat[draw_idx[k], ], clip_quantile, scale_mediator)
    
    dat_tmp <- analysis_df
    dat_tmp$mediator <- mediator
    fit <- tryCatch(
      fit_mediation(dat_tmp, cmest_args),
      error = function(e) list(error_msg = conditionMessage(e))
    )
    
    if (is.list(fit) && "error_msg" %in% names(fit)) {
      res_list[[k]] <- empty_result(draw_idx[k], fit$error_msg)
    } else {
      res_list[[k]] <- cbind(draw = draw_idx[k],
                             extract_estimates(summary(fit)$summarydf),
                             warning = NA)
    }
  }
  
  results <- do.call(rbind, res_list)
  effects <- c("TE", "TNDE", "TNIE", "PNDE", "PNIE", "PM")
  
  summary_df <- purrr::map_dfc(effects, function(e) {
    out <- summarize_effect(results[[paste0(e, "_est")]])
    names(out) <- c(e, paste0(e, "_lwr"), paste0(e, "_upr"))
    out
  })
  summary_df$n_success <- sum(!is.na(results$TE_est))
  summary_df$n_fail    <- sum(is.na(results$TE_est))
  summary_df$fail_rate <- mean(is.na(results$TE_est))
  
  return(list(summary = summary_df, draws = results))
}

#> Main function
run_mediation <- function(slope_mat, analysis_df, n_draws = 2000, nboot = 500, 
                          basec = c("AGE", "SEX", "LTA4H", "TBMGRADE", "CSFWBC", "CSFPROTEIN"),
                          scale_mediator = FALSE, clip_quantile = NULL, seed = 1234) {
  library(CMAverse)
  set.seed(seed)
  # Prepare data
  subjects <- colnames(slope_mat)
  analysis_df <- analysis_df |>
    dplyr::distinct(Subject_id, .keep_all = TRUE) |> 
    dplyr::mutate(arm = factor(analysis_df$arm, levels = c("Placebo", "Dexamethasone")))
  analysis_df <- analysis_df[match(subjects, analysis_df$Subject_id), ]
  stopifnot(identical(analysis_df$Subject_id, subjects))
  
  n_draws <- min(n_draws, nrow(slope_mat))
  draw_idx <- sample.int(nrow(slope_mat), size = n_draws, replace = FALSE)
  
  res_list <- vector("list", length(draw_idx))
  pb <- txtProgressBar(min = 0, max = n_draws, style = 3)
  
  for (k in seq_along(draw_idx)) {
    setTxtProgressBar(pb, k)
    mediator <- slope_mat[draw_idx[k], ]
    if (!is.null(clip_quantile)) {
      q <- quantile(mediator, c(clip_quantile, 1 - clip_quantile), na.rm = T)
      mediator <- pmin(pmax(mediator, q[1]), q[2])
    }
    
    if (scale_mediator) {
      mediator <- as.numeric(scale(mediator))
    }
    dat_tmp <- analysis_df
    dat_tmp$mediator <- mediator
    mval <- list(mean(dat_tmp$mediator, na.rm = T))
    fit <- tryCatch({
      old_opt <- options(warn = 2)
      on.exit(options(old_opt), add = T)
      CMAverse::cmest(
        data = dat_tmp,
        model = "gformula",
        exposure = "arm",
        mediator = "mediator",
        EMint = TRUE,
        basec = basec,
        mreg = list("linear"),
        yreg = "coxph",
        outcome = "tt_death_12m",
        event = "ev_death_12m",
        estimation = "imputation",
        inference = "bootstrap",
        nboot = nboot,
        a = "Dexamethasone",
        astar = "Placebo",
        mval = mval
      )
    }, error = function(e) {
      return(list(error_msg = conditionMessage(e)))
    })
    if (is.list(fit) && "error_msg" %in% names(fit)) {
      res_list[[k]] <- data.frame(
        draw = draw_idx[k],
        TE_est = NA,
        TNDE_est = NA,
        PNDE_est = NA,
        TNIE = NA,
        PNIE = NA,
        PM_est = NA,
        warning = fit$error_msg
      )
    } else {
      res <- summary(fit)$summarydf
      res_list[[k]] <- data.frame(
        draw = draw_idx[k],
        TE_est  = safe_extract(res, "Rte", "Estimate"),
        TNDE_est = safe_extract(res, "Rtnde", "Estimate"),
        PNDE_est = safe_extract(res, "Rpnde", "Estimate"),
        TNIE_est = safe_extract(res, "Rtnie", "Estimate"),
        PNIE_est = safe_extract(res, "Rpnie", "Estimate"),
        PM_est  = safe_extract(res, "pm", "Estimate"),
        warning = NA
      )
    }
  }
  close(pb)
  
  results <- do.call(rbind, res_list)
  summary_df <- data.frame(
    TE  = median(results$TE_est, na.rm = TRUE),
    TE_lwr  = quantile(results$TE_est, 0.025, na.rm = TRUE),
    TE_upr  = quantile(results$TE_est, 0.975, na.rm = TRUE),
    
    TNDE = median(results$TNDE_est, na.rm = TRUE),
    TNDE_lwr = quantile(results$TNDE_est, 0.025, na.rm = TRUE),
    TNDE_upr = quantile(results$TNDE_est, 0.975, na.rm = TRUE),
    
    TNIE = median(results$TNIE_est, na.rm = TRUE),
    TNIE_lwr = quantile(results$TNIE_est, 0.025, na.rm = TRUE),
    TNIE_upr = quantile(results$TNIE_est, 0.975, na.rm = TRUE),
    
    PNDE = median(results$PNDE_est, na.rm = TRUE),
    PNDE_lwr = quantile(results$PNDE_est, 0.025, na.rm = TRUE),
    PNDE_upr = quantile(results$PNDE_est, 0.975, na.rm = TRUE),
    
    PNIE = median(results$PNIE_est, na.rm = TRUE),
    PNIE_lwr = quantile(results$PNIE_est, 0.025, na.rm = TRUE),
    PNIE_upr = quantile(results$PNIE_est, 0.975, na.rm = TRUE),
    
    PM  = median(results$PM_est, na.rm = TRUE),
    fail_rate = mean(is.na(results$TE_est))
    
  )
  return(list(summary = summary_df, draws = results))
}
