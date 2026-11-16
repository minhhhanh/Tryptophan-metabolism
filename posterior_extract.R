extract_longitudinal <- function(jm,
                                 metabolite_idx,
                                 data.long,
                                 time = c(0, 30),
                                 fixed_formula = ~ day * arm + AGE + SEX,
                                 random_formula = ~ day,
                                 quantity = c("prediction", "slope"),
                                 summary = c("median", "mean"),
                                 effect = c("conditional", "random"),
                                 merge_obs = TRUE,
                                 keep.iter = NULL) {
  
  quantity <- match.arg(quantity)
  effect <- match.arg(effect)
  summary <- match.arg(summary)
  
  summarize_draws <- function(draw_matrix, summary) {
    estimate <- switch(
      summary,
      median = apply(draw_matrix, 1, median),
      mean = apply(draw_matrix, 1, mean)
    )
    
    data.frame(
      estimate = estimate,
      lwr = apply(draw_matrix, 1, quantile, .025),
      upr = apply(draw_matrix, 1, quantile, .975)
    )
  }
  
  subjects <- unique(as.character(jm$model_data$idT))
  n_subjects <- length(subjects)
  n_chains <- length(jm$mcmc$b)
  
  if (is.null(keep.iter))
    keep.iter <- nrow(jm$mcmc$betas1[[1]])
  iter_use <- seq_len(keep.iter)
  
  subj_info <- data.long %>%
    group_by(Subject_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(Subject_id, AGE, SEX, arm);
  subj_info <- subj_info[match(subjects, subj_info$Subject_id), ]
  
  met_name_model <- names(jm$model_data$extra_parms)[metabolite_idx]
  met_name_data <- names(data.long)[make.names(names(data.long)) == met_name_model]
  
  if (length(met_name_data) != 1)
    stop("Could not uniquely map metabolite name.")
  
  if (quantity == "slope") {
    slope_draws <- simplify2array(lapply(iter_use, function(iter) {
      do.call(cbind, lapply(seq_len(n_chains), function(chain) {
        
        b_it <- jm$mcmc$b[[chain]][, slope_idx, iter] # n_subjects-length vector
        
        if (effect == "conditional") {
          beta <- jm$mcmc[[paste0("betas", metabolite_idx)]][[chain]][iter, ]
          beta_names <- names(beta)
          
          beta_day <- beta[grep("^day$", beta_names)]
          beta_day_arm <- beta[grep("^day:arm|^arm.*:day", beta_names)]
          
          if (length(beta_day) != 1)
            stop("Couldn't uniquely identify the day coefficient.")
          
          if (length(beta_day_arm) == 0) {
            fixed_part <- rep(beta_day, n_subjects)
          } else {
            dex <- subj_info$arm == "Dexamethasone"
            fixed_part <- beta_day + beta_day_arm * as.numeric(dex)
          } # n_subjects-length vector
          pred <- b_it + fixed_part
        } else {
          pred <- b_it
        }
        names(pred) <- subjects
        pred
        
      }))
    }))
    
    slope_mat <- matrix(aperm(slope_draws, c(1, 3, 2)), nrow = n_subjects);
    colnames(slope_mat) <- paste0("draw_", seq_len(ncol(slope_mat)))
    
    slope_est <- summarize_draws(slope_mat, summary)
    
    slope_df <- data.frame(
      Subject_id = subjects,
      slope = slope_est$estimate,
      lwr = slope_est$lwr,
      upr = slope_est$upr)
    
    return(list(slope = slope_df, slope_draws = slope_mat))
    
  } else {
    
    X_pred <- expand.grid(Subject_id = subjects, day = time, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE);
    X_pred <- X_pred %>% left_join(subj_info, by = "Subject_id");
    X_pred$arm <- factor(X_pred$arm, levels = levels(data.long$arm));
    X_pred$SEX <- factor(X_pred$SEX, levels = levels(data.long$SEX))
    
    Zmat <- model.matrix(random_formula, X_pred)
    Xmat <- model.matrix(fixed_formula, X_pred)
    
    K <- ncol(Zmat) # the numbers of random effects
    RE_idx <- ((metabolite_idx - 1) * K + 1):(metabolite_idx * K)
    subject_index <- match(X_pred$Subject_id, subjects)
    
    pred_draws <- simplify2array(lapply(iter_use, function(iter) {
      do.call(cbind, lapply(seq_len(n_chains), function(chain) {
        b_it <- jm$mcmc$b[[chain]][, RE_idx, iter, drop = FALSE] # drop = FALSE preserve the third dimension
        b_it <- matrix(b_it, ncol = K)
        rand_part <- rowSums(Zmat * b_it[subject_index, , drop = FALSE])
        
        beta_it <- jm$mcmc[[paste0("betas", metabolite_idx)]][[chain]][iter, ]
        fix_part <- as.vector(Xmat %*% beta_it)
        fix_part + rand_part
      }))
    }))
    
    # Convert: prediction × chain × iteration -> prediction × posterior_draw
    pred_draws <- matrix(aperm(pred_draws, c(1, 3, 2)), nrow = nrow(X_pred))
    value_est <- summarize_draws(pred_draws, summary)
    
    traj_df <- data.frame(
      Subject_id = X_pred$Subject_id,
      day = X_pred$day,
      arm = X_pred$arm,
      value = value_est$estimate,
      lwr = value_est$lwr,
      upr = value_est$upr
    )
    if (merge_obs) {
      obs_df <- data.long %>% dplyr::select(Subject_id, day, arm, all_of(met_name_data))
      final_df <- traj_df %>% left_join(obs_df, by = c("Subject_id", "day", "arm")) %>%
        dplyr::rename(pred = value, 
                      obs = all_of(met_name_data)) %>%
        dplyr::mutate(
          value = if_else(!is.na(obs), obs, pred),
          source = if_else(!is.na(obs), "observed", "JM_imputed"),
          metabolite = met_name_data
        )
    } else {
      cat(
        "\n==================",
        "\nNot merged with observational data.",
        "\nReturning JM-inferred data only.",
        "\n==================\n"
      )
      final_df <- NULL
    }
    
    pop_grid <- expand.grid(day = sort(unique(time)), arm = levels(X_pred$arm));
    pop_grid$AGE <- mean(subj_info$AGE, na.rm = TRUE);
    pop_grid$SEX <- factor(names(which.max(table(subj_info$SEX))), levels = levels(subj_info$SEX))
    
    X_pop <- model.matrix(fixed_formula, pop_grid)
    
    pop_draws <- simplify2array(lapply(iter_use, function(iter) {
      do.call(cbind, lapply(seq_len(n_chains), function(chain) {
        beta_it <- jm$mcmc[[paste0("betas", metabolite_idx)]][[chain]][iter, ]
        as.vector(X_pop %*% beta_it)
      }))
    }))
    
    pop_draws <- matrix(aperm(pop_draws, c(1, 3, 2)), nrow = nrow(X_pop))
    pop_est <- summarize_draws(pop_draws, summary)
    
    pop_df <- data.frame(
      day = pop_grid$day,
      arm = pop_grid$arm,
      value = pop_est$estimate,
      lwr = pop_est$lwr,
      upr = pop_est$upr,
      metabolite = met_name_data
    )
    return(
      list(
        prediction_grid = X_pred,
        traj_df = traj_df,
        final_df = final_df,
        pop_df = pop_df,
        pred_draws = pred_draws,
        pop_draws = pop_draws
      )
    )
  }
}

