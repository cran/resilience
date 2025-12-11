#' Prep-Resilience Analysis with Pre-Imputed Datasets
#'
#' This function takes a list of already imputed datasets and runs resilience analysis
#' on each, then pools the results using Rubin's rules. Users are responsible for
#' the imputation process, allowing maximum flexibility in imputation methods.
#'
#' @param data_list List of imputed datasets (data frames) or a `mids` object from `mice`
#' @param formula Formula for the prepost model (e.g., post_score ~ pre_score + age + bmi)
#' @param k Resilience parameter (default = 1.2)
#' @param nboot Number of bootstrap samples (default = 200)
#' @param pool_method Pooling method: "detailed" (Barnard & Rubin) or "simple" (default = "detailed")
#' @param coef_labels Optional named vector for renaming coefficients (e.g., c("age" = "Age (years)"))
#' @param verbose Whether to print progress messages (default = TRUE)
#' @param ... Additional arguments passed to prepost()
#'
#' @return A list containing:
#'   - pooled_results: The pooled estimates
#'   - results_table: Formatted results table
#'   - individual_results: List of results from each imputation
#'   - model_info: Information about the model specification
#'   - imputation_info: Summary of imputation characteristics
#'
#'@seealso 
#' \code{\link{prepost}} for complete-case analysis without missing data
#' \code{\link[mice]{mice}} for multiple imputation
#' 
#' @examples
#' \dontrun{
#' # Example 1: List of data frames from any imputation method
#' library(mice)
#' data <- data.frame(
#'   post_score = rnorm(100),
#'   pre_score = rnorm(100),
#'   age = rnorm(100, 50, 10),
#'   bmi = rnorm(100, 25, 5)
#' )
#' 
#' # Create some missing data
#' data$age[1:10] <- NA
#' data$bmi[5:15] <- NA
#' 
#' # User does their own imputation (could be from any package/method)
#' # Method 1: Using mice
#' imp_mice <- mice(data, m = 5, printFlag = FALSE)
#' imp_list <- list()
#' for (i in 1:5) {
#'   imp_list[[i]] <- complete(imp_mice, i)
#' }
#' 
#' # Method 2: Using Amelia
#' if (requireNamespace("Amelia", quietly = TRUE)) {
#'   imp_amelia <- Amelia::amelia(data, m = 5)
#'   imp_list <- imp_amelia$imputations
#' }
#' 
#' # Method 3: Using aregImpute from Hmisc
#' if (requireNamespace("Hmisc", quietly = TRUE)) {
#'   set.seed(123)
#'   imp_areg <- Hmisc::aregImpute(~ post_score + pre_score + age + bmi, data = data, n.impute = 5)
#'   imp_list <- list()
#'   for (i in 1:5) {
#'     imp_data <- data
#'     imp_data$age[is.na(imp_data$age)] <- imp_areg$imputed$age[, i]
#'     imp_data$bmi[is.na(imp_data$bmi)] <- imp_areg$imputed$bmi[, i]
#'     imp_list[[i]] <- imp_data
#'   }
#' }
#' 
#' # Run resilience analysis on the pre-imputed list
#' result <- prepost_resilience_mi_list(
#'   data_list = imp_list,
#'   formula = post_score ~ pre_score + age + bmi,
#'   k = 1.2,
#'   nboot = 200
#' )
#' 
#' # Example 2: Directly using a mice mids object
#' result2 <- prepost_resilience_mi_list(
#'   data_list = imp_mice,  # mids object
#'   formula = post_score ~ pre_score + age + bmi
#' )
#' }
#'
#' @export
prepost_resilience_mi_list <- function(data_list, 
                                       formula, 
                                       k = 1.0, 
                                       nboot = 500,
                                       pool_method = "detailed",
                                       coef_labels = NULL,
                                       verbose = TRUE,
                                       ...) {
  
  if (!requireNamespace("resilience", quietly = TRUE)) {
    stop("Package 'resilience' is required but not installed.")
  }
  
  library(resilience)
  
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object.")
  }
  
  if (!pool_method %in% c("detailed", "simple")) {
    stop("pool_method must be either 'detailed' or 'simple'")
  }
  
  # Check if data_list is a mids object from mice
  if (inherits(data_list, "mids")) {
    if (verbose) cat("Detected 'mids' object from mice package\n")
    m <- data_list$m
    
    if (verbose) {
      cat(sprintf("Converting %d imputed datasets from mice object\n", m))
    }
    
    # Convert mids object to list of data frames
    imp_list <- list()
    for (i in 1:m) {
      imp_list[[i]] <- mice::complete(data_list, i)
    }
    data_list <- imp_list
    
    mids_object <- data_list
  } else if (is.list(data_list)) {
    if (length(data_list) == 0) {
      stop("data_list is an empty list")
    }
    
    # Check if all elements are data frames
    are_data_frames <- sapply(data_list, is.data.frame)
    if (!all(are_data_frames)) {
      stop("All elements in data_list must be data frames")
    }
    
    mids_object <- NULL
  } else {
    stop("data_list must be either a list of data frames or a 'mids' object from mice")
  }
  
  m <- length(data_list)
  
  if (verbose) {
    cat("Starting prepost resilience analysis on pre-imputed datasets...\n")
    cat(sprintf("Number of imputed datasets: %d\n", m))
    cat(sprintf("Formula: %s\n", deparse(formula)))
    cat(sprintf("Resilience parameter (k): %.2f\n", k))
    cat(sprintf("Bootstrap samples per dataset: %d\n", nboot))
  }
  
  # Extract variables from formula
  formula_vars <- all.vars(formula)
  
  # Check if variables exist in the first dataset
  missing_vars <- setdiff(formula_vars, names(data_list[[1]]))
  if (length(missing_vars) > 0) {
    stop(sprintf("Variables not found in imputed datasets: %s", 
                 paste(missing_vars, collapse = ", ")))
  }
  
  # Step 1: Run prepost() on each imputed dataset
  if (verbose) cat("Running prepost analysis on imputed datasets...\n")
  
  results_list <- list()
  successful_imputations <- 0
  dataset_characteristics <- list()
  
  for (i in 1:m) {
    if (verbose) cat(sprintf("  Dataset %d/%d\n", i, m))
    
    current_data <- data_list[[i]]
    
    dataset_characteristics[[i]] <- list(
      n_obs = nrow(current_data),
      n_vars = ncol(current_data),
      missing_count = sum(is.na(current_data[, formula_vars])),
      complete_cases = sum(complete.cases(current_data[, formula_vars]))
    )
    
    tryCatch({
      result <- prepost(formula, 
                        data = current_data, 
                        k = k, 
                        nboot = nboot,
                        ...)
      
       if (!is.null(result$corrected.beta) && nrow(result$corrected.beta) > 0) {
        # Determine which k value to use (takes the first one)
        ci_name <- names(result$CI)[1]
        ci_matrix <- result$CI[[ci_name]]
        
        results_list[[i]] <- list(
          imputation = i,
          corrected_beta = result$corrected.beta,
          CI = ci_matrix,
          k_value = ci_name,
          dataset_summary = dataset_characteristics[[i]]
        )
        successful_imputations <- successful_imputations + 1
      }
    }, error = function(e) {
      warning(sprintf("Error in imputed dataset %d: %s", i, e$message))
    })
  }
  
  results_list <- results_list[!sapply(results_list, is.null)]
  
  if (length(results_list) == 0) {
    stop("No valid results from any imputed dataset. Check your data and model.")
  }
  
  if (verbose) {
    cat(sprintf("Successfully processed %d of %d imputed datasets\n", 
                successful_imputations, m))
  }
  
  # Step 2: Pool results using Rubin's rules
  if (verbose) cat("Pooling results using Rubin's rules...\n")
  
  reference_data <- data_list[[1]]
  
  if (pool_method == "detailed") {
    pooled_results <- rubin_pool_list(results_list, reference_data)
  } else {
    pooled_results <- rubin_pool_simple_list(results_list)
  }
  
  # Step 3: Create results table
  results_table <- create_results_table_list(pooled_results, coef_labels)
  
  # Step 4: Calculate imputation diagnostics
  imputation_diagnostics <- calculate_imputation_diagnostics(results_list, reference_data)
  
  # Step 5: Return comprehensive results
  output <- list(
    pooled_results = pooled_results,
    results_table = results_table,
    individual_results = results_list,
    imputation_info = list(
      n_imputations = m,
      successful_imputations = successful_imputations,
      dataset_characteristics = dataset_characteristics,
      diagnostics = imputation_diagnostics,
      mids_object = if(exists("mids_object")) mids_object else NULL
    ),
    model_info = list(
      formula = formula,
      formula_vars = formula_vars,
      k = k,
      nboot = nboot,
      pool_method = pool_method
    ),
    call = match.call()
  )
  
  class(output) <- "prepost_resilience_mi_list"
  
  if (verbose) {
    cat("Analysis complete!\n")
    cat("\nPooled Results:\n")
    print(results_table[, 1:6])
    cat("\nImputation Diagnostics:\n")
    cat(sprintf("  Fraction of Missing Information (average): %.3f\n",
                mean(pooled_results$fraction_missing_info, na.rm = TRUE)))
    cat(sprintf("  Relative Increase in Variance (average): %.3f\n",
                mean((pooled_results$between_var / pooled_results$within_var), na.rm = TRUE)))
  }
  
  return(output)
}

#' Rubin's pooling for list of results
#' @keywords internal
#' @noRd
rubin_pool_list <- function(results_list, reference_data) {
  m <- length(results_list)
  
  # Get coefficient names from first result
  coef_names <- rownames(results_list[[1]]$corrected_beta)
  n_coef <- length(coef_names)
  
  Q <- matrix(NA, nrow = n_coef, ncol = m)  # Point estimates
  U <- matrix(NA, nrow = n_coef, ncol = m)  # Within-imputation variances
  rownames(Q) <- rownames(U) <- coef_names
  
  for (i in 1:m) {
    Q[, i] <- results_list[[i]]$corrected_beta[, 1]
    
    # Calculate variances from CIs
    for (j in 1:n_coef) {
      coef_name <- coef_names[j]
      ci_data <- results_list[[i]]$CI
      
      if (coef_name %in% rownames(ci_data)) {
        ci_row <- ci_data[coef_name, ]
        
        # Handle different CI column names
        if ("upper" %in% names(ci_row) && "lower" %in% names(ci_row)) {
          ci_width <- ci_row["upper"] - ci_row["lower"]
        } else if ("97.5%" %in% names(ci_row) && "2.5%" %in% names(ci_row)) {
          ci_width <- ci_row["97.5%"] - ci_row["2.5%"]
        } else if (length(ci_row) >= 2) {
          # Take last column minus first column as fallback
          ci_width <- ci_row[length(ci_row)] - ci_row[1]
        } else {
          ci_width <- NA
        }
        
        if (!is.na(ci_width)) {
          U[j, i] <- (ci_width / (2 * qnorm(0.975)))^2
        } else {
          U[j, i] <- NA
        }
      } else {
        U[j, i] <- NA
      }
    }
  }
  
  # Handle any NA values
  for (j in 1:n_coef) {
    if (any(is.na(U[j, ]))) {
      U[j, is.na(U[j, ])] <- mean(U[j, ], na.rm = TRUE)
    }
    if (any(is.na(Q[j, ]))) {
      Q[j, is.na(Q[j, ])] <- mean(Q[j, ], na.rm = TRUE)
    }
  }
  
  # Apply Rubin's rules
  Q_bar <- rowMeans(Q)
  U_bar <- rowMeans(U)
  B <- apply(Q, 1, var)
  
  T_total <- U_bar + (1 + 1/m) * B
  pooled_se <- sqrt(T_total)
  
  # Calculate statistics for each coefficient
  n_obs <- nrow(reference_data)
  ci_lower <- numeric(n_coef)
  ci_upper <- numeric(n_coef)
  df <- numeric(n_coef)
  lambda <- numeric(n_coef)
  
  names(ci_lower) <- names(ci_upper) <- names(df) <- names(lambda) <- coef_names
  
  for (j in 1:n_coef) {
    lambda[j] <- ((1 + 1/m) * B[j]) / T_total[j]
    
    # Degrees of freedom (Barnard & Rubin, 1999)
    df_complete <- n_obs - n_coef - 1
    
    if (lambda[j] > 0 && !is.na(lambda[j])) {
      df_old <- (m - 1) / lambda[j]^2
      df_obs <- ((df_complete + 1) / (df_complete + 3)) * df_complete * (1 - lambda[j])
      df[j] <- (df_old * df_obs) / (df_old + df_obs)
    } else {
      df[j] <- df_complete
    }
    
    df[j] <- max(df[j], 1)
    
    # 95% confidence interval
    t_val <- qt(0.975, df = df[j])
    ci_lower[j] <- Q_bar[j] - t_val * pooled_se[j]
    ci_upper[j] <- Q_bar[j] + t_val * pooled_se[j]
  }
  
  return(list(
    coefficients = Q_bar,
    se = pooled_se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    df = df,
    within_var = U_bar,
    between_var = B,
    total_var = T_total,
    fraction_missing_info = lambda,
    n_imputations = m
  ))
}

#' Simple Rubin's rules for list of results
#' @keywords internal
#' @noRd
rubin_pool_simple_list <- function(results_list) {
  m <- length(results_list)
  coef_names <- rownames(results_list[[1]]$corrected_beta)
  n_coef <- length(coef_names)
  
  Q <- matrix(NA, nrow = n_coef, ncol = m)
  SE <- matrix(NA, nrow = n_coef, ncol = m)
  rownames(Q) <- rownames(SE) <- coef_names
  
  for (i in 1:m) {
    Q[, i] <- results_list[[i]]$corrected_beta[, 1]
    
    for (j in 1:n_coef) {
      coef_name <- coef_names[j]
      ci_data <- results_list[[i]]$CI
      
      if (coef_name %in% rownames(ci_data)) {
        ci_row <- ci_data[coef_name, ]
        
        if ("upper" %in% names(ci_row) && "lower" %in% names(ci_row)) {
          ci_width <- ci_row["upper"] - ci_row["lower"]
        } else if ("97.5%" %in% names(ci_row) && "2.5%" %in% names(ci_row)) {
          ci_width <- ci_row["97.5%"] - ci_row["2.5%"]
        } else if (length(ci_row) >= 2) {
          ci_width <- ci_row[length(ci_row)] - ci_row[1]
        } else {
          ci_width <- NA
        }
        
        if (!is.na(ci_width)) {
          SE[j, i] <- ci_width / (2 * qnorm(0.975))
        } else {
          SE[j, i] <- NA
        }
      } else {
        SE[j, i] <- NA
      }
    }
  }
  
  # Handle NAs
  for (j in 1:n_coef) {
    if (any(is.na(SE[j, ]))) {
      SE[j, is.na(SE[j, ])] <- mean(SE[j, ], na.rm = TRUE)
    }
    if (any(is.na(Q[j, ]))) {
      Q[j, is.na(Q[j, ])] <- mean(Q[j, ], na.rm = TRUE)
    }
  }
  
  # Rubin's rules
  Q_bar <- rowMeans(Q)
  U <- SE^2
  U_bar <- rowMeans(U)
  B <- apply(Q, 1, var)
  T_total <- U_bar + (1 + 1/m) * B
  pooled_se <- sqrt(T_total)
  
  # Degrees of freedom (simpler formula)
  df_old <- (m - 1) * (1 + (U_bar / ((1 + 1/m) * B)))^2
  df_old[is.na(df_old) | is.infinite(df_old)] <- m - 1
  df_old <- pmax(df_old, 1)
  
  # Confidence intervals
  t_vals <- qt(0.975, df = df_old)
  ci_lower <- Q_bar - t_vals * pooled_se
  ci_upper <- Q_bar + t_vals * pooled_se
  
  lambda <- ((1 + 1/m) * B) / T_total
  
  return(list(
    coefficients = Q_bar,
    se = pooled_se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    df = df_old,
    within_var = U_bar,
    between_var = B,
    total_var = T_total,
    fraction_missing_info = lambda,
    n_imputations = m
  ))
}

#' Create results table for list-based results
#' @keywords internal
#' @noRd
create_results_table_list <- function(pooled_results, coef_labels = NULL) {
  
  results_table <- data.frame(
    Coefficient = names(pooled_results$coefficients),
    Estimate = round(pooled_results$coefficients, 4),
    SE = round(pooled_results$se, 4),
    `95% CI Lower` = round(pooled_results$ci_lower, 4),
    `95% CI Upper` = round(pooled_results$ci_upper, 4),
    df = round(pooled_results$df, 1),
    Within_Var = round(pooled_results$within_var, 6),
    Between_Var = round(pooled_results$between_var, 6),
    FMI = round(pooled_results$fraction_missing_info, 3),
    stringsAsFactors = FALSE
  )
  
  # Apply custom coefficient labels if provided
  if (!is.null(coef_labels)) {
    if (!is.character(coef_labels) || is.null(names(coef_labels))) {
      warning("coef_labels should be a named character vector. Skipping label application.")
    } else {
      results_table$Coefficient <- ifelse(
        results_table$Coefficient %in% names(coef_labels),
        coef_labels[results_table$Coefficient],
        results_table$Coefficient
      )
    }
  }
  
  rownames(results_table) <- NULL
  return(results_table)
}

#' Calculate imputation diagnostics
#' @keywords internal
#' @noRd
calculate_imputation_diagnostics <- function(results_list, reference_data) {
  m <- length(results_list)
  n_coef <- length(rownames(results_list[[1]]$corrected_beta))
  
  # Extract estimates from all imputations
  Q <- matrix(NA, nrow = n_coef, ncol = m)
  coef_names <- rownames(results_list[[1]]$corrected_beta)
  rownames(Q) <- coef_names
  
  for (i in 1:m) {
    Q[, i] <- results_list[[i]]$corrected_beta[, 1]
  }
  
  # Calculate diagnostics
  diagnostics <- list(
    convergence = list(
      # Check for convergence in estimates across imputations
      max_relative_change = apply(Q, 1, function(x) {
        if (all(is.na(x))) return(NA)
        # Relative change from first to last imputation
        if (!is.na(x[1]) && x[1] != 0) {
          abs((x[m] - x[1]) / x[1])
        } else {
          NA
        }
      }),
      # Coefficient of variation across imputations
      cv_across_imputations = apply(Q, 1, function(x) {
        if (all(is.na(x))) return(NA)
        sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE))
      })
    ),
    between_imputation_variation = list(
      # Proportion of coefficients with substantial between-imputation variation
      # (CV > 0.2 is often considered substantial)
      high_variation_coefs = sum(apply(Q, 1, function(x) {
        if (all(is.na(x))) return(FALSE)
        cv <- sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE))
        return(cv > 0.2)
      }), na.rm = TRUE) / n_coef
    )
  )
  
  return(diagnostics)
}

#' S3 print method for prepost_resilience_mi_list
#' @title Print Results from Multiply Imputed Resilience Analysis
#' @param x An object of class \code{prepost_resilience_mi_list}.
#' @param ... Additional arguments passed to plotting functions.
#' @export
print.prepost_resilience_mi_list <- function(x, ...) {
  cat("Prepost Resilience Analysis on Pre-Imputed Datasets\n")
  cat("===================================================\n")
  cat(sprintf("Formula: %s\n", deparse(x$model_info$formula)))
  cat(sprintf("Imputed datasets: %d (successful: %d)\n", 
              x$imputation_info$n_imputations,
              x$imputation_info$successful_imputations))
  cat(sprintf("Resilience parameter (k): %.2f\n", x$model_info$k))
  cat(sprintf("Pooling method: %s\n", x$model_info$pool_method))
  cat("\nPooled Results:\n")
  print(x$results_table[, 1:6])
  invisible(x)
}

#' S3 summary method for prepost_resilience_mi_list
#' @title Summary Results from Multiply Imputed Resilience Analysis
#' @param object An object of class \code{prepost_resilience_mi_list}.
#' @param ... Additional arguments passed to plotting functions.
#' @export
summary.prepost_resilience_mi_list <- function(object, ...) {
  cat("Summary: Prepost Resilience Analysis on Pre-Imputed Datasets\n")
  cat("============================================================\n")
  
  cat("Model Information:\n")
  cat(sprintf("  Formula: %s\n", deparse(object$model_info$formula)))
  cat(sprintf("  Variables: %s\n", paste(object$model_info$formula_vars, collapse = ", ")))
  
  cat("\nImputation Information:\n")
  cat(sprintf("  Total imputations: %d\n", object$imputation_info$n_imputations))
  cat(sprintf("  Successful imputations: %d\n", object$imputation_info$successful_imputations))
  
  cat("\nPooling Diagnostics:\n")
  cat(sprintf("  Average Fraction of Missing Information: %.3f\n",
              mean(object$pooled_results$fraction_missing_info, na.rm = TRUE)))
  cat(sprintf("  Minimum FMI: %.3f\n",
              min(object$pooled_results$fraction_missing_info, na.rm = TRUE)))
  cat(sprintf("  Maximum FMI: %.3f\n",
              max(object$pooled_results$fraction_missing_info, na.rm = TRUE)))
  
  cat("\nKey Results:\n")
  print(object$results_table[, 1:5])
  
  invisible(object)
}

#' Plot method for visualizing imputation diagnostics
#' @title Plot Results from Multiply Imputed Resilience Analysis
#' @param x An object of class \code{prepost_resilience_mi_list}.
#' @param type Type of plot: \code{"estimates"} (default) or \code{"fmi"}.
#' @param ... Additional arguments passed to plotting functions.
#' @export
plot.prepost_resilience_mi_list <- function(x, type = "coefficients", ...) {
  if (type == "coefficients") {
    # Plot coefficient estimates across imputations
    n_coef <- length(x$pooled_results$coefficients)
    coef_names <- names(x$pooled_results$coefficients)
    n_imp <- x$imputation_info$n_imputations
    
    # Extract estimates into a matrix: rows = imputations, cols = coefficients
    est_matrix <- matrix(NA, nrow = n_imp, ncol = n_coef)
    rownames(est_matrix) <- paste("Imp", 1:n_imp)
    colnames(est_matrix) <- coef_names
    
    for (i in 1:min(n_imp, length(x$individual_results))) {
      est <- x$individual_results[[i]]$corrected_beta[, 1]
      est_matrix[i, names(est)] <- est
    }
    
    # Plot using base R
    matplot(
      x = 1:n_imp, 
      y = est_matrix, 
      type = "b", 
      pch = 19, 
      lty = 1, 
      col = 1:n_coef,
      xlab = "Imputation Number",
      ylab = "Coefficient Estimate",
      main = "Coefficient Estimates Across Imputations",
      ...
    )
    legend(
      "bottom", 
      legend = coef_names, 
      col = 1:n_coef, 
      lty = 1, 
      pch = 19, 
      ncol = 2, 
      cex = 0.8
    )
    
  } else if (type == "fmi") {
    # Plot fraction of missing information
    fmi_vals <- x$pooled_results$fraction_missing_info
    coef_names <- names(fmi_vals)
    
    # Order by FMI (descending)
    ord <- order(fmi_vals, decreasing = TRUE)
    fmi_vals <- fmi_vals[ord]
    coef_names <- coef_names[ord]
    
    # Horizontal bar plot
    par(mar = c(5, 8, 4, 2))  # more space for labels
    barplot(
      fmi_vals, 
      names.arg = coef_names, 
      horiz = TRUE, 
      col = "steelblue", 
      border = NA,
      xlab = "Fraction of Missing Information",
      main = "Fraction of Missing Information by Coefficient",
      ...
    )
    abline(v = 0.5, lty = 2, col = "red")
    mtext("High FMI (>0.5)", side = 3, line = -1, adj = 0.9, col = "red")
    
    par(mar = c(5, 4, 4, 2))  # reset margins
  }
  
  invisible(x)
}
utils::globalVariables(c("Imputation", "Estimate", "Coefficient", "FMI"))
########################################################################################
