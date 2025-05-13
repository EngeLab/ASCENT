# Function to generate genotypes and expected BAF
model_gts <- function(cn, err) {
  # Make genotypes (A/B)
  nB <- 0:cn
  genotypes <- sapply(nB, function(b) {
    paste0(paste(rep("A", cn - b), collapse = ""),
           paste(rep("B", b), collapse = ""))
  })
  # Expected BAFs
  true_proportion <- nB / cn
  expected_BAF <- true_proportion * (1 - 2*err) + err
  return(list(genotypes = genotypes, expected_BAF = expected_BAF))
}

# Function to calculate dynamic error rate based on total count
calc_dynamic_error <- function(total_count, base_error, 
                               min_count = 10, max_count = 500,
                               max_multiplier = 5) {
  if (total_count <= min_count) {
    return(base_error * max_multiplier)
  }
  if (total_count >= max_count) {
    return(base_error)
  }
  scale_factor <- 1 + (max_multiplier - 1) * 
    exp(-(total_count - min_count) / (max_count - min_count) * 3)
  return(base_error * scale_factor)
}

# Function to calculate genotype probabilities
calc_gt_prob <- function(A_count, B_count, cn, err, 
                         dynamic_error = FALSE,
                         min_count = 10, max_count = 500, 
                         max_multiplier = 5) {
  
  total <- A_count + B_count
  
  # Optional: read-depth adjusted error rate
  if (dynamic_error) {
    err <- calc_dynamic_error(total, err, min_count, max_count, max_multiplier)
  }
  
  model <- model_gts(cn, err)
  # Calculate log likelihoods using binomial distribution
  log_likelihoods <- dbinom(B_count, total, model$expected_BAF, log = TRUE)
  
  # Convert to probabilities (normalized)
  max_log_lik <- max(log_likelihoods)
  probs <- exp(log_likelihoods - max_log_lik)
  probs <- probs / sum(probs)
  
  # Return results
  return(data.frame(
    gt = model$genotypes,
    expected_BAF = model$expected_BAF,
    log_likelihood = log_likelihoods,
    probability = probs,
    error_rate = err
  ))
}

# Function to get most likely genotype
get_gt <- function(probs_df) {
  result <- probs_df[order(-probs_df$probability),]
  
  # Return the best genotype and its probability
  return(list(
    best_genotype = as.character(result$gt[1]),
    probability = result$probability[1],
    log_likelihood = result$log_likelihood[1],
    error_rate = result$error_rate[1],
    results = result
  ))
}

# Create heatmap function
create_heatmap <- function(cn, base_err = 0.02, dynamic = FALSE,
                           max_count = 500, min_count = 10, max_multiplier = 5,
                           log = TRUE) {
  # Create a finer grid with fixed steps (like the older version)
  baf_values <- seq(0, 1, by = 0.02)  # More fine-grained BAF steps
  
  # Use appropriate spacing for count values based on log parameter
  if (log) {
    # Create logarithmically spaced points, with more points at lower values
    log_min <- log(min_count)
    log_max <- log(max_count)
    # Use more points for better resolution
    total_counts <- exp(seq(log_min, log_max, length.out = 150))
  } else {
    # Use linear spacing with fixed step size
    total_counts <- seq(min_count, max_count, by = (max_count - min_count) / 80)
  }
  
  # Ensure total_counts don't exceed max_count
  total_counts <- pmin(total_counts, max_count)
  
  # Create a complete grid
  grid <- expand.grid(BAF = baf_values, Total = total_counts)
  grid$Genotype <- NA_character_
  grid$Probability <- NA_real_
  grid$Uncertainty <- NA_real_  # Add uncertainty measure from older version
  
  # Generate genotype names
  model <- model_gts(cn, base_err)
  genotypes <- model$genotypes
  
  # Process the grid with improved boundary handling and better error handling
  for (i in 1:nrow(grid)) {
    baf <- grid$BAF[i]
    total <- grid$Total[i]
    
    # Calculate counts
    B_count <- round(total * baf)
    A_count <- total - B_count
    
    # Calculate error rate (fixed or dynamic)
    err <- base_err
    if (dynamic) {
      err <- calc_dynamic_error(total, base_err, min_count, max_count, max_multiplier)
    }
    
    # Calculate expected BAFs directly
    expected_BAF <- (0:cn)/cn * (1-2*err) + err
    
    # Ensure values are within bounds - use slightly wider bounds
    expected_BAF <- pmax(pmin(expected_BAF, 1-1e-14), 1e-14)
    
    # Calculate log-likelihoods with proper error handling
    log_liks <- rep(-Inf, length(expected_BAF))
    
    # Use tryCatch to suppress warnings and handle errors
    for (j in 1:length(expected_BAF)) {
      log_liks[j] <- tryCatch({
        dbinom(B_count, total, expected_BAF[j], log = TRUE)
      }, warning = function(w) {
        return(-Inf)
      }, error = function(e) {
        return(-Inf)
      })
    }
    
    # Additional check for numerical issues
    log_liks[is.nan(log_liks) | is.na(log_liks) | is.infinite(log_liks)] <- -Inf
    
    # Skip if all values are -Inf
    if (all(is.infinite(log_liks))) {
      next
    }
    
    # Calculate probabilities with improved numerical stability
    max_ll <- max(log_liks[is.finite(log_liks)])
    probs <- exp(log_liks - max_ll)
    probs <- probs / sum(probs, na.rm = TRUE)
    
    # Find max likelihood genotype
    best_idx <- which.max(probs)
    
    # Calculate uncertainty (1 - difference between top two probabilities)
    # This is key for better boundary representation
    sorted_probs <- sort(probs, decreasing = TRUE)
    uncertainty <- 1
    if (length(sorted_probs) >= 2) {
      uncertainty <- 1 - (sorted_probs[1] - sorted_probs[2])
    }
    
    # Store results
    grid$Genotype[i] <- genotypes[best_idx]
    grid$Probability[i] <- probs[best_idx]
    grid$Uncertainty[i] <- uncertainty
  }
  
  # Remove rows with NA values
  grid <- grid[!is.na(grid$Genotype), ]
  
  # Define colors based on copy number - using more pastel colors like older version
  colors <- switch(as.character(cn),
                   "2" = c("#6495ED", "#90EE90", "#FA8072"), # blue, light green, salmon
                   "3" = c("#6495ED", "#00CED1", "#90EE90", "#FA8072"), # blue, cyan, light green, salmon
                   "4" = c("#6495ED", "#00CED1", "#90EE90", "#FFFFE0", "#FFA07A", "#EE82EE"), # blue, cyan, light green, light yellow, light salmon, violet
                   "5" = c("#6495ED", "#00CED1", "#90EE90", "#FFFFE0", "#FFA07A", "#EE82EE"),
                   colorRampPalette(c("#6495ED", "#00CED1", "#90EE90", "#FFFFE0", "#FFA07A", "#EE82EE"))(cn + 1))
  
  # Create color map - order from all B to all A
  color_map <- setNames(rev(colors[1:(cn+1)]), genotypes)
  
  # Create plot with improved aesthetics similar to older version
  title <- if (dynamic) {
    paste0("CN = ", cn, " with dynamic error rate\nbase error rate = ", 
           base_err, " , max multiplier = ", max_multiplier)
  } else {
    paste0("CN = ", cn, " with fixed error rate (", base_err*100, "%)")
  }
  
  # Create the plot with adjustments for better visualization
  p <- ggplot(grid, aes(x = BAF, y = Total, fill = Genotype)) +
    # Use uncertainty for alpha to create smoother transitions at boundaries
    geom_tile(aes(alpha = 1 - (0.4 * Uncertainty))) +  # Adjust alpha based on uncertainty
    scale_fill_manual(values = color_map, name = "Most likely\ngenotype") +
    # Adjust alpha range for better visualization
    scale_alpha_continuous(range = c(0.5, 1), name = "Probability") +
    labs(
      title = title,
      x = "B Allele Frequency (BAF)",
      y = "Total Allele Count"
    ) +
    theme_minimal() +
    # Apply log scale to y-axis only if requested
    {if (log) scale_y_log10() else scale_y_continuous()} +
    geom_vline(xintercept = (0:cn)/cn, linetype = "dashed", color = "gray30") +
    geom_vline(xintercept = c(0.1,0.9), linetype = "dashed", color = "pink") +
    coord_cartesian(expand = FALSE) +
    # Improve theme for better visualization
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", size = 0.2),
      plot.title = element_text(size = 11)
    )
  
  return(p)
}

# library(tidyverse)
# library(ggplot2)
library(gridExtra)
# # Parameters
base_err <- 0.02
min_count <- 5
max_count <- 500
max_multiplier <- 5
copy_numbers <- c(2, 3, 4, 5)

# # Create all plots - fixed error rate
# fixed_plots_list <- list()
# for (i in 1:length(copy_numbers)) {
#   cn <- copy_numbers[i]
#   fixed_plots_list[[i]] <- create_heatmap(cn=cn, base_err=base_err, dynamic = FALSE,
#                                           max_count = max_count, min_count = min_count,
#                                           max_multiplier = max_multiplier, log=F)
# }
# # Arrange fixed error plots
# grid.arrange(grobs = fixed_plots_list, ncol = 2, nrow = 2,
#              top = paste0("Fixed Error Rate Models (",base_err*100,"%)"))


# Create all plots - dynamic error rate
dynamic_plots_list <- list()
for (i in 1:length(copy_numbers)) {
  cn <- copy_numbers[i]
  dynamic_plots_list[[i]] <- create_heatmap(cn=cn, base_err=base_err, dynamic = TRUE,
                                            max_count = max_count, min_count = min_count,
                                            max_multiplier = max_multiplier, log=F)
}


# Arrange dynamic error plots
pdf("/wrk/data/solrun/AnalysisPipelinePaper/ALL40/Figure4/Supplementary.pdf", width=15, height=12)
grid.arrange(grobs = dynamic_plots_list, ncol = 2, nrow = 2)
dev.off()
