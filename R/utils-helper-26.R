

#' Generate Age Measurement Schedules for Individuals
#'
#' This function generates a list of age vectors for a given number of
#' individuals, ensuring each individual has measurements that:
#' - Start at a specified minimum age
#' - End at a specified maximum age
#' - Include at least one measurement in a specified target age
#' range (e.g., pubertal window)
#' - Have a total number of measurements sampled from a defined range, 
#' with an approximate target median
#' - Enforce minimum counts of individuals with the lowest and highest 
#' possible measurement counts
#' 
#' @param n Integer. Number of individuals.
#' @param meas_range Integer vector. Range of total measurement counts per
#'   individual (e.g., `5:15`).
#' @param min_age Numeric. Minimum age (start of measurement range, e.g., `6`).
#' @param max_age Numeric. Maximum age (end of measurement range, e.g., `20`).
#' @param median_target Integer. Desired median number of measurements (must be
#'   in `meas_range`).
#' @param min_n_at_min Integer. Minimum number of individuals with the lowest
#'   measurement count (i.e., `min(meas_range)`).
#' @param min_n_at_max Integer. Minimum number of individuals with the highest
#'   measurement count (i.e., `max(meas_range)`).
#' @param must_include_range Numeric vector of length 2. A continuous age window
#'   (e.g., `c(11, 15)`) where each individual must have at least one
#'   measurement.
#' 
#' @return A list with two elements:
#' \describe{
#'   \item{ages_list}{A list of numeric vectors, one per individual, containing
#'   age values sorted in ascending order. Each vector starts at `min_age`, ends
#'   at `max_age`, and includes at least one value in `must_include_range`.}
#'   \item{num_measurements}{An integer vector of length `n`, showing the total
#'   number of measurements per individual.}
#' }
#' 
#' @examples
#' set.seed(2025)
#' res <- generate_age_list(
#'   n = 50,
#'   meas_range = 5:15,
#'   min_age = 6,
#'   max_age = 20,
#'   median_target = 10,
#'   min_n_at_min = 3,
#'   min_n_at_max = 4,
#'   must_include_range = c(11, 15)
#' )
#' str(res$ages_list[[1]])  # View age schedule for one individual
#' table(res$num_measurements)  # Check distribution
#'
#' @keywords internal
#' @noRd
generate_age_list <- function(n,
                              meas_range = 5:15,
                              min_age = 6,
                              max_age = 20,
                              median_target = 10,
                              min_n_at_min = 3,
                              min_n_at_max = 4,
                              must_include_range = c(11, 15),
                              seed = NULL) {
  
  if(is.null(seed)) {
    set.seed(123)
  } else {
    set.seed(seed)
  }
  
  
  stopifnot(length(meas_range) >= 2)
  stopifnot(median_target %in% meas_range)
  stopifnot(min_n_at_min + min_n_at_max <= n)
  
  # Triangular distribution centered at median
  x <- meas_range
  center <- which(x == median_target)
  probs <- sapply(1:length(x), function(i) 1 - abs(i - center) / max(center - 1,
                                                                     length(x) - 
                                                                       center))
  probs <- probs / sum(probs)
  
  # Forced assignments
  min_val <- min(meas_range)
  max_val <- max(meas_range)
  fixed_values <- c(rep(min_val, min_n_at_min), rep(max_val, min_n_at_max))
  remaining_n <- n - length(fixed_values)
  remaining_vals <- sample(meas_range, remaining_n, replace = TRUE, 
                           prob = probs)
  num_meas_per_id <- sample(c(fixed_values, remaining_vals))  # shuffle
  
  # Function to generate individual age vector
  generate_ages <- function(k) {
    if (k < 3) stop2c("Need at least 3 measurements to satisfy 
                      constraints (start, middle, end)")
    
    # Always include start, middle (in must_include_range), and end
    middle_age <- runif(1, min = must_include_range[1], 
                        max = must_include_range[2])
    
    remaining_k <- k - 3
    other_ages <- if (remaining_k > 0) {
      runif(remaining_k, min = min_age, max = max_age)
    } else {
      numeric(0)
    }
    
    sort(c(min_age, max_age, middle_age, other_ages))
  }
  
  # Generate list of age vectors
  ages_list <- lapply(num_meas_per_id, generate_ages)
  
  list(ages_list = ages_list,
       num_measurements = num_meas_per_id)
}


###############################################################################
###############################################################################

# This below code 'generate_age_list_new' gives control over min/max age gap 
# but fails sometimes

# Helper: Create random gaps under constraints
create_gaps_with_constraints <- function(n_gaps, 
                                         total_range, 
                                         min_gap, 
                                         max_gap, 
                                         max_iter = 500) {
  for (i in 1:max_iter) {
    weights <- runif(n_gaps)
    gaps <- weights / sum(weights) * total_range
    if (all(gaps >= min_gap) && all(gaps <= max_gap)) {
      return(gaps)
    }
  }
  return(NULL)
}

# Helper: Generate age vector with constraints and flexible fallback
generate_ages <- function(k, 
                          min_age, 
                          max_age, 
                          must_include_range, 
                          min_gap, 
                          max_gap, 
                          relax = TRUE) {
  total_range <- max_age - min_age
  n_gaps <- k - 1
  
  attempt <- 1
  max_attempts <- 100
  relax_step <- 0.01
  
  mgap <- min_gap
  xgap <- max_gap
  
  while (attempt <= max_attempts) {
    gaps <- create_gaps_with_constraints(n_gaps, total_range, mgap, xgap)
    if (!is.null(gaps)) {
      age_seq <- cumsum(c(min_age, gaps))
      if (!any(age_seq >= must_include_range[1] & 
               age_seq <= must_include_range[2])) {
        replace_idx <- sample(2:(k - 1), 1)
        age_seq[replace_idx] <- runif(1, must_include_range[1],
                                      must_include_range[2])
        age_seq <- sort(age_seq)
      }
      return(age_seq)
    }
    if (relax) {
      mgap <- max(mgap - relax_step, 0.01)
      xgap <- xgap + relax_step
    }
    attempt <- attempt + 1
  }
  return(NULL)
}

# Main: Generate age list
generate_age_list_new <- function(n,
                              meas_range = 5:15,
                              min_age = 6,
                              max_age = 20,
                              median_target = 10,
                              min_n_at_min = 3,
                              min_n_at_max = 4,
                              must_include_range = c(11, 15),
                              impose_agegap = TRUE,
                              min_gap = 0.1,
                              max_gap = 1.0,
                              seed = 123) {
  
  set.seed(seed)
  stopifnot(length(meas_range) >= 2)
  stopifnot(median_target %in% meas_range)
  stopifnot(min_n_at_min + min_n_at_max <= n)
  
  # Triangle-shaped probabilities
  x <- meas_range
  center <- which(x == median_target)
  probs <- sapply(1:length(x), 
                  function(i) 1 - abs(i - center) / 
                    max(center - 1, length(x) - center))
  probs <- probs / sum(probs)
  
  min_val <- min(meas_range)
  max_val <- max(meas_range)
  fixed_values <- c(rep(min_val, min_n_at_min), rep(max_val, min_n_at_max))
  remaining_n <- n - length(fixed_values)
  remaining_vals <- sample(meas_range, remaining_n, replace = TRUE, 
                           prob = probs)
  num_meas_per_id <- sample(c(fixed_values, remaining_vals))
  
  ages_list <- vector("list", length = n)
  for (i in seq_len(n)) {
    k <- num_meas_per_id[i]
    age_vec <- NULL
    
    if (impose_agegap) {
      age_vec <- generate_ages(k, 
                               min_age, 
                               max_age,
                               must_include_range, 
                               min_gap,
                               max_gap, 
                               relax = TRUE)
    } else {
      mid_age <- runif(1, must_include_range[1], must_include_range[2])
      other_ages <- if (k > 3) runif(k - 3, min_age, max_age) else numeric(0)
      age_vec <- sort(c(min_age, max_age, mid_age, other_ages))
    }
    
    if (is.null(age_vec)) {
      warning2c(sprintf("Individual %d: Failed to generate age
                        vector. Relax gap constraints or age range.", i))
      stop2c("Could not generate valid age vector even after fallback.")
    }
    
    ages_list[[i]] <- age_vec
  }
  
  list(ages_list = ages_list, num_measurements = num_meas_per_id)
}






# # ---- Test Run ----
# result <- generate_age_list(
#   n = 50,
#   meas_range = 5:15,
#   min_age = 6,
#   max_age = 20,
#   median_target = 10,
#   min_n_at_min = 3,
#   min_n_at_max = 4,
#   must_include_range = c(11, 15),
#   impose_agegap = TRUE,
#   min_gap = 0.1,
#   max_gap = 4.0,
#   seed = 42
# )
# 
# # ---- Check Gaps ----
# ages_list <- result$ages_list
# 
# 
# gap_summary <- map_dfr(seq_along(ages_list), function(i) {
#   age_vec <- sort(ages_list[[i]])
#   gaps <- diff(age_vec)
#   data.frame(
#     id = factor(i),
#     min_gap = min(gaps),
#     max_gap = max(gaps),
#     mean_gap = mean(gaps)
#   )
# })
# 
# gap_summary %>% dplyr::reframe(mean_min_gap = mean(min_gap),
#                                mean_mean_gap = mean(mean_gap),
#                                mean_max_gap = mean(max_gap)
# ) %>% 
#   mutate(across(where(is.numeric), function(x) round(x, 2)))
# 
# # ---- Caterpillar Plot ----
# x_min <- max(min(gap_summary$min_gap) - 0.5, 0)
# x_max <- max(gap_summary$max_gap) + 0.5
# 
# ggplot(gap_summary, aes(y = id)) +
#   geom_linerange(aes(xmin = min_gap, xmax = max_gap), color = "gray60") +
#   geom_point(aes(x = min_gap), color = "red", size = 2) +
#   geom_point(aes(x = mean_gap), color = "blue", size = 2) +
#   geom_point(aes(x = max_gap), color = "darkgreen", size = 2) +
#   geom_text(aes(x = min_gap, label = sprintf("%.2f", min_gap)),
#             vjust = -1, nudge_x = 0.05, size = 2.8, color = "red") +
#   geom_text(aes(x = mean_gap, label = sprintf("%.2f", mean_gap)),
#             vjust = -1, nudge_x = 0.05, size = 2.8, color = "blue") +
#   geom_text(aes(x = max_gap, label = sprintf("%.2f", max_gap)),
#             vjust = -1, nudge_x = 0.05, size = 2.8, color = "darkgreen") +
#   labs(title = "Caterpillar Plot of Age Gaps per Individual",
#        x = "Age Gap (years)", y = "Individual ID") +
#   scale_x_continuous(expand = c(0.05, 0.05)) +
#   coord_cartesian(xlim = c(x_min, x_max)) +
#   theme_minimal() +
#   theme(axis.text.y = element_text(size = 6))


