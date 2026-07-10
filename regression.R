library(stats)
library(qad)
library(copula)
library(ggplot2)
library(R6)

CB <- R6Class("CB",
                private = list(
                n = NULL,
                N = NULL,
                
                F_n = NULL,
                G_n = NULL,
                
                prob_mass_matrix = NULL,
                K_CB_matrix = NULL,
                r_CB_matrix = NULL,
                q_CB_matrix = NULL,
                var_CB_matrix = NULL,
                
                calculate_r_CB_matrix = function() {
                  N <- private$N
                  n <- private$n
                  
                  X <- (1:N) / N
                  
                  y_pos <- sort(self$y[which(self$y > 0)])
                  diff_y_pos <- diff(c(0, y_pos))
                  
                  M_summand <- matrix(
                    data = 1 - self$K_CB(rep(X , each = n), 
                                         rep(0:(n-1) / n , times = N)), 
                    nrow = N, 
                    ncol = n, 
                    byrow = TRUE
                  )
                  
                  M_product <- M_summand * matrix(rep(diff_y_pos, each = N), N, n)
                  results_vector <- rowSums(M_product)
                  
                  return(results_vector)
                },
                
                calculate_2M_CB_matrix = function() {
                  N <- private$N
                  n <- private$n
                  
                  X <- (1:N) / N
                  
                  y_pos <- sort(self$y[which(self$y > 0)])
                  diff_y_pos <- diff(c(0, y_pos^2))
                  M_summand <- matrix(
                    data = 1 - self$K_CB(rep(X , each = n), 
                                         rep(0:(n-1) / n , times = N)), 
                    nrow = N, 
                    ncol = n, 
                    byrow = TRUE
                  )
                  
                  M_product <- M_summand * matrix(rep(diff_y_pos, each = N), N, n)
                  results_vector <- rowSums(M_product)
                  
                  return(results_vector)
                },
                
                calculate_q_CB_matrix = function(quantile_level) {
                  N <- private$N
                  n <- private$n
                  
                  X <- (1:N) / N
                  row_nums <- 1:N
                  
                  K_CB_M <- matrix(
                    data = self$K_CB(rep(X , each = n+1),
                                     rep(0:n / n , times = N)), 
                    nrow = N, 
                    ncol = n+1, 
                    byrow = TRUE
                  )
                  
                  logical_matrix <- K_CB_M >= quantile_level
                  indices <- max.col(logical_matrix, ties.method = "first")
                  
                  interpolating_low <- K_CB_M[cbind(row_nums, indices-1)]
                  interpolating_high <- K_CB_M[cbind(row_nums, indices)]
                  
                  interpolating_factor <- (quantile_level - interpolating_low) / (interpolating_high - interpolating_low)
                  
                  interpolated_val <- ((indices - 1) + interpolating_factor)/n
                  
                  results_vector <- unname(quantile(private$G_n, interpolated_val))
                  
                  return(results_vector)
                }
              ),
              
              public = list(
                x = NULL,
                y = NULL,
                s = NULL,
                quantile_level = NULL,
                
                K_CB = function(X, Y) {
                  X = pmin(pmax(X, 0), 1)
                  Y = pmin(pmax(Y, 0), 1)
                  
                  i_0 = ceiling(X * private$N)
                  j_0 = ceiling(Y * private$N)
                  
                  i_0[i_0 == 0] <- 1
                  j_0[j_0 == 0] <- 1
                  
                  interpolate_low  <- private$K_CB_matrix[cbind(i_0, j_0)]
                  interpolate_high <- private$K_CB_matrix[cbind(i_0, j_0+1)]
                  interpolate_val  <- Y * private$N - (j_0 - 1)
                  result = interpolate_low + (interpolate_high - interpolate_low) * interpolate_val
                  
                  return(result)
                },
                
                initialize = function(x, y, s, quantile_level = NULL, variance = FALSE) {
                  
                  self$x <- x
                  self$y <- y
                  self$s <- s
                  self$quantile_level <- quantile_level
                  
                  private$n <- length(self$x)
                  private$N <- floor(private$n^s)
                  
                  private$F_n <- ecdf(self$x)
                  private$G_n <- ecdf(self$y)
                  
                  private$prob_mass_matrix <- ECBC(X = private$F_n(self$x), 
                                                   Y = private$G_n(self$y), 
                                                   resolution = private$N)
                  
                  private$K_CB_matrix <- matrix(0, nrow = private$N, ncol = private$N + 1)
                  private$K_CB_matrix[, 2:(private$N + 1)] <- t(apply(private$prob_mass_matrix, 1, cumsum)) * private$N
                  
                  private$r_CB_matrix <- private$calculate_r_CB_matrix()
                  
                  if (!is.null(quantile_level)) {
                    private$q_CB_matrix <- private$calculate_q_CB_matrix(quantile_level)
                  }
                  
                  if (variance == TRUE) {
                    private$var_CB_matrix <- private$calculate_var_CB_matrix()
                  }
                  
                },
                
                r_CB = function(X) {
                  i_0 <- pmax(1, ceiling(private$F_n(X) * private$N))
                  return (private$r_CB_matrix[i_0])
                },
                
                q = function(quantile_level) {
                  private$q_CB_matrix <- private$calculate_q_CB_matrix(quantile_level)
                },
                
                q_CB = function(X) {
                  i_0 = ceiling(private$F_n(X) * private$N)
                  i_0[i_0 == 0] <- 1
                  return (private$q_CB_matrix[i_0])
                },
                
                var = function(){
                  private$var_CB_matrix <- private$calculate_2M_CB_matrix() - private$r_CB_matrix^2
                },
                
                var_CB = function(X) {
                  i_0 = ceiling(private$F_n(X) * private$N)
                  i_0[i_0 == 0] <- 1
                  return (private$var_CB_matrix[i_0])
                }
              )
)




