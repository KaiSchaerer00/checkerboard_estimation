source("../regression.R", echo = TRUE)
library(dplyr)
library(nw)
library(ggplot2)

shape <- function(x) {
  return(pmax(0.5, sqrt(x)))
}

scale <- function(x) {
  return(pmin(pmax(1, x), 6))
}

set.seed(0)

n_list <- c(100, 250, 1000, 2500, 10000, 25000)
s <- 0.45
quantile_level <- 0.5
sim_rep <- 2000
m <- 2000

C <- 100

alpha_param <- 1
beta_param <- 1

df_cb_max  <- data.frame()
df_cb_mean <- data.frame()
df_alt_max  <- data.frame()
df_alt_mean <- data.frame()

for ( n in n_list ) {
  N <- floor(n^s)
  
  x_mat <- matrix(
    rbeta(sim_rep * n, alpha_param, beta_param),
    nrow = n
  ) * 10
  
  shape_vec <- apply(x_mat, 2, shape)
  scale_vec <- apply(x_mat, 2, scale)
  
  alpha_vec <- shape_vec * (-shape_vec * scale_vec - scale_vec + C) / C
  beta_vec  <- (shape_vec * scale_vec - C) *
    (shape_vec * scale_vec + scale_vec - C) / (C * scale_vec)
  
  y_mat <- matrix(
    rbeta(sim_rep * n, alpha_vec, beta_vec),
    nrow = n
  ) * C
  
  x_list_mat <- matrix(
    rbeta(sim_rep * m, alpha_param, beta_param),
    nrow = m
  ) * 10
  
  shape_list_mat <- apply(x_list_mat, 2, shape)
  scale_list_mat <- apply(x_list_mat, 2, scale)
  
  results_vector <- matrix(0, 4, sim_rep) 
  
  
  
  for (i in seq_len(sim_rep)) {
    if (i %% 200 == 0) {
      print(i)
    }
    
    x <- x_mat[, i]
    y <- y_mat[, i]
    
    ord_x_list <- order(x_list_mat[, i])
    
    x_list <- x_list_mat[, i][order(x_list_mat[, i])]
    shape_list <- shape_list_mat[, i][order(x_list_mat[, i])]
    scale_list <- scale_list_mat[, i][order(x_list_mat[, i])]
    
    CB_N <- CB$new(x = x, y = y, s = s, quantile_level = quantile_level)
    
    q_CB_list <- CB_N$q_CB(x_list)
    alpha_list <- shape_list*(-shape_list*scale_list - scale_list + C)/C
    beta_list <- (shape_list*scale_list - C)*(shape_list*scale_list + scale_list - C)/(C*scale_list)
    q_true <- qbeta(quantile_level, alpha_list, beta_list) * C
    q_alt <- nw_q(x, y, quantile_level, kernel = "normal", bandwidth = sd(x) * n ^(-1/5), x.points = x_list)$y
    
    cb_diff <- abs(q_CB_list - q_true)
    alt_diff <- abs(q_alt - q_true)
    
    results_vector[,i] <- c(max(cb_diff), mean(cb_diff), max(alt_diff), mean(alt_diff))
  }
  
  df_cb_max  <- df_cb_max  |> rbind(data.frame(x = results_vector[1,], n = n))
  df_cb_mean <- df_cb_mean |> rbind(data.frame(x = results_vector[2,], n = n))
  df_alt_max  <- df_alt_max  |> rbind(data.frame(x = results_vector[3,], n = n))
  df_alt_mean <- df_alt_mean |> rbind(data.frame(x = results_vector[4,], n = n))
  print(c(median(results_vector[1,]),
          median(results_vector[2,]),
          median(results_vector[3,]),
          median(results_vector[4,])))
  
}



for (method in c("max", "mean")) {
  if (method == "max") {
    df_cb <- df_cb_max
    df_alt <- df_alt_max
    y_title <- "Max Error"
  } else {
    df_cb <- df_cb_mean
    df_alt <- df_alt_mean
    y_title <- "Mean Error"
  }
  
  df_cb$type <- "Checkerboard Quantile Regression"
  df_alt$type <- "Nadaraya-Watson Quantile Regression"
  
  df <- df_cb |> rbind(df_alt)
  
  df$type <- factor(df$type, levels = c("Checkerboard Quantile Regression", "Nadaraya-Watson Quantile Regression"))
  
  ggplot(df) + geom_boxplot(aes(x=factor(n),y=x,colour=factor(type), fill=factor(type)), alpha = 0.25) + 
    labs(
         x="n", 
         y=y_title, 
         colour ="Estimation", fill = "Estimation") + 
    scale_colour_manual(values = c("Checkerboard Quantile Regression" = "blue", "Nadaraya-Watson Quantile Regression" ="orange" )) +
    scale_fill_manual(values = c("Checkerboard Quantile Regression" = "blue", "Nadaraya-Watson Quantile Regression" ="orange" )) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          
          axis.text=element_text(size=12), 
          axis.title=element_text(size=14),
          
          legend.position = c(0.8, 0.8),
          legend.box.background = element_rect(colour = "black"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
  
  save_name <- sprintf("quantile_regression_convergence_%s.pdf", method)
  ggsave(save_name, width = 16, height = 4, dpi = 400)
}

saveRDS(df_cb_max, file = "quantile_regression_convergence_df_cb_max.rds")
saveRDS(df_cb_mean, file = "quantile_regression_convergence_df_cb_mean.rds")
saveRDS(df_alt_max, file = "quantile_regression_convergence_df_nw_max.rds")
saveRDS(df_alt_mean, file = "quantile_regression_convergence_df_nw_mean.rds")

df_cb_max <- readRDS("quantile_regression_convergence_df_cb_max.rds")
df_cb_mean <- readRDS("quantile_regression_convergence_df_cb_mean.rds")
df_alt_max <- readRDS("quantile_regression_convergence_df_nw_max.rds")
df_alt_mean <- readRDS("quantile_regression_convergence_df_nw_mean.rds")
