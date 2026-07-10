source("../regression.R", echo = TRUE)
library(dplyr)
library(nw)
library(ggplot2)

shape <- function(x) {
  return(pmax(1, sqrt(x)) * (1 + sin(x*10) / 4))
}

scale <- function(x) {
  return(pmin(pmax(1, x), 6))
}

set.seed(0)

n_list <- c(100, 250, 1000, 2500, 10000, 25000, 100000)
s <- 0.45
m <- 2000

alpha_param <- 1
beta_param <- 1

df_cb_max  <- data.frame()
df_cb_mean <- data.frame()
df_alt_max  <- data.frame()
df_alt_mean <- data.frame()

for ( n in n_list ) {
  N <- floor(n^s)
  
  results_vector <- sapply(1:2000, function(i) {
    x <- rbeta(n, alpha_param, beta_param) * 10
    shape_vec <- shape(x)
    scale_vec <- scale(x)
    alpha_vec <- shape_vec*(-shape_vec*scale_vec - scale_vec + 100)/100
    beta_vec <- (shape_vec*scale_vec - 100)*(shape_vec*scale_vec + scale_vec - 100)/(100*scale_vec)
    y <- rbeta(n, alpha_vec, beta_vec) * 100
    
    CB_N <- CB$new(x = x, y = y, s = s)
    
    x_list <- sort(rbeta(m, 1, 1) * 10)
    shape_list <- shape(x_list)
    scale_list <- scale(x_list)
    
    r_CB_list <- CB_N$r_CB(x_list)
    r_true <- shape_list * scale_list
    r_alt <- nw_e(x, y, kernel = "normal", bandwidth = sd(x) * n ^(-1/5), x.points = x_list)$y
    
    cb_diff <- abs(r_CB_list - r_true)
    alt_diff <- abs(r_alt - r_true)
    
    return (c(max(cb_diff), mean(cb_diff), max(alt_diff), mean(alt_diff)))
    
  })
  
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
  
  df_cb$type <- "Checkerboard Mean Regression"
  df_alt$type <- "Nadaraya-Watson Regression"
  
  df <- df_cb |> rbind(df_alt)
  
  df$type <- factor(df$type, levels = c("Checkerboard Mean Regression", "Nadaraya-Watson Regression"))
  
  ggplot(df) + geom_boxplot(aes(x=factor(n),y=x,colour=factor(type), fill=factor(type)), alpha = 0.25) + 
    labs(
         x="n", 
         y=y_title, 
         colour ="Estimation", fill = "Estimation") + 
    scale_colour_manual(values = c("Checkerboard Mean Regression" = "blue", "Nadaraya-Watson Regression" ="orange" )) +
    scale_fill_manual(values = c("Checkerboard Mean Regression" = "blue", "Nadaraya-Watson Regression" ="orange" )) +
    theme_minimal() + 
    theme(
          axis.text=element_text(size=18), 
          axis.title=element_text(size=16),
          
          legend.position = c(0.8, 0.8),
          legend.box.background = element_rect(colour = "black"),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18))
  
  save_name <- sprintf("mean_regression_sin_convergence_%s.pdf", method)
  ggsave(save_name, width = 16, height = 4, dpi = 400)
}

saveRDS(df_cb_max, file = "mean_regression_sin_convergence_df_cb_max.rds")
saveRDS(df_cb_mean, file = "mean_regression_sin_convergence_df_cb_mean.rds")
saveRDS(df_alt_max, file = "mean_regression_sin_convergence_df_nw_max.rds")
saveRDS(df_alt_mean, file = "mean_regression_sin_convergence_df_nw_mean.rds")

df_cb_max <- readRDS("mean_regression_sin_convergence_df_cb_max.rds")
df_cb_mean <- readRDS("mean_regression_sin_convergence_df_cb_mean.rds")
df_alt_max <- readRDS("mean_regression_sin_convergence_df_nw_max.rds")
df_alt_mean <- readRDS("mean_regression_sin_convergence_df_nw_mean.rds")
