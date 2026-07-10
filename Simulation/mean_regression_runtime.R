source("../regression.R", echo = TRUE)
library(dplyr)
library(nw)
library(ggplot2)

nw_r <- function(x_train, y_train, x_eval, h) {
  
  h = h * 0.3706506
  gaussian_kernel <- function(u) (1/sqrt(2*pi)) * exp(-0.5 * u^2)
  
  sapply(x_eval, function(x) {
    weights <- gaussian_kernel((x - x_train) / h)
    
    if (sum(weights) == 0) return(NA)
    return(sum(weights * y_train) / sum(weights))
  })
}

shape <- function(x) {
  return(pmax(0.5, sqrt(x)))
}

scale <- function(x) {
  return(pmin(pmax(1, x), 6))
}

convergence <- function(n_list, m, s) {
  
  
}

shape <- function(x) {
  return(pmax(0.5, sqrt(x)))
}

scale <- function(x) {
  return(pmin(pmax(1, x), 6))
}

set.seed(0)

n_list <- c(100, 250, 1000, 2500, 10000, 25000)
s <- 1/3
sim_rep <- 2000
m <- 2000

alpha_param <- 1
beta_param <- 1

df_cb_time = c()
df_nw_time = c()

for ( n in n_list ) {
  N <- floor(n^s)
  m <- floor(n^(1/2))
  
  print("")
  print(c("n = ", n))
  print(c("N = ", N))
  print(c("m = ", m))
  
  x_mat <- matrix(
    rbeta(sim_rep * n, alpha_param, beta_param),
    nrow = n
  ) * 10
  
  shape_vec <- apply(x_mat, 2, shape)
  scale_vec <- apply(x_mat, 2, scale)
  
  alpha_vec <- shape_vec * (-shape_vec * scale_vec - scale_vec + 100) / 100
  beta_vec  <- (shape_vec * scale_vec - 100) *
    (shape_vec * scale_vec + scale_vec - 100) / (100 * scale_vec)
  
  y_mat <- matrix(
    rbeta(sim_rep * n, alpha_vec, beta_vec),
    nrow = n
  ) * 100
  
  x_list_mat <- matrix(
    rbeta(sim_rep * m, alpha_param, beta_param),
    nrow = m
  ) * 10
  
  shape_list_mat <- apply(x_list_mat, 2, shape)
  scale_list_mat <- apply(x_list_mat, 2, scale)
  

  
  results_vector <- pbsapply(1:sim_rep, function(i) {
    x <- x_mat[, i]
    y <- y_mat[, i]
    
    ord_x_list <- order(x_list_mat[, i])
    
    x_list <- x_list_mat[, i][order(x_list_mat[, i])]
    shape_list <- shape_list_mat[, i][order(x_list_mat[, i])]
    scale_list <- scale_list_mat[, i][order(x_list_mat[, i])]
    
    r_true <- shape_list * scale_list
    
    start_cb <- Sys.time()
    CB_N <- CB$new(x = x, y = y, s = s)
    
    r_CB_list <- CB_N$r_CB(x_list)
    end_cb <- Sys.time()
    cb_time <- end_cb - start_cb
    
    start_nw <- Sys.time()
    r_alt <- nw_r(x, y, x_list, sd(x) * n ^(-1/5))
    end_nw <- Sys.time()
    nw_time <- end_nw - start_nw
    
    return (c(cb_time, nw_time))
  })
  
  df_cb_time  <- df_cb_time  |> rbind(data.frame(x = results_vector[1,], n = n))
  df_nw_time  <- df_nw_time  |> rbind(data.frame(x = results_vector[2,], n = n))
  print(c(median(results_vector[1,]),
          median(results_vector[2,])))
  
}


mean_df_cb_time <- df_cb_time %>% group_by(n) %>% summarise(cb_time = mean(x))
mean_df_nw_time <- df_nw_time %>% group_by(n) %>% summarise(nw_time = mean(x))

mean_df_time <- inner_join(mean_df_cb_time, mean_df_nw_time, by="n")

df_long <- mean_df_time %>%
  pivot_longer(
    cols = c(cb_time, nw_time), 
    names_to = "time_type", 
    values_to = "seconds"
  )

ggplot(df_long, aes(x=n, y=seconds, color=time_type)) + 
  geom_line(linewidth=1) + 
  geom_point(size=3) + 
  
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(
    values = c("cb_time" = "blue", "nw_time" = "orange"),
    labels = c("cb_time" = "Checkerboard Mean Regression", "nw_time" = "Nadaraya-Watson Regression") 
  ) +
  
  labs(
    x = "n",
    y = "Time (seconds)", 
    color = "Measurements"
  ) + 
  theme_bw()

save_name <- "mean_regression_runtime.pdf"
ggsave(save_name, width = 16, height = 6, dpi = 400)

saveRDS(df_cb_time, "mean_regression_runtime_cb.rds")
saveRDS(df_nw_time, "mean_regression_runtime_nw.rds")
