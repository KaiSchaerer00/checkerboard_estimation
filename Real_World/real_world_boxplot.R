source("../regression.R", echo = TRUE)
library(nw)
library(copula)


data(loss)

set.seed(0)

filtered_loss <- loss[loss$censored != 1, ]
X <- log(filtered_loss$loss)
Y <- log(filtered_loss$alae)

n <- length(X)
s <- 0.45
m = 10000

res_vec <- sapply(1:m, function(i) {
  
  sample <- sample(c(TRUE, FALSE), n, 
                   replace=TRUE, prob=c(0.8,0.2))
  
  X_train <- X[sample]
  Y_train <- Y[sample]
  X_test <- X[!sample]
  Y_test <- Y[!sample]
  
  CB_N <- CB$new(x = X_train, y = Y_train, s = s)
  
  r_CB <- CB_N$r_CB(X_test)

  r_NW <- nw_e(X_train, Y_train, kernel = "normal", bandwidth = sd(X_test) * n ^(-1/5), x.points = X_test)$y
  
  loss_CB <- abs(r_CB - Y_test)
  loss_NW <- abs(r_NW - Y_test)

  return (c(max(loss_CB), mean(loss_CB), max(loss_NW), mean(loss_NW)))
})

max_loss_CB <- list(res_vec[1,])
mean_loss_CB <- list(res_vec[2,])
max_loss_NW <- list(res_vec[3,])
mean_loss_NW <- list(res_vec[4,])


df_cb_max = data.frame(x = res_vec[1,])
df_cb_max$type <- "Checkerboard\nMean Regression"
df_nw_max = data.frame(x = res_vec[3,])
df_nw_max$type <- "Nadaraya-Watson\nRegression"
df_max <- df_cb_max |> rbind(df_nw_max)

ggplot(df_max) + geom_boxplot(aes(x = factor(n), y=x,colour=factor(type), fill=factor(type)), alpha = 0.25) + 
  labs(
       x="n", 
       y="Max Error", 
       colour ="Estimation", fill = "Estimation") + 
  scale_colour_manual(values = c("Checkerboard\nMean Regression" = "blue", "Nadaraya-Watson\nRegression" ="orange" )) +
  scale_fill_manual(values = c("Checkerboard\nMean Regression" = "blue", "Nadaraya-Watson\nRegression" ="orange" )) +
  theme_minimal() + 
  theme(
        
        axis.text=element_text(size=18), 
        axis.title=element_text(size=20),
        
        plot.margin=unit(c(1,7,1,0.5),"cm"),
        
        legend.position = c(1.2, 0.8),
        legend.box.background = element_rect(colour = "black"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20))

save_name = "real_world_mean_boxplot_max.pdf"
ggsave(save_name, width = 8, height = 6, dpi = 400)

saveRDS(df_max, file = "real_world_mean_df_max")
df_max <- readRDS("real_world_mean_df_max")


df_cb_mean = data.frame(x = res_vec[2,])
df_cb_mean$type <- "Checkerboard\nMean Regression"
df_nw_mean = data.frame(x = res_vec[4,])
df_nw_mean$type <- "Nadaraya-Watson\nRegression"
df_mean <- df_cb_mean |> rbind(df_nw_mean)

ggplot(df_mean) + geom_boxplot(aes(x = factor(n), y=x,colour=factor(type), fill=factor(type)), alpha = 0.25) + 
  labs(
       x="n", 
       y="Mean Error", 
       colour ="Estimation", fill = "Estimation") + 
  scale_colour_manual(values = c("Checkerboard\nMean Regression" = "blue", "Nadaraya-Watson\nRegression" ="orange" )) +
  scale_fill_manual(values = c("Checkerboard\nMean Regression" = "blue", "Nadaraya-Watson\nRegression" ="orange" )) +
  theme_minimal() + 
  theme(
        
        axis.text=element_text(size=18), 
        axis.title=element_text(size=20),
        
        plot.margin=unit(c(1,7,1,0.5),"cm"),
        
        legend.position = c(1.2, 0.8),
        legend.box.background = element_rect(colour = "black"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20))


save_name = "real_world_mean_boxplot_mean.pdf"
ggsave(save_name, width = 8, height = 6, dpi = 400)

saveRDS(df_mean, file = "real_world_quant_df_mean")
df_mean <- readRDS("real_world_quant_df_mean")
