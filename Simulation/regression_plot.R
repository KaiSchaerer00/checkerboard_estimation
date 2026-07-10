library(nw)
source("../regression.R", echo = TRUE)

plot = function(method, x, y, x_list, true_val, estim_val, alt_val, save_name) {
  if (method == "E") {
    line_labels <- c("Training Data", 
                     "True Conditional Expectation",
                     "Checkerboard Mean Regression",
                     "Nadaraya-Watson Regression")
    title <- sprintf("Mean Regression (n = %s)", length(x))
    
    
  } else {
    line_labels <- c("Training Data",
                     "True Conditional Quantile",
                     "Checkerboard Quantile Regression",
                     "Naradaya-Watson Quantile Regression")
    title <- sprintf("Quantile Regression (n = %s)", length(x))
  }
  
  manual_colours <- c("grey", "red", "blue", "orange")
  names(manual_colours) <- line_labels 
  
  plot_data_wide <- data.frame(
    True_Value = true_val,
    Estimator_1 = estim_val, 
    Estimator_2 = alt_val
  )
  
  ggplot() +
    geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y, color = line_labels[1]),
               alpha = 0.2, size = 1) +
    
    geom_line(data = plot_data_wide, aes(x = x_list, y = true_val, color = line_labels[2]),
              linewidth = 1) +
    
    geom_step(data = plot_data_wide, aes(x = x_list, y = estim_val, color = line_labels[3]),
              linewidth = 1) +
    
    geom_line(data = plot_data_wide, aes(x = x_list, y = alt_val, color = line_labels[4]),
              linewidth = 1) +
    
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 35)) +
    
    labs(
         x = "X",
         y = "Y",
         color = "Legend") +
    
    scale_color_manual(breaks = line_labels,
                       values = manual_colours) + 
    
    theme_bw() + 
    theme(
          
          legend.position = c(0.2, 0.8),
          legend.box.background = element_rect(colour = "black"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          
          axis.text=element_text(size=18), 
          axis.title=element_text(size=20),
          
          panel.grid.major = element_line(color = "gray80", linetype = "dashed", linewidth = 0.5),
          panel.grid.minor = element_line(color = "gray90", linetype = "dotted", linewidth = 0.2),
          panel.background = element_rect(fill = 'white', colour = 'black'))
  
  
  
  ggsave(save_name, width = 16, height = 6, dpi = 400)
}

shape <- function(x) {
  return(pmax(1, sqrt(x)))
}

scale <- function(x) {
  return(pmin(pmax(1, x), 6))
}


set.seed(0)
n = 10000
s = 0.45

alpha_param <- 1
beta_param <- 1

C <- 100

x <- rbeta(n, alpha_param, beta_param) * 10
shape_vec <- shape(x)
scale_vec <- scale(x)
alpha_vec <- shape_vec*(-shape_vec*scale_vec - scale_vec + C)/C
beta_vec <- (shape_vec*scale_vec - C)*(shape_vec*scale_vec + scale_vec - C)/(C*scale_vec)
y <- rbeta(n, alpha_vec, beta_vec) * C

CB_N <- CB$new(x = x, y = y, s = s)

x_list <- seq(0, 10, by = 0.05)
shape_list <- shape(x_list)
scale_list <- scale(x_list)


# Expectation
r_CB_list <- CB_N$r_CB(x_list)
r_true <- shape_list * scale_list
r_alt <- nw_e(x, y, kernel = "normal", bandwidth = sd(x) * n ^(-1/5), x.points = x_list)$y

save_name <- sprintf("mean_regression_%06d.pdf", length(x))
plot("E", x, y, x_list, r_true, r_CB_list, r_alt, save_name)

# Quantiles
quantile_level <- 0.5
CB_N$q(quantile_level)

q_CB_list <- CB_N$q_CB(x_list)

alpha_list <- shape_list*(-shape_list*scale_list - scale_list + C)/C
beta_list <- (shape_list*scale_list - C)*(shape_list*scale_list + scale_list - C)/(C*scale_list)
q_true <- qbeta(quantile_level, alpha_list, beta_list) * C

q_alt <- nw_q(x, y, quantile_level, kernel = "normal", bandwidth = sd(x) * n ^(-1/5), x.points = x_list)$y

save_name <- sprintf("quantile_regression_%06d.pdf", length(x))
plot("Q", x, y, x_list, q_true, q_CB_list, q_alt, save_name)

# Variance
CB_N$var()
var_CB_list <- CB_N$var_CB(x_list)
var_true <- alpha_list * beta_list / ((alpha_list + beta_list)^2 * (alpha_list + beta_list + 1)) * 10000

line_labels <- c("True Conditional Variance",
                 "Checkerboard Variance Regression")

manual_colours <- c("red", "blue")
names(manual_colours) <- line_labels 


plot_data_wide <- data.frame(
  True_Value = var_true,
  Estimator_1 = var_CB_list
)

ggplot() +
  
  geom_line(data = plot_data_wide, aes(x = x_list, y = var_true, color = line_labels[1]),
            linewidth = 1) +
  
  geom_step(data = plot_data_wide, aes(x = x_list, y = var_CB_list, color = line_labels[2]),
            linewidth = 1) +
  
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 120)) +
  
  labs(
    x = "X",
    y = "Y",
    color = "Legend") +
  
  scale_color_manual(breaks = line_labels,
                     values = manual_colours) + 
  
  theme_bw() + 
  theme(
    
    legend.position = c(0.2, 0.8),
    legend.box.background = element_rect(colour = "black"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    
    axis.text=element_text(size=18), 
    axis.title=element_text(size=20),
    
    panel.grid.major = element_line(color = "gray80", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor = element_line(color = "gray90", linetype = "dotted", linewidth = 0.2),
    panel.background = element_rect(fill = 'white', colour = 'black'))


save_name <- sprintf("variance_regression_%06d.pdf", length(x))
ggsave(save_name, width = 16, height = 6, dpi = 400)

