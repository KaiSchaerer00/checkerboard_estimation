library(copula)
library(qad)
library(nw)
library(ggplot2)
source("../regression.R", echo = TRUE)


plot = function(method, x, y, x_list, estim_val, alt_val, save_name) {
  if (method == "E") {
    line_labels <- c("Training Data", 
                     "Checkerboard Mean Regression",
                     "Nadaraya-Watson Regression")
    title <- sprintf("Mean Regression (n = %s)", length(x))

  } else {
    line_labels <- c("Training Data",
                     "Checkerboard Quantile Regression",
                     "Nadaraya Watson Quantile Regression")
    title <- sprintf("Quantile Regression (n = %s)", length(x))
  }
  
  manual_colours <- c("grey", "blue", "orange")
  names(manual_colours) <- line_labels 
  
  plot_data_wide <- data.frame(
    Estimator_1 = estim_val, 
    Estimator_2 = alt_val
  )
  
  ggplot() +
    geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y, color = line_labels[1]),
               alpha = 0.2, size = 2) +
    
    geom_line(data = plot_data_wide, aes(x = x_list, y = estim_val, color = line_labels[2]),
              linewidth = 1.2) +
    
    geom_line(data = plot_data_wide, aes(x = x_list, y = alt_val, color = line_labels[3]),
              linewidth = 1.2) +
    
    coord_cartesian(xlim = c(min(x), max(x)), ylim = c(min(y), max(y))) +
    
    labs(
         x = "log-idemnity payment",
         y = "log-ALAE",
         color = "Legend") +
    
    scale_color_manual(breaks = line_labels,
                       values = manual_colours) + 
    
    theme_minimal() + 
    theme(
          
          legend.position = c(0.85, 0.25),
          legend.box.background = element_rect(colour = "black"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          
          axis.text=element_text(size=18), 
          axis.title=element_text(size=20),
          
          panel.grid.major = element_line(color = "gray80", linetype = "dashed", linewidth = 0.5),
          panel.grid.minor = element_line(color = "gray90", linetype = "dotted", linewidth = 0.2),
          panel.background = element_rect(fill = 'white', colour = 'black'))
  
  
  
  ggsave(save_name, width = 16, height = 4, dpi = 400)
  
}



data(loss) 
filtered_loss <- loss[loss$censored != 1, ]
X <- log(filtered_loss$loss)
Y <- log(filtered_loss$alae)

n = length(X)

CB_N <- CB$new(x = X, y = Y, s = 0.45)

ggplot(data.frame(x = X, y = Y), aes(x=X, y=Y)) +
  geom_point(shape=1)

x_list <- seq(min(X), max(X), by=0.05)



# Expectation
r_CB_list <- CB_N$r_CB(x_list)
r_alt <- nw_e(X, Y, kernel = "normal", bandwidth = sd(X) * n ^(-1/5), x.points = x_list)$y

save_name = "C:/Users/b1119245/00_Kai/03_Code/20251014_uniform_convergence_K_CB/regression/real_world_mean_reg.pdf"
plot("E", X, Y, x_list, r_CB_list, r_alt, save_name)

# Quantiles
quantile_level <- 0.5
CB_N$q(quantile_level)

q_CB_list <- CB_N$q_CB(x_list)
q_alt <- nw_q(X, Y, quantile_level, kernel = "normal", bandwidth = sd(X) * n ^(-1/5), x.points = x_list)$y

save_name = "real_world_quant_reg.pdf"
plot("Q", X, Y, x_list, q_CB_list, q_alt, save_name)

# VaR
quantile_level <- 0.9
CB_N$q(quantile_level)

VaR_CB_list <- CB_N$q_CB(x_list)

VaR_alt <- nw_q(X, Y, quantile_level, kernel = "normal", bandwidth = sd(X) * n ^(-1/5), x.points = x_list)$y

save_name = "real_world_VaR_reg.pdf"
plot("Q", X, Y, x_list, VaR_CB_list, VaR_alt, save_name)

