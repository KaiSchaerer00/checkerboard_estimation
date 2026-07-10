library(copula)
library(dplyr)
library(ggplot2)

library(pbapply)
library(parallel)

source("regression.R", echo = TRUE)

K_copula <- function(XY, copula) {
  x = XY[,1]
  y = XY[,2]
  
  if (class(copula) == "claytonCopula") {
    theta = copula@parameters
    
    if (theta == 0) {
      result = y
    } else {
      condition = x^(-theta) + y^(-theta) > 1
      value_true = x^(-theta-1) * (x^(-theta) + y^(-theta) - 1)^(-1/theta - 1)
      value_false = 0
      result = ifelse(condition, value_true, value_false)
    }
  } else if (class(copula) == "amhCopula") {
    theta = copula@parameters
    result = (theta * (y-1) * y + y) / (1 - theta * (1-x) * (1-y))^2
  } else if (class(copula) == "moCopula") {
    alpha = copula@parameters[1]
    beta = copula@parameters[2]
    
    condition = (x^alpha > y^beta)
    value_true = (1 - alpha) * y * (x^(-alpha))
    value_false = y^(1 - alpha)
    result = ifelse(condition, value_true, value_false)
  } else if (class(copula) == "indepCopula") {
    result = y
  } else if (class(copula) == "frankCopula") {
    theta = copula@parameters 
    result <- exp( - theta * x ) * ( exp( - theta * y) - 1 ) / ( exp( - theta ) - 1 + ( exp( - theta * x) - 1 ) * ( exp( - theta * y) - 1 ) )
  }
  return(result)
}

###############################################################################
### Main ######################################################################
###############################################################################

### Setup Copula
theta <- 0.75
copula <- amhCopula(param=theta, dim=2)

# theta <- 2
# copula <- claytonCopula(param=theta, dim=2)

# theta <- c(0.5, 0.5)
# copula <- moCopula(theta)

# theta <- -2
# copula <- frankCopula(theta)

copula_name <- ""
if (class(copula) == "claytonCopula") {
  copula_name <- "Clayton"
} else if (class(copula) == "amhCopula") {
  copula_name <- "Ali-Mikhail-Haq"
} else if (class(copula) == "MOCopula") {
  copula_name <- "Marshall-Olkin"
} else if (class(copula) == "frankCopula") {
  copula_name <- "Frank"
}

n_values <- list(100, 250, 1000, 2500, 10000, 25000)
s <- 0.45

set.seed(0)
num_simulations <- 2000

all_results <- list()

cl <- makeCluster(parallel::detectCores() - 1)
clusterExport(cl, c("CB", "s", "ECBC", "rCopula", "copula", "K_copula"))

for (n in n_values) {
  N <- floor(n^s)
  m <- 2 * N^2
  
  clusterExport(cl, c("n", "N", "m"))
  
  results_vector <- pbsapply(1:num_simulations, function(i) {
    copula_samples <- rCopula(n, copula) 
    x_test <- matrix(runif(m * 2), m, 2)
    
    K_copula <- K_copula(x_test, copula)
    CB_N <- CB$new(x = copula_samples[,1], y = copula_samples[,2], s = s)
    K_CB_emp_cop <- CB_N$K_CB(x_test[,1], x_test[,2])
    
    max(abs(K_copula - K_CB_emp_cop))
  }, cl = cl)
  
  all_results[[as.character(n)]] <- results_vector
  print(median(results_vector))
}

library(ggplot2)
library(patchwork)

verylightgray <- rgb(0.94, 0.94, 0.94)

plot_data <- stack(all_results) %>%
  rename(error = values, n = ind) %>%
  mutate(n = factor(n, levels = as.character(n_values), 
                    labels = paste0("n=", n_values)))

ggplot(plot_data, aes(x = n, y = error)) +
  geom_boxplot(fill = "orange", color = "black", alpha=0.5) +
  labs(
    y = "Max Error",
    x = ""
  ) +
  theme_bw() +
  theme(
    axis.text=element_text(size=18), 
    axis.title=element_text(size=20),
    
    legend.position = "none",

  )

plotname = sprintf("%s\\%s_%s_convergence.pdf", class(copula), class(copula), gsub('\\.', '', paste(copula@parameters, collapse = "_")))
ggsave(plotname,
       width =12,
       height = 4, 
       units = "in")

save_name <- sprintf("%s\\%s_%s_convergence_all_results.rds", class(copula), class(copula), gsub('\\.', '', paste(copula@parameters, collapse = "_")))
saveRDS(all_results, save_name)

all_results <- readRDS(save_name)

