library(ggplot2)
source('./project/code/func_bb.R')

set.seed(42)

# Define the function with fluctuations and plateaus
# f <- function(x) {
#   y <- -x + 1
#   y <- y + rnorm(length(x), mean = 0.5, sd = 100)  # Add Gaussian noise with mean 0.5 and standard deviation 100
#   plateau_x_values <- c(0.2, 0.5, 0.9)
#   plateau_heights <- c(0.1, 1.5, 0.2)
#   for (i in seq_along(plateau_x_values)) {
#     x_val <- plateau_x_values[i]
#     height <- plateau_heights[i]
#     y[abs(x - x_val) < 0.05] <- y[abs(x - x_val) < 0.05] + height
#   }
#   return(y)
# }

# Generate random numbers in the interval [0, 1]
# num_samples <- 10000
# samples <- runif(num_samples)
#
# # Apply the modified function to the samples
# sampled_values <- f(samples)
#
# # Define the number of bins
# num_bins <- as.integer(sqrt(num_samples)/2)
# bins <- seq(0, 1, length.out = num_bins)
#
# dataframe <- read.csv('./project/code/data_bbh.csv')
# val <- na.omit(dataframe$val)
# hbins <- na.omit(dataframe$hbins)
# data <- dataframe$data
#
# x <- rep(0, length(data))
# x[sample(seq_along(data), length(data) %/% 10)] <- 1
#
# bb_edges <- bayesian_blocks(data, data_type='Events')


# Define parameters
# dt <- 0.05
# data <- dt * seq(1, 1000)
# x <- rep(0, length(data))
# x[sample(seq_along(data), length(data)%/%10)] <- 1

data <- 100 * runif(100)
x <- exp(-0.5 * (data - 50) ^ 2)
sigma <- 0.1
x_obs <- rnorm(length(x), mean = x, sd = sigma)

# Perform Bayesian blocks analysis
bb_edges <- bayesian_blocks(data, x_obs, data_type="PointMeasures", sigma=sigma)


bb_num_bins <- length(bb_edges)
bb_steps <- rep(0, length(bb_num_bins))

for (K in seq(1:bb_num_bins-1)){
  bb_steps[K] <- sum(x[data>=bb_edges[K] & data<bb_edges[K+1]])/sum(x*(bb_edges[K+1] - bb_edges[K]))
}

bb_steps[bb_num_bins] <- bb_steps[bb_num_bins-1]

ggplot() +
    geom_histogram(aes(x=data, y=after_stat(density), color=as.factor(1), fill=as.factor(1)), alpha=0.5, bins=num_bins, position='identity') +
    geom_step(aes(bb_edges, bb_steps, color=as.factor(2), fill=as.factor(2)), linewidth=1) +
    # geom_vline(xintercept=bb_edges)+
    # scale_x_reverse() +
     geom_histogram(aes(x=data, y=after_stat(density), color=as.factor(3)), alpha=0.3, breaks=bb_edges, position='identity')+
    # labs(x='q', y='Counts', title='Distribution of BBHs q')+
    theme(legend.position.inside = c(0.8, 0.8))

