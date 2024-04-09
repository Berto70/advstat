source('./project/code/func_bb.R')
library(ggplot2)

# data
set.seed(42)
dt <- 0.05
t <- dt * seq(0,999)
x <- rep(0, length(t))
x[sample(seq_along(t), length(t) %/% 10)] <- 1
df <- data.frame(t, x)
edges.uniform <- bayesian_blocks(t=t, x=x, data_type='Events', dt=dt, p0=0.01)
cat("edges dataframe:\n")
print(edges.uniform)

N.uniform <- length(t)
N.edges.uniform <- length(edges.uniform)
steps.uniform <- rep(0, length(N.edges.uniform))

    for (K in seq(1: N.edges.uniform)){

        steps.uniform[K] <- sum(x[t>=edges.uniform[K] & t<edges.uniform[K+1]])/sum(x*(edges.uniform[K+1] -edges.uniform[K]))

        if (K==N.edges.uniform-1){break}
    }
steps.uniform[N.edges.uniform] <- steps.uniform[N.edges.uniform-1]

cat("steps dataframe:\n")
print(steps.uniform)

p.unif <- ggplot()+
geom_histogram(aes(df$t, y=after_stat(density), color=as.factor(1), fill=as.factor(1)), alpha=0.5, bins=31)+
geom_step(aes(edges.uniform, steps.uniform, color=as.factor(2), fill=as.factor(2)), linewidth=1)
print(p.unif)


# t <- rnorm(5000)
# df <- data.frame(t)
# edges <- bayesian_blocks(t, data_type="Events")
# cat("edges dataframe:\n")
# print(edges)

# t <- 100 * runif(100)
# x <- exp(-0.5*(t-50)^2)
# sigma <- 100
# x_obs <- rnorm(x, sigma)
# edges <- bayesian_blocks(t=t, x=x_obs, sigma=sigma, data_type='PointMeasures', p0=0.01)
#
# print(edges)

# dt <- 0.05
# t <- dt * rnorm(1000)
# x <- rep(0, length(t))
# x[sample(1:length(t), length(t) %/% 10)] <- 1
# edges <- bayesian_blocks(t, x, data_type='RegularEvents', dt=dt)
# print(edges)