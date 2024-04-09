# sink(stdout(), type = "message")
# options(error = recover)
options(warn = -1)

# Install and load the R6 package
if (!require(R6)) {
    install.packages("R6")
}
library(R6) # object-oriented programming in R

if (!require(tidyverse)) {
    install.packages("tidyverse")

}
library(tidyverse)

bayesian_blocks <- function(data, gamma=NULL, dt = NULL, fitness='events') {
    FITNESS_DICT <- list(
        "events" = Events$new(data, fitness),
        "regularevents" = RegularEvents$new(data, dt, fitness),
        'pointmeasures' = PointMeasures$new(data, fitness)
    )

    fitfunc <- FITNESS_DICT[[fitness]]

    if (is.null(fitfunc)) {
        stop("Invalid fitness type")
    }

    return(fitfunc$fit(data, fitness))
}

FitnessFunc <- R6Class(
    "FitnessFunc",

    public = list(
        p0 = NULL,
        gamma = NULL,
        ncp_prior = NULL,

        initialize = function(p0 = 0.05, gamma = NULL, ncp_prior = NULL) {
            self$p0 <- p0
            self$gamma <- gamma
            self$ncp_prior <- ncp_prior
        },

        validate_input = function(data, fitness='events') {
            cells <- data
            if(!is.data.frame(cells)){
                cells <- as.data.frame(cells)
            }

            if (fitness=='events'){
                cells <- cells |> count(cells)
                return(cells[order(cells[1], )])
            }

            if (fitness=='regularevents'){
                
                return(cells[order(cells[1], )])
            }

            if (fitness=='pointmeasures'){
                
                return(cells[order(cells[1], )])
            }
        },

        # fitness = function(...) {
        #     stop("Not implemented")
        # },

        p0_prior = function(N) {
            return(4 - log(73.53 * self$p0 * (N ** -0.478)))
        },

        compute_ncp_prior = function(N) {
            if (!is.null(self$gamma)) {
                return(-log(self$gamma))
            } else if (!is.null(self$p0)) {
                return(self$p0_prior(N))
            } else {
                stop("ncp_prior cannot be computed as neither gamma nor p0 is defined")
            }
        },

        fit = function(data, fitness, gamma=NULL) {

            cells <- data
            if(!is.data.frame(cells)){
                cells <- as.data.frame(cells)
            }

            if (fitness=='events'){
                cells <- cells |> count(cells)
                return(cells[order(cells[1], )])
            }

            if (fitness=='regularevents'){
                
                return(cells[order(cells[1]), ])
            }

            if (fitness=='pointmeasures'){
                
                return(cells[order(cells[1], )])
            }
            
            t <- cells[[1]]
            x <- cells[[2]]
            if (length(data)==3){
                sigma <- cells[[3]]
            } else {
                sigma <- 1
            }

            if ("a_k" %in% names(formals(self$fitness))) {
                ak_raw <- rep(1, times = length(x)) / sigma ** 2
            }
            if ("b_k" %in% names(formals(self$fitness))) {
                bk_raw <- x / sigma ** 2
            }

            N <- length(t)
            edges <- c(t[1], 0.5 * (t[2:N] + t[1:(N - 1)]), t[N])
            block_length <- t[N] - edges

            best <- rep(0, times = N)
            last <- rep(1, times = N)

            if (is.null(self$ncp_prior)) {
                ncp_prior <- self$compute_ncp_prior(N)
            } else {
                ncp_prior <- self$ncp_prior
            }

            for (K in seq(1:N)) {
                kwds <- list()

                if ("T_k" %in% names(formals(self$fitness))) {
                    kwds$T_k <- block_length[1:K] - block_length[K + 1]
                }

                if ("N_k" %in% names(formals(self$fitness))) {
                    kwds$N_k <- rev(cumsum(rev(x[1:K])))
                }

                if ("a_k" %in% names(formals(self$fitness))) {
                    kwds$a_k <- 0.5 * rev(cumsum(rev(ak_raw[1:K])))
                }

                if ("b_k" %in% names(formals(self$fitness))) {
                    kwds$b_k <- -1 * rev(cumsum(rev(bk_raw[1:K])))
                }

                fit_vec <- do.call(self$fitness, kwds)

                A_R <- fit_vec - ncp_prior
                A_R[2:length(A_R)] <- A_R[2:length(A_R)] + best[1:(K - 1)]

                i_max <- which.max(A_R)
                last[K] <- i_max
                best[K] <- A_R[i_max]
            }

            change_points <- rep(1, times = N)
            i_cp <- N
            ind <- N + 1
            while (i_cp > 1) {
                change_points[i_cp] <- ind
                if (ind == 1) { break }
                i_cp <- i_cp - 1
                ind <- last[ind - 1]
            }

            if (i_cp == 1) {
                change_points[i_cp] <- 1
            }
            change_points <- change_points[i_cp:length(change_points)]

            return(edges[change_points])
        }
    )
)

Events <- R6Class(
    "Events",
    inherit = FitnessFunc,

    public = list(

        initialize = function(p0 = 0.05, gamma = NULL, ncp_prior = NULL) {
            super$initialize(p0, gamma, ncp_prior)
        },

        fitness = function(N_k, T_k) {
            return(N_k * (log(N_k / T_k)))
        },

        validate_input = function(t, x, sigma) {
            validation <- super$validate_input(t, x, sigma)
            t <- validation[[1]]
            x <- validation[[2]]
            sigma <- validation[[3]]
            if (!is.null(x) && any(x %% 1 > 0)) {
                stop("x must be integer counts for fitness='events'")
            }
            return(list(t, x, sigma))
        }
    )
)

RegularEvents <- R6Class(
    "RegularEvents",
    inherit = FitnessFunc,

    public = list(
        dt = NULL,

        initialize = function(dt = dt, p0 = 0.05, gamma = NULL, ncp_prior = NULL) {
            self$dt <- dt
            super$initialize(p0, gamma, ncp_prior)
        },

        validate_input = function(data) {
            cells <- super$validate_input(data)
            t <- cells[[1]]
            x <- cells[[2]]
            if (length(data)==3){
                sigma <- cells[[3]]
            } else {
                sigma <- 1
            }

            if (!all(x %in% c(0, 1))) {
                stop("Regular events must have only 0 and 1 in x")
            }
            return(cells[order(cells[1], )])
        },

        fitness = function(T_k, N_k) {
            M_k <- T_k / self$dt
            N_over_M <- N_k / M_k

            eps <- 1e-8
            if (any(N_over_M > 1 + eps)) {
                warning("regular events: N/M > 1.  Is the time step correct?")
            }

            one_m_NM <- 1 - N_over_M
            N_over_M[N_over_M <= 0] <- 1
            one_m_NM[one_m_NM <= 0] <- 1

            return(N_k * log(N_over_M) + (M_k - N_k) * log(one_m_NM))
        }
    )
)

PointMeasures <- R6Class(
    "PointMeasures",
    inherit = FitnessFunc,

    public = list(

        initialize = function(p0 = 0.05, gamma = NULL, ncp_prior = NULL) {
            super$initialize(p0, gamma, ncp_prior)
        },

        fitness = function(a_k, b_k){
            return((b_k * b_k) / (4 * a_k))
        },

        validate_input = function(t, x, sigma) {
            if (is.null(x)){
                stop("x must be specified for fitness=`PointMeasures`")
            }
            validation <- super$validate_input(t, x, sigma)
            t <- validation[[1]]
            x <- validation[[2]]
            sigma <- validation[[3]]

            return(list(t, x, sigma))
        }
    )
)