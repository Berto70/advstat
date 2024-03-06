bayesian_blocks <- function(t, x=NULL, p0=0.5, sigma=NULL, dt=NULL, gamma=NULL, ncp_prior=NULL, data_type='Events'){
    # Data validation
    validation <- validate_input(t, x=NULL, sigma=NULL, data_type)

    t <- validation$t
    x <- validation$x
    sigma <- validation$sigma

    edges <- c(t[1], 0.5 * (t[2:length(t)] + t[1:(length(t) - 1)]), t[length(t)])
    block_length <- t[length(t)] - edges

    N <- length(t)
    best <- rep(0, times = N)
    last <- rep(1, times = N)

    if (is.null(ncp_prior)) {
        ncp_prior <- compute_ncp_prior(N, p0, gamma)
    } else {
        ncp_prior <- ncp_prior
    }

    #-----------------------------------------------------------------
    # Start with first data cell; add one cell at each iteration
    #-----------------------------------------------------------------

    for (K in seq(1:N)) {
        kwds <- list()

        if (data_type=='Events' | data_type=='RegularEvents'){
            kwds$T_k <- block_length[1:K] - block_length[K + 1]
            kwds$N_k <- rev(cumsum(rev(x[1:K])))

            fit_vec <- fitness_func(N_k=kwds$N_k, T_k=kwds$T_k, dt=dt)

        } else if (data_type=='PointMeasures'){
            ak_raw <- rep(1, times=length(x)) / sigma^2
            bk_raw <- x / sigma^2

            kwds$a_k <- 0.5 * rev(cumsum(rev(ak_raw[1:K])))
            kwds$b_k <- -1 * rev(cumsum(rev(bk_raw[1:K])))

            fit_vec <- fitness_func(a_k=kwds$a_k, b_k=kwds$b_k)
        }


        A_R <- fit_vec - ncp_prior
        A_R[2:length(A_R)] <- A_R[2:length(A_R)] + best[1:(K - 1)]

        i_max <- which.max(A_R)
        last[K] <- i_max
        best[K] <- A_R[i_max]
    }

    #-----------------------------------------------------------------
    # Recover changepoints by iteratively peeling off the last block
    #-----------------------------------------------------------------

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

validate_input <- function(t, x = NULL, sigma = NULL, data_type='Events'){
    t <- as.numeric(t)

    # Get the order of 't' in ascending order
    order_t <- order(t)

    # Sort 't' in ascending order
    t <- t[order_t]

    # If 'x' and 'sigma' are not NULL, sort them according to 't'
    if (!is.null(x)) {
        x <- x[order_t]
    }

    if (!is.null(sigma)) {
        sigma <- sigma[order_t]
    }


    unq_t <- unique(t)
    unq_ind <- match(unq_t, t)
    unq_inv <- match(t, unq_t)

    if (is.null(x)) {
        if (!is.null(sigma)) {
            stop("If sigma is specified, x must be specified")
        } else {
            sigma <- 1
        }

        if (length(unq_t) == length(t)) {
            x <- rep(1, times = length(t))
        } else {
            x <- tabulate(unq_inv)
        }

        t <- unq_t
    } else {
        x <- as.numeric(x)

        x <- x + rep(0, times = length(t))

        t <- unq_t
        x <- x[unq_ind]
    }

    if (is.null(sigma)) {
        sigma <- 1
    } else {
        sigma <- as.numeric(sigma)
        if (!(length(sigma) %in% c(0, 1, length(t)))) {
            stop("sigma does not match the length of x")
        }
    }

    if (data_type=='Events'){
        if (!is.null(x) && any(x %% 1 > 0)) {
                stop("x must be integer counts for fitness='Events'")
        }
        return(list(t = t, x = x, sigma = sigma))

    } else if (data_type=='RegularEvents'){
        if (!is.null(x) && any(x %% 1 > 0)) {
                stop("x must be integer counts for fitness='RegularEvents'")
        }
        return(list(t = t, x = x, sigma = sigma))

    } else if (data_type=='PointMeasures'){
        if (is.null(x)){
            stop("x must be specified for fitness=`PointMeasures`")
        }
        return(list(t = t, x = x, sigma = sigma))
    }
}

compute_p0_prior <- function(N, p0=0.05){
            return(4 - log(73.53 * p0 * (N^-0.478)))
}

compute_ncp_prior <- function(N, p0=0.05, gamma=NULL) {
            if (!is.null(gamma)) {
                return(-log(gamma))
            } else if (is.null(gamma) & !is.null(p0)) {
                return(compute_p0_prior(N))
            } else {
                stop("ncp_prior cannot be computed as neither gamma nor p0 is defined")
            }
}

fitness_func <- function(N_k=NULL, T_k=NULL, dt=NULL, a_k=NULL, b_k=NULL, data_type='Events'){

    if (data_type=='Events'){
        return(N_k * (log(N_k / T_k)))

    } else if (data_type=='RegularEvents'){
        M_k <- T_k / dt
        N_over_M <- N_k / M_k

        eps <- 1e-8
        if (any(N_over_M > 1 + eps)) {
            warning("Regular Events: N/M > 1.  Is the time step correct?")
        }

        one_m_NM <- 1 - N_over_M
        N_over_M[N_over_M <= 0] <- 1
        one_m_NM[one_m_NM <= 0] <- 1

        return(N_k * log(N_over_M) + (M_k - N_k) * log(one_m_NM))

    } else if (data_type=='PointMeasures'){
        return((b_k * b_k) / (4 * a_k))
    }
}