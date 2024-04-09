# Define the list of symbols accessible when importing the module
exported_symbols <- c("FitnessFunc", "events", "reg_events", "real_data", "bayesian_blocks", "Events", "RegEvents", "RealData")

bayesian_blocks <- function(time, data = NULL, err = NULL, fitness = "events", ...) {
  # Define a dictionary-like structure for fitness functions
  FITNESS_DICT <- list(
    events = "Events",
    regevents = "RegularEvents",
    realdata = "RealData"
  )
  
  # Get the appropriate fitness function class name
  fitness <- ifelse(fitness %in% names(FITNESS_DICT), FITNESS_DICT[[fitness]], fitness)

  print(fitness)
  
  # Check if the fitness parameter is a valid class name
  if (!exists(fitness) || !is.function(get(fitness))) {
    stop("fitness parameter not understood")
  }
  
  # Create an instance of the fitness function class
  fitfunc <- do.call(fitness, list(...))
  
  # Call the fit method of the fitness function
  return(fitfunc$fit(time, data, err))


    # The FITNESS_DICT dictionary is created using the list()
    # function, with the keys as character strings and the values 
    # as the corresponding R objects (Events, RegEvents, RealData).
    # To retrieve the value associated with the fitness key from
    # the dictionary, we use double square brackets [[ and assign it to the
    # variable fitness.
    # Using is() and inherits() functions to check if fitness is a 
    # class and a subclass of FitnessFunc. If it is, we create an instance 
    # of the class using do.call() and assign it to fitfunc. 
    # If fitness is already an instance of FitnessFunc, we assign it 
    # directly to fitfunc.
    # Finally, we call the fit() method on fitfunc with the arguments 
    # t, x, and err.
}

# Define the FitnessFunc class
FitnessFunc <- list( 
    p0 = 0.05,
    gamma = NULL,
    ncp_prior = NULL,
    validate_input = function(time, data = NULL, err = NULL) {
        time <- as.numeric(time)
        unq_time <- unique(time)
        unq_index <- match(unq_time, time)
        unq_inverse <- match(time, unq_time)

        if (is.null(data)) {
            if (!is.null(err)) {
                stop("If err is specified, data must be specified")
            } else {
                err <- 1
            }

            if (length(unq_time) == length(time)) {
                data <- rep(1, length(time))
            } else {
                data <- tabulate(unq_inverse)
            }

            time <- unq_time
        } else {
            data <- as.numeric(data)
            data <- data + rep(0, length(time))
            time <- unq_time
            data <- data[unq_index]
        }

        if (is.null(err)) {
            err <- 1
        } else {
            err <- as.numeric(err)
        }

        return(list(time = time, data = data, err = err))
        },
    p0_prior = function(N) {
        return(4 - log(73.53 * self$p0 * (N^-0.478)))
    },
    fitness_args = function() {
        return(names(formals(self$fit)))
    },
    c_ncp_prior = function(N) {
        if (!is.null(self$gamma)) {
            return (-log(self$gamma))
        } else if (!is.null(self$p0)) {
            return (self$p0_prior(N))
        }
    },
    fit = function(time, data=NULL, err=NULL) {
        input <- self$validate_input(time, data, err)
        time <- input$time
        data <- input$data
        err <- input$err
        
        # compute values needed for computation
        if ("a_k" %in% super$fitness_args()) {
            ak_raw <- rep(1, length(data))/ err^2
        }
        if ("b_k" %in% super$fitness_args()) {
            bk_raw <- data / err^2
        }

        # create length-(N+1) array of cell edges
        edges <- c(time[1], 0.5 * (time[-length(time)] + time[-1]), time[length(time)])
        block_len <- time[length(time)] - edges

        # arrays to store the best configuration
        N <- length(time)
        best <- rep(0, N)
        last <- rep(0, N)

        # compute ncp_prior if not defined
        if (is.null(super$ncp_prior)) {
            ncp_prior <- super$c_ncp_prior(N)
        } else {
            ncp_prior <- super$ncp_prior
        }

        # ----------------------------------------------------------------
        # Start with first data cell; add one cell at each iteration
        # ----------------------------------------------------------------
        for (R in 0:(N-1)) {
            kwds <- list()

            if ("T_k" %in% super$fitness_args()) {
                kwds$T_k <- block_len[1:(R+1)] - block_len[R+1]
            }

            if ("N_k" %in% super$fitness_args()){
                kwds$N_k <- rev(cumsum(rev(data[1:(R+1)])))
            }

            if ("a_k" %in% super$fitness_args()) {
                kwds$a_k <- 0.5 * rev(cumsum(rev(ak_raw[1:(R+1)])))
            }

            if ("b_k" %in% super$fitness_args()){
                kwds$b_k <- -rev(cumsum(rev(bk_raw[1:(R+1)])))
            }

            if (inherits(super, "Events")) {
                fit_vec <- super$fitness(N_k = kwds$N_k, T_k = kwds$T_k)
            } else if (inherits(super, "RegEvents")) {
                fit_vec <- super$fitness(T_k = kwds$T_k, N_k = kwds$N_k)
            } else if (inherits(super, "RealData")) {
                fit_vec <- super$fitness(a_k = kwds$a_k, b_k = kwds$b_k)
            } else {
                stop("Unknown subclass")
            }

            A_R <- fit_vec - ncp_prior
            A_R[2:length(A_R)] <- A_R[2:length(A_R)] + best[1:R]

            i_max <- which.max(A_R)
            last[R] <- i_max
            best[R] <- A_R[i_max]
        }

        #-------------------------------------------------------------
        # find changepoints by iteratively peeling off the last block
        #-------------------------------------------------------------
        chg_pnts <- rep(0, N)
        i_cp <- N
        index <- N
        while (i_cp > 1) {
            i_cp <- i_cp-1
            chg_pnts[i_cp] <- index
            if (index == 1) {
                break
            }
            index <- last[index-1]
        }
        if (i_cp == 1) {
            chg_pnts[i_cp] <- 1
        }
        chg_pnts <- chg_pnts[i_cp:length(chg_pnts)]

        return(edges[chg_pnts])
    }
)

# Define the Events class
Events <- list(
    fitness = function(N_k, T_k) {
        return(N_k * (log(N_k / T_k)))
    },
    validate_input = function(time, data, err) {
        parent_result <- callSuper("validate_input", time = time, data = data, err = err)
        time <- parent_result$time
        data <- parent_result$data
        err <- parent_result$err
        if (!is.null(data) && any(data %% 1 > 0)) {
            stop("Events fitness function requires integer data")
        }
        return(list(time = time, data = data, err = err))
    }
)

# Define the RegEvents class
RegEvents <- list(
    initialize = function(dt, p0 = 0.05, gamma = NULL, ncp_prior = NULL) {
        self$dt <- dt
        callSuper(p0, gamma, ncp_prior)
    },
    validate_input = function(time, data, err) {
        parent_result <- callSuper("validate_input", time = time, data = data, err = err)
        time <- parent_result$time
        data <- parent_result$data
        err <- parent_result$err
        if (!all(data == 0 | data == 1)) {
            stop("Regular events must have only 0 and 1 in data")
        } 
        return(list(time = time, data = data, err = err))
    },
    fitness = function(T_k, N_k) {
        M_k <- T_k/self$dt
        N_over_M <- N_k/M_k

        eps <- 1e-8
        if (any(N_over_M > 1+eps)) {
            stop("regular events: N/M > 1.  Is the time step correct?")
        }
        one_m_NM <- 1 - N_over_M
        N_over_M[N_over_M <= 0] <- 1
        one_m_NM[one_m_NM <= 0] <- 1

        return(N_k * log(N_over_M) + (M_k - N_k) * log(one_m_NM))
    }
)      

# Define the RealData class
RealData <- list(
    initialize = function(p0 = 0.05, gamma = NULL, ncp_prior = NULL) {
        callSuper(p0 = p0, gamma = gamma, ncp_prior = ncp_prior)
    },
    fitness = function(a_k, b_k) {
        return((b_k * b_k) / (4 * a_k))
    },
    validate_input = function(time, data, err) { 
        if (is.null(data)) {
        stop("data must be specified for point measures")
    }
    return(callSuper(time = time, data = data, err = err))
    }
)

# Create an instance of Events class (optional)
events_instance <- Events()

# Create an instance of RegularEvents class (optional)
regevents_instance <- RegEvents(dt = 0.1)  # Provide the dt parameter if needed

# Create an instance of RealData class (optional)
realdata_instance <- RealData()

# Call bayesian_blocks with the desired fitness function instance
bb <- bayesian_blocks(time = x, fitness = "events")

library(ggplot2)
set.seed(42)
x <- c(stats::rcauchy(500, -5, 1.8),
       stats::rcauchy(2000, -4, 0.8),
       stats::rcauchy(500, -1, 0.3),
       stats::rcauchy(1000, 2, 0.8),
       stats::rcauchy(500, 4, 1.5))

# Truncate values to a reasonable range
x <- x[x > -15 & x < 15]

hist(x, breaks = 200, density = TRUE)

bb <- bayesian_blocks(time=x, fitness = "events")

