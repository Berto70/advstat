options(warn=-1)

# Install and load the R6 package
if (!require(R6)) {
  install.packages("R6")
}
library(R6) # object-oriented programming in R

bayesian_blocks <- function(data, gamma=NULL, ncp_prior=NULL, dt=NULL, fitness='Events')

FitnessFunc <- R6Class(
    "FitnessFunc",
    
    public = list(
        p0 = NULL,
        gamma = NULL,
        ncp_prior = NULL,
        fitness = NULL,
        
        init = function(p0 = 0.05, gamma = NULL, ncp_prior = NULL, fitness = 'Events') {
            self$p0 <- p0
            self$gamma <- gamma
            self$ncp_prior <- ncp_prior
            self$fitness <- fitness
        },

        validate_input = function(data) {
            # Validate inputs to the model.
            
            cells <- data

            #using the data.frame structure
            if (!is.data.frame(cells)){            
                cells <- as.data.frame(cells)
            }
            
            if (self$fitness=='Events'){
                #count the number of times an event is registered
                cells <- cells |> count(cells)
                #sort along times
                return(cells[order(cells[1]),])
            }
            
            if (self$fitness=='RegularEvents'){
            # In each tick, there are either zero or one counts. 
                #just sorting along times
                return(cells[order(cells[1]), ])
            }
            
            if (self$fitness=='PointMeasures'){
                #just sorting along times            
                return(cells[order(cells[1]), ])
            }
        },

        p0_prior = function(N){
            
            # Empirical prior, parametrized by the false alarm probability ``p0``
            
            return (4 - log(73.53 * self$p0 * (N^-0.478)))
        },

        compute_ncp_prior = function(N) {
            if (!is.null(self$gamma)) {
            return(-log(self$gamma))
            } 
            else if (!is.null(self$p0)) {
            return(self$p0_prior(N))
            } 
            else {
            stop("``ncp_prior`` cannot be computed as neither ``gamma`` nor ``p0`` is defined.")
            }
        }

    )
)

# Define the Events class
Events <- R6Class(
  "Events",
  inherit = FitnessFunc,
  
  public = list(
    fitness = function(N_k, T_k) {
      return(N_k * (log(N_k / T_k)))
    },
    
    validate_input = function(t, x, sigma) {
      list(t = t, x = x, sigma = sigma) <- super$validate_input(t, x, sigma)
      if (!is.null(x) && any(x %% 1 > 0)) {
        stop("x must be integer counts for fitness='events'")
      }
      return(list(t = t, x = x, sigma = sigma))
    }
  )
)

# Define the RegularEvents class
RegularEvents <- R6Class(
    "RegularEvents",
    inherit = FitnessFunc,
    
    public = list(
        dt = NULL,
        p0 = NULL,
        gamma = NULL,
        ncp_prior = NULL,

        init = function(dt, p0=0.05, gamma=NULL, ncp_prior=NULL){
            self$dt <- dt
            self$p0 <- p0
            self$gamma <- gamma
            self$ncp_prior <- ncp_prior
        }
    )
)