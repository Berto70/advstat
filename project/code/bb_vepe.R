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

validate.input = function(data, data_type='Events'){
    
    cells = data
    
    #using the data.frame structure
    if (!is.data.frame(cells)){
        
        cells = as.data.frame(cells)
    }
    
    if (data_type=='Events'){
        #count the number of times an event is registered
        cells = cells |> count(cells)
        #sort along times
        return(cells[order(cells[1]), ])
    }
    
    if (data_type=='RegularEvents'){
       # In each tick, there are either zero or one counts. 
        #just sorting along times
        return(cells[order(cells[1]), ])
    }
    
    if (data_type=='PointMeasures'){
        #just sorting along times
        
        return(cells[order(cells[1]), ])
    }
    
}

compute.p0_prior = function(N, p0=0.05){
    return(4 - log(73.53 * p0 * (N**-0.478))) 
}

compute.ncp_prior = function(gamma=NULL, N){
    if (is.null(gamma)){
        return(compute.p0_prior(N))
    }
    else{
        return(-log(gamma))
    }
}

compute.fitness = function(count_vec, width, dt=NULL, a_k=NULL, b_k=NULL, data_type='Events'){
    
    if (data_type=='Events'){
        return(count_vec*log(count_vec/width))
    }
    else if (data_type=='RegularEvents'){
        M = width/dt
        count_vec_over_M = count_vec/M
        
        esp = 1e-8
        if (sum(count_vec_over_M > 1+esp) > 0){
            cat("regular events: N/M > 1.  Is the time step correct?")
            break
        }
        
        one_m_NM = 1 - count_vec_over_M
        count_vec_over_M[count_vec_over_M <= 0] = 1
        one_m_NM[one_m_NM <= 0] = 1
        
        return(count_vec * log(count_vec_over_M) + (M  - count_vec)*log(one_m_NM))
    }
    else if(data_type=='PointMeasures'){
        return( (b_k*b_k)/(4*a_k) )
    }
    
}

bayesian_blocks = function(t, x=NULL, sigma=NULL, dt=NULL, gamma=NULL, ncp_prior=NULL, data_type='Events'){
    
    cells = validate.input(data, data_type)
    
    t = cells[[1]]
    x = cells[[2]]
    if (length(data)==3){sigma = cells[[3]]} #length(data) returns the number of columns of the dataframe 
    else{sigma = 1}
    N = length(t)
    

   
    #quantities needed for the computation of a_k, b_k
    ak_raw = rep(x = 1, times=length(x))/sigma**2
    bk_raw = x/sigma**2
    ck_raw = (x*x)/sigma**2
    
    edges = c(t[1],  0.5 * (t[2:length(t)] + t[1:(length(t)-1)] ), t[length(t)])
    block_length = t[length(t)] - edges
    
    best = rep(0, times = N)
    last = rep(1, times = N)
    
    #-----------------------------------------------------------------
    # Start with first data cell; add one cell at each iteration
    #-----------------------------------------------------------------
    
    for (K in seq(1:N)){
        # Compute the width and count of the final bin for all possible
        # locations of the K^th changepoint
        width = block_length[1:K] - block_length[K + 1] #T_K
        count_vec = rev(cumsum(rev(x[1:K]))) #N_K
        
        a_k = NULL
        b_k = NULL
        
        if (is.null(ncp_prior)){
            ncp_prior = compute.ncp_prior(gamma, N)
        }
        
        if (data_type !='Events' ){
            a_k = 0.5*rev(cumsum(rev(ak_raw[1:K])))
            b_k = -1 * rev(cumsum(rev(bk_raw[1:K])))
        }
        
        
        fit_vec = compute.fitness(count_vec = count_vec, width = width, a_k=a_k, b_k=b_k, dt=dt, data_type)
        
        post_vec = fit_vec - ncp_prior
        
        post_vec[2:length(post_vec)] = post_vec[2:length(post_vec)] + best[1:(K-1)]
        
        i_max = which.max(post_vec)
        last[K] = i_max
        best[K] = post_vec[i_max]
    }
    
    #-----------------------------------------------------------------
    # Recover changepoints by iteratively peeling off the last block
    #-----------------------------------------------------------------

    change_points =  rep(1, times=N)
    i_cp = N
    ind = N+1
    
    while (i_cp>1){
        change_points[i_cp] = ind
        if (ind == 1){ break }
        i_cp = i_cp-1
        ind = last[ind - 1]
    }
    
    if (i_cp==1){change_points[i_cp] = 1}
    change_points = change_points[i_cp:length(change_points)]

    return(edges[change_points])
}