import numpy as np

def bayesian_blocks(t, x=None, sigma=None, fitness="events", p0=0.05, gamma=None, ncp_prior=None):
    FITNESS_DICT = {
        "events": fit_events_bayesian_blocks,
        "regular_events": fit_regevents_bayesian_blocks,
        "measures": fit_realdata_bayesian_blocks,
    }
    fitness = FITNESS_DICT.get(fitness, None)

    if fitness is None:
        raise ValueError("Invalid fitness function specified")

    return fitness(t, x, sigma, p0=p0, gamma=gamma, ncp_prior=ncp_prior)

def validate_input(t, x=None, sigma=None):
    """Validate inputs to the model.

    Parameters
    ----------
    t : array-like
        times of observations
    x : array-like, optional
        values observed at each time
    sigma : float or array-like, optional
        errors in values x

    Returns
    -------
    t, x, sigma : array-like, float or None
        validated and perhaps modified versions of inputs
    """
    # validate array input
    t = np.asarray(t, dtype=float)

    # find unique values of t
    t = np.array(t)
    
    unq_t, unq_ind, unq_inv = np.unique(t, return_index=True, return_inverse=True)

    # if x is not specified, x will be counts at each time
    if x is None:
        if sigma is not None:
            raise ValueError("If sigma is specified, x must be specified")
        else:
            sigma = 1

        if len(unq_t) == len(t):
            x = np.ones_like(t)
        else:
            x = np.bincount(unq_inv)

        t = unq_t

    # if x is specified, then we need to simultaneously sort t and x
    else:
        x = np.asarray(x, dtype=float)

        x += np.zeros_like(t)

        t = unq_t
        x = x[unq_ind]

    # verify the given sigma value
    if sigma is None:
        sigma = 1
    else:
        sigma = np.asarray(sigma, dtype=float)
        if sigma.shape not in [(), (1,), (t.size,)]:
            raise ValueError("sigma does not match the shape of x")

    return t, x, sigma

def p0_prior(N, p0=0.05):
    """
    Empirical prior, parametrized by the false alarm probability ``p0``
    """
    return 4 - np.log(73.53 * p0 * (N**-0.478))

def compute_ncp_prior(N, gamma=None, p0=0.05):
    """
    If ``ncp_prior`` is not explicitly defined, compute it from ``gamma``
    or ``p0``.
    """

    if gamma is not None:
        return -np.log(gamma)
    elif p0 is not None:
        return p0_prior(N, p0)
    else:
        raise ValueError(
            "``ncp_prior`` cannot be computed as neither "
            "``gamma`` nor ``p0`` is defined."
        )

def fit_bayesian_blocks(t, x=None, sigma=None, p0=0.05, gamma=None, ncp_prior=None):
    """Fit the Bayesian Blocks model given the specified fitness function.

    Parameters
    ----------
    t : array-like
        data times (one dimensional, length N)
    x : array-like, optional
        data values
    sigma : array-like or float, optional
        data errors
    p0 : float, optional
        false-alarm probability
    gamma : float, optional
        gamma parameter
    ncp_prior : float, optional
        prior for the ncp parameter

    Returns
    -------
    edges : ndarray
        array containing the (M+1) edges defining the M optimal bins
    """
    t, x, sigma = validate_input(t, x, sigma)

    # compute values needed for computation
    ak_raw = np.ones_like(x) / sigma**2
    bk_raw = x / sigma**2

    # create length-(N + 1) array of cell edges
    edges = np.concatenate([t[:1], 0.5 * (t[1:] + t[:-1]), t[-1:]])
    block_length = t[-1] - edges

    # arrays to store the best configuration
    N = len(t)
    best = np.zeros(N, dtype=float)
    last = np.zeros(N, dtype=int)

    # Compute ncp_prior if not defined
    if ncp_prior is None:
        ncp_prior = compute_ncp_prior(N, gamma, p0)

    # Start with first data cell; add one cell at each iteration
    for R in range(N):
        # Compute fit_vec: fitness of putative last block (end at R)
        kwds = {}

        # T_k: width/duration of each block
        kwds["T_k"] = block_length[: (R + 1)] - block_length[R + 1]

        # N_k: number of elements in each block
        kwds["N_k"] = np.cumsum(x[: (R + 1)][::-1])[::-1]

        # a_k: eq. 31
        kwds["a_k"] = 0.5 * np.cumsum(ak_raw[: (R + 1)][::-1])[::-1]

        # b_k: eq. 32
        kwds["b_k"] = -np.cumsum(bk_raw[: (R + 1)][::-1])[::-1]

        # evaluate fitness function
        fit_vec = fitness(**kwds, p0=p0, gamma=gamma, ncp_prior=ncp_prior)

        A_R = fit_vec - ncp_prior
        A_R[1:] += best[:R]

        i_max = np.argmax(A_R)
        last[R] = i_max
        best[R] = A_R[i_max]

    # Now find changepoints by iteratively peeling off the last block
    change_points = np.zeros(N, dtype=int)
    i_cp = N
    ind = N
    while i_cp > 0:
        i_cp -= 1
        change_points[i_cp] = ind
        if ind == 0:
            break
        ind = last[ind - 1]
    if i_cp == 0:
        change_points[i_cp] = 0
    change_points = change_points[i_cp:]

    return edges[change_points]


def events_fitness(N_k, T_k):
    # eq. 19 from Scargle 2013
    return N_k * (np.log(N_k / T_k))

def events_validate_input(t, x, sigma):
    t, x, sigma = validate_input(t, x, sigma)
    if x is not None and np.any(x % 1 > 0):
        raise ValueError("x must be integer counts for fitness='events'")
    return t, x, sigma

def fit_events_bayesian_blocks(t, x=None, sigma=None, p0=0.05, gamma=None, ncp_prior=None):
    t, x, sigma = events_validate_input(t, x, sigma)

    # compute values needed for computation
    ak_raw = np.ones_like(x) / sigma**2
    bk_raw = x / sigma**2

    # create length-(N + 1) array of cell edges
    edges = np.concatenate([t[:1], 0.5 * (t[1:] + t[:-1]), t[-1:]])
    block_length = t[-1] - edges

    # arrays to store the best configuration
    N = len(t)
    best = np.zeros(N, dtype=float)
    last = np.zeros(N, dtype=int)

    # Compute ncp_prior if not defined
    if ncp_prior is None:
        ncp_prior = compute_ncp_prior(N, gamma, p0)

    # Start with first data cell; add one cell at each iteration
    for R in range(N):
        # Compute fit_vec: fitness of putative last block (end at R)
        kwds = {}

        # T_k: width/duration of each block
        kwds["T_k"] = block_length[: (R + 1)] - block_length[R + 1]

        # N_k: number of elements in each block
        kwds["N_k"] = np.cumsum(x[: (R + 1)][::-1])[::-1]

        # evaluate fitness function
        fit_vec = events_fitness(kwds["N_k"], kwds["T_k"])

        A_R = fit_vec - ncp_prior
        A_R[1:] += best[:R]

        i_max = np.argmax(A_R)
        last[R] = i_max
        best[R] = A_R[i_max]

    # Now find changepoints by iteratively peeling off the last block
    change_points = np.zeros(N, dtype=int)
    i_cp = N
    ind = N
    while i_cp > 0:
        i_cp -= 1
        change_points[i_cp] = ind
        if ind == 0:
            break
        ind = last[ind - 1]
    if i_cp == 0:
        change_points[i_cp] = 0
    change_points = change_points[i_cp:]

    return edges[change_points]

import numpy as np
import warnings

def regular_events_fitness(dt, T_k, N_k):
    # Eq. C23 of Scargle 2013
    M_k = T_k / dt
    N_over_M = N_k / M_k

    eps = 1e-8
    
    one_m_NM = 1 - N_over_M
    N_over_M[N_over_M <= 0] = 1
    one_m_NM[one_m_NM <= 0] = 1

    return N_k * np.log(N_over_M) + (M_k - N_k) * np.log(one_m_NM)

def regular_events_validate_input(t, x, sigma):
    t, x, sigma = validate_input(t, x, sigma)
    if not np.all((x == 0) | (x == 1)):
        raise ValueError("Regular events must have only 0 and 1 in x")
    return t, x, sigma

def fit_regevents_bayesian_blocks(t, x=None, sigma=None, dt=1.0, p0=0.05, gamma=None, ncp_prior=None):
    t, x, sigma = regular_events_validate_input(t, x, sigma)

    # compute values needed for computation
    ak_raw = np.ones_like(x) / sigma**2
    bk_raw = x / sigma**2

    # create length-(N + 1) array of cell edges
    edges = np.concatenate([t[:1], 0.5 * (t[1:] + t[:-1]), t[-1:]])
    block_length = t[-1] - edges

    # arrays to store the best configuration
    N = len(t)
    best = np.zeros(N, dtype=float)
    last = np.zeros(N, dtype=int)

    # Compute ncp_prior if not defined
    if ncp_prior is None:
        ncp_prior = compute_ncp_prior(N, gamma, p0)

    # Start with first data cell; add one cell at each iteration
    for R in range(N):
        # Compute fit_vec: fitness of putative last block (end at R)
        kwds = {}

        # T_k: width/duration of each block
        kwds["T_k"] = block_length[: (R + 1)] - block_length[R + 1]

        # N_k: number of elements in each block
        kwds["N_k"] = np.cumsum(x[: (R + 1)][::-1])[::-1]

        # evaluate fitness function
        fit_vec = regular_events_fitness(dt, kwds["T_k"], kwds["N_k"])

        A_R = fit_vec - ncp_prior
        A_R[1:] += best[:R]

        i_max = np.argmax(A_R)
        last[R] = i_max
        best[R] = A_R[i_max]

    # Now find changepoints by iteratively peeling off the last block
    change_points = np.zeros(N, dtype=int)
    i_cp = N
    ind = N
    while i_cp > 0:
        i_cp -= 1
        change_points[i_cp] = ind
        if ind == 0:
            break
        ind = last[ind - 1]
    if i_cp == 0:
        change_points[i_cp] = 0
    change_points = change_points[i_cp:]

    return edges[change_points]

import numpy as np

def point_measures_fitness(a_k, b_k):
    # eq. 41 from Scargle 2013
    return (b_k * b_k) / (4 * a_k)

def point_measures_validate_input(t, x, sigma):
    if x is None:
        raise ValueError("x must be specified for point measures")
    return validate_input(t, x, sigma)

def fit_realdata_bayesian_blocks(t, x, sigma=None, p0=0.05, gamma=None, ncp_prior=None):
    t, x, sigma = point_measures_validate_input(t, x, sigma)

    # compute values needed for computation
    ak_raw = np.ones_like(x) / sigma**2
    bk_raw = x / sigma**2

    # create length-(N + 1) array of cell edges
    edges = np.concatenate([t[:1], 0.5 * (t[1:] + t[:-1]), t[-1:]])
    block_length = t[-1] - edges

    # arrays to store the best configuration
    N = len(t)
    best = np.zeros(N, dtype=float)
    last = np.zeros(N, dtype=int)

    # Compute ncp_prior if not defined
    if ncp_prior is None:
        ncp_prior = compute_ncp_prior(N, gamma, p0)

    # Start with first data cell; add one cell at each iteration
    for R in range(N):
        # Compute fit_vec: fitness of putative last block (end at R)
        kwds = {}

        # T_k: width/duration of each block
        kwds["T_k"] = block_length[: (R + 1)] - block_length[R + 1]

        # a_k: eq. 31
        kwds["a_k"] = 0.5 * np.cumsum(ak_raw[: (R + 1)][::-1])[::-1]

        # b_k: eq. 32
        kwds["b_k"] = -np.cumsum(bk_raw[: (R + 1)][::-1])[::-1]

        # evaluate fitness function
        fit_vec = point_measures_fitness(kwds["a_k"], kwds["b_k"])

        A_R = fit_vec - ncp_prior
        A_R[1:] += best[:R]

        i_max = np.argmax(A_R)
        last[R] = i_max
        best[R] = A_R[i_max]

    # Now find changepoints by iteratively peeling off the last block
    change_points = np.zeros(N, dtype=int)
    i_cp = N
    ind = N
    while i_cp > 0:
        i_cp -= 1
        change_points[i_cp] = ind
        if ind == 0:
            break
        ind = last[ind - 1]
    if i_cp == 0:
        change_points[i_cp] = 0
    change_points = change_points[i_cp:]

    return edges[change_points]
