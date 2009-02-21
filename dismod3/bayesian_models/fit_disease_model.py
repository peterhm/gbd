import numpy as np
import pymc as mc

import probabilistic_utils
from rate_single_binomial import setup_rate_model

def map_fit(dm):
    """
    tuck the pymc vars into the disease_model.vars, if they don't
    exist and then fit them with map
    """
    setup_disease_model(dm)

    map = mc.MAP(dm.vars)

    print "searching for maximum likelihood point estimate"
    iterlim = 5
    method = 'fmin'
    map.fit(verbose=10, iterlim=iterlim, method=method)

    probabilistic_utils.save_map(dm.i)
    probabilistic_utils.save_map(dm.r)
    probabilistic_utils.save_map(dm.f)
    probabilistic_utils.save_map(dm.p)

def mcmc_fit(dm):
    map_fit(dm)
    
    trace_len, thin, burn = 500, 10, 5000
    trace_len, thin, burn = 1, 1, 1
    mcmc = mc.MCMC(dm.vars)
    mcmc.sample(trace_len*thin+burn, burn, thin, verbose=1)

    probabilistic_utils.save_mcmc(dm.i)
    probabilistic_utils.save_mcmc(dm.r)
    probabilistic_utils.save_mcmc(dm.f)
    #import pdb; pdb.set_trace()
    probabilistic_utils.save_mcmc(dm.p)

def initialized_rate_vars(rf, rate_stoch=None):
    rf.fit['out_age_mesh'] = range(probabilistic_utils.MAX_AGE)

    # do normal approximation first, to generate a good starting point
    M,C = probabilistic_utils.normal_approx(rf)
    setup_rate_model(rf, rate_stoch)
    return rf.vars.items(), rf.rate_stoch
    

#################### Code for generating a hierarchical bayesian model of the asrf interactions
def setup_disease_model(dm):
    """
    
    Generate a list of the PyMC variables for the generic disease
    model, and store it as dm.vars

    See comments in the code for exact details on the model.
    """
    dm.vars = {}
    out_age_mesh = range(probabilistic_utils.MAX_AGE)

    dm.vars['i'], i = initialized_rate_vars(dm.i_in())
    dm.vars['r'], r = initialized_rate_vars(dm.r_in())
    dm.vars['f'], f = initialized_rate_vars(dm.f_in())
    # TODO: create m, from all-cause mortality
    m = np.zeros(len(out_age_mesh))
    
    # TODO: make error in C_0 a semi-informative stochastic variable
    C_0 = mc.Uniform('C_0', 0.0, 1.0)
    @mc.deterministic
    def S_0(C_0=C_0):
        return max(0.0, 1.0 - C_0)
    dm.vars['bins'] = [S_0,C_0]
    
    # iterative solution to difference equations to obtain bin sizes for all ages
    @mc.deterministic
    def S_C_D_M(S_0=S_0, C_0=C_0, i=i, r=r, f=f, m=m):
        S = np.zeros(len(out_age_mesh))
        C = np.zeros(len(out_age_mesh))
        D = np.zeros(len(out_age_mesh))
        M = np.zeros(len(out_age_mesh))
        
        S[0] = S_0
        C[0] = C_0
        D[0] = 0.0
        M[0] = 0.0
        
        for a in range(len(out_age_mesh)-1):
            S[a+1] = S[a]*(1-i[a]-m[a]) + C[a]*r[a]
            C[a+1] = S[a]*i[a]          + C[a]*(1-r[a]-m[a]-f[a])
            D[a+1] =                      C[a]*f[a]               + D[a]
            M[a+1] = S[a]*m[a]          + C[a]*m[a]                      + M[a]
                
        return S,C,D,M
    dm.vars['bins'] += [S_C_D_M]

    # prevalence = # with condition / (# with condition + # without)
    @mc.deterministic
    def p(S_C_D_M=S_C_D_M, tau_p=1./.01**2):
        S,C,D,M = S_C_D_M
        return C/(S+C)
    dm.vars['p'], p = initialized_rate_vars(dm.p_in(), rate_stoch=p)

