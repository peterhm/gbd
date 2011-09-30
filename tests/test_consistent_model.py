""" Test data model

These tests are use randomized computation, so they might fail
occasionally due to stochastic variation
"""

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

import data
import data_model
import consistent_model
reload(consistent_model)
reload(data_model)
import data_simulation

def test_consistent_model_forward():
    m = data.ModelData()
    vars = consistent_model.consistent_model(m, 'all', 'total', 'all', {})

    vars['i']['gamma_bar'].value = pl.log(.01)
    vars['r']['gamma_bar'].value = pl.log(.0001)
    vars['f']['gamma_bar'].value = pl.log(.0001)
    print vars['p']['mu_age'].value[::10].round(3)

    vars['i']['gamma_bar'].value = pl.log(.02)
    vars['r']['gamma_bar'].value = pl.log(.0001)
    vars['f']['gamma_bar'].value = pl.log(.0001)
    print vars['p']['mu_age'].value[::10].round(3)

    vars['i']['gamma_bar'].value = pl.log(2.)
    vars['r']['gamma_bar'].value = pl.log(30.)
    vars['f']['gamma_bar'].value = pl.log(.0001)
    print vars['p']['mu_age'].value[::10].round(3)


def test_consistent_model_sim():
    m = data.ModelData()

    # generate simulated data
    n = 50
    sigma_true = .025
    a = pl.arange(0, 100, 1)
    pi_age_true = .0001 * (a * (100. - a) + 100.)

    m.input_data = data_simulation.simulated_age_intervals('p', n, a, pi_age_true, sigma_true)
    m.input_data['data_type'][-1] = 'r'  # make sure that there are multiple data types in the data set

    # create model and priors
    vars = consistent_model.consistent_model(m, 'all', 'total', 'all', {})

    # fit model
    m = mc.MCMC(vars)
    m.sample(1)

    return vars

if __name__ == '__main__':
    import nose
    nose.runmodule()
    