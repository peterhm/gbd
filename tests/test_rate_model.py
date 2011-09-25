""" Test Rate Model

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

import data
import rate_model
reload(rate_model)

def test_neg_binom_model_sim():
    # simulate negative binomial data
    pi_true = .01
    delta_true = 50

    n = pl.array(pl.exp(mc.rnormal(10, 1**-2, size=16)), dtype=int)
    k = pl.array(mc.rnegative_binomial(n*pi_true, delta_true), dtype=float)
    p = k/n

    # create and fit NB model
    pi = mc.Uniform('pi', 0., 1., value=.01)
    delta = mc.Uniform('delta', 0., 10000., value=1000.)
    vars = dict(pi=pi, delta=delta)
    vars.update(rate_model.neg_binom_model('sim', pi, delta, p, n))

    m = mc.MCMC(vars)
    m.sample(10000, 5000, 5)

    # compare estimate to ground truth
    assert pl.allclose(m.pi.stats()['mean'], pi_true, rtol=.1)
    lb, ub = m.delta.stats()['95% HPD interval']
    assert lb <= delta_true <= ub

def test_neg_binom_model_in_sample():
    # load real data
    d = data.ModelData.from_gbd_json('tests/hep_c_europe_western.json')

    # create and fit NB model
    pi = mc.Uniform('pi', 0., 1., value=.01)
    delta = mc.Uniform('delta', 0., 10000., value=1000.)
    vars = dict(pi=pi, delta=delta)
    vars.update(rate_model.neg_binom_model('sim', pi, delta,
                                           d.input_data['value'], d.input_data['effective_sample_size']))

    m = mc.MCMC(vars)
    m.sample(10000, 5000, 5)

    # generate graphical comparison of data and estimates
    k, n = m.p_pred.trace().shape
    d.input_data['standard_error'] = pl.sqrt(d.input_data['value']*(1-d.input_data['value'])/d.input_data['effective_sample_size'])

    sorted_indices = pl.argsort(d.input_data['value'].__array__())
    pl.errorbar(range(n),
                d.input_data['value'].__array__()[sorted_indices],
                yerr=d.input_data['standard_error'].__array__()*1.96,
                fmt='ks', mew=1, mec='w', ms=4,
                barsabove=True, zorder=10)
    pl.plot((pl.outer(pl.ones(k), range(n)) + pl.randn(k, n)*.1).flatten(),
            m.p_pred.trace()[:, sorted_indices].flatten(),
            'k,')

    pi_stats = m.pi.stats()
    pl.hlines(pi_stats['mean'], 0, n, 'w', linewidth=2, zorder=100)
    pl.hlines(pi_stats['mean'], 0, n, 'k', zorder=100)
    pl.hlines(pi_stats['95% HPD interval'], 0, n, 'w', linewidth=2, zorder=100)
    pl.hlines(pi_stats['95% HPD interval'], 0, n, 'k', linestyle='dashed', zorder=100)
    pl.axis([0,n,0,.3])

    pl.savefig('neg_binom_hep_c.png')

if __name__ == '__main__':
    import nose
    nose.runmodule()
    