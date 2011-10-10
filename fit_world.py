#!/usr/bin/python2.5
""" Generate empirical prior of specified parameter type

Expects the disase model json to be saved already.
"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 


import dismod3

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

import consistent_model
import covariate_model
import fit_model
import graphics

reload(consistent_model)
reload(fit_model)


def fit_world(dm, map_only=False):
    """ Fit consistent for all data in world

    Parameters
    ----------
    id : int
      The model id number for the job to fit

    Example
    -------
    >>> import fit_world
    >>> dm = fit_world.dismod3.load_disease_model(1234)
    >>> fit_world.fit_world(dm)
    """

    dir = dismod3.settings.JOB_WORKING_DIR % dm.id

    ## load the model from disk or from web
    import simplejson as json
    import data
    reload(data)

    model = data.ModelData.from_gbd_jsons(json.loads(dm.to_json()))

    ## next block fills in missing covariates with zero
    for col in model.input_data.columns:
        if col.startswith('x_'):
            model.input_data[col] = model.input_data[col].fillna(0.)
    # also fill all covariates missing in output template with zeros
    model.output_template = model.output_template.fillna(0)

    ### For testing:
    ## speed up computation by reducing number of knots
    ## for t in 'irf':
    ##    model.parameters[t]['parameter_age_mesh'] = [0, 20, 40, 60, 80, 100]

    vars = consistent_model.consistent_model(model,
                                             root_area='all', root_sex='total', root_year='all',
                                             priors={})

    ## fit model to data
    if map_only:
        dm.map, dm.mcmc = fit_model.fit_consistent_model(vars, 105, 0, 1, 100)
    else:
        dm.map, dm.mcmc = fit_model.fit_consistent_model(vars, 4040, 2000, 20, 100)

    dm.model = model
    dm.vars = vars

    for t in 'i r f p rr pf'.split():
        param_type = dict(i='incidence', r='remission', f='excess-mortality', p='prevalence', rr='relative-risk', pf='prevalence_x_excess-mortality')[t]
        #graphics.plot_one_type(model, vars[t], {}, t)
        for a in model.hierarchy['all'].keys() + [dismod3.utils.clean(a) for a in dismod3.settings.gbd_regions]:
            print 'generating empirical prior for %s' % a
            for s in dismod3.settings.gbd_sexes:
                for y in dismod3.settings.gbd_years:
                    key = dismod3.utils.gbd_key_for(param_type, a, y, s)
                    emp_priors = covariate_model.predict_for(model.output_template, model.hierarchy,
                                                             'all', 'total', 'all',
                                                             a, dismod3.utils.clean(s), int(y),
                                                             vars[t])
                    n = len(emp_priors)
                    emp_priors.sort(axis=0)
                    dm.set_mcmc('emp_prior_mean', key, emp_priors.mean(0))
                    dm.set_mcmc('emp_prior_std', key, emp_priors.std(0))
    
                    #pl.plot(model.parameters['ages'], dm.get_mcmc('emp_prior_mean', key), 'r-')

        from fit_emp_prior import store_effect_coefficients
        store_effect_coefficients(dm, vars[t], param_type)

    
        if 'p_pred' in vars[t]:
            graphics.plot_one_ppc(vars[t], t)
            pl.savefig(dir + '/prior-%s-ppc.png'%param_type)

            graphics.plot_one_effects(vars[t], t, model.hierarchy)
            pl.savefig(dir + '/prior-%s-effects.png'%param_type)

    graphics.plot_fit(dm.model, dm.vars, {}, {})
    pl.savefig(dir + '/prior.png')

    graphics.plot_convergence_diag(vars)
    pl.savefig(dir + '/prior-convergence.png')

    graphics.plot_trace(vars)
    pl.savefig(dir + '/prior-trace.png')
    
    # save results (do this last, because it removes things from the disease model that plotting function, etc, might need
    try:
        dm.save('dm-%d-prior-%s.json' % (dm.id, param_type))
    except IOError, e:
        print e

    return dm

def main():
    import optparse

    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    parser.add_option('-f', '--fast', default='False',
                      help='use MAP only')

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    dm = dismod3.load_disease_model(id)
    fit_world(dm, options.fast == 'True')
    return dm
      

if __name__ == '__main__':
    dm = main()
