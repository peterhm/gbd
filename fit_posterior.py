#!/usr/bin/python2.5
""" Generate a posterior estimate for a specific region, sex, and year

Examples
--------

$ python fit_posterior.py 3828 -r australasia -s male -y 2005 

>>> # ipython example
>>> from fit_posterior import *
>>> dm = dismod3.get_disease_model(3828)
>>> mort = dismod3.get_disease_model('all-cause_mortality')
>>> dm.data += mort.data
>>> import dismod3.gbd_disease_model as model
>>> model.fit(dm, method='map', keys=keys)
>>> model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=5, burn=5000, verbose=1)
>>> dismod3.post_disease_model(dm)
"""


# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 

import dismod3

def fit_posterior(id, region, sex, year):
    """ Fit posterior of specified region/sex/year for specified model

    Parameters
    ----------
    id : int
      The model id number for the job to fit
    region : str
      From dismod3.settings.gbd_regions, but clean()-ed
    sex : str, from dismod3.settings.gbd_sexes
    year : str, from dismod3.settings.gbd_years

    Example
    -------
    >>> import fit_posterior
    >>> fit_posterior.fit_posterior(2552, 'asia_east', 'male', '2005')
    """
    #print 'updating job status on server'
    #dismod3.log_job_status(id, 'posterior', '%s--%s--%s' % (region, sex, year), 'Running')

    dm = dismod3.load_disease_model(id)
    #dm.data = []  # for testing, remove all data
    keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])

    # fit the model
    import dismod3.gbd_disease_model as model
    model.fit(dm, method='map', keys=keys, verbose=1)     ## first generate decent initial conditions
    model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=5, burn=5000, verbose=1)     ## then sample the posterior via MCMC

    # generate plots of results
    dismod3.tile_plot_disease_model(dm, keys, defaults={})
    dm.savefig('dm-%d-posterior-%s-%s-%s-%s.png' % (id, 'all', region, sex, year))
    for param_type in dismod3.settings.output_data_types:
        keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex], type_list=[param_type])
        dismod3.tile_plot_disease_model(dm, keys, defaults={})
        dm.savefig('dm-%d-posterior-%s-%s-%s-%s.png' % (id, dismod3.utils.clean(param_type), region, sex, year))


    # save results (do this last, because it removes things from the disease model that plotting function, etc, might need
    dm.save('dm-%d-posterior-%s-%s-%s.json' % (id, region, sex, year), keys_to_save=keys)

    # update job status file
    #print 'updating job status on server'
    #dismod3.log_job_status(id, 'posterior',
    #                       '%s--%s--%s' % (region, sex, year), 'Completed')


def main():
    import optparse

    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    parser.add_option('-s', '--sex', default='male',
                      help='only estimate given sex (valid settings ``male``, ``female``, ``all``)')
    parser.add_option('-y', '--year', default='2005',
                      help='only estimate given year (valid settings ``1990``, ``2005``)')
    parser.add_option('-r', '--region', default='australasia',
                      help='only estimate given GBD Region')

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    fit_posterior(id, options.region, options.sex, options.year)
        

if __name__ == '__main__':
    main()