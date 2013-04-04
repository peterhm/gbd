import dismod3
import fit_posterior 
import json

import pandas 
import pymc as mc

from time import clock
import pylab as pl
import sys
import os

import dmco_methods as dmco
reload(dmco)

# model id number for data
data_num = int(sys.argv[1])
prior_num = int(sys.argv[2])

# download data to j drive
os.system('python download_model.py %s'%(data_num))

# run settings
iter=500
burn=1
thin=2
# iter=20000
# burn=10000
# thin=10

# load mortality estimates
mortality = dmco.open_mortality()

country = 'USA'
sex='male'
# load country model and add country-specific mortality estimates
model = dmco.load_new_model(data_num, country, sex, cov='average')
dmco.add_data(model, mortality, country, year)
model, model_priors, model_t, model_mare = dmco.mvn(model, prior_num, 'consistent', country, sex, year, iter, burn, thin)

# generate estimates
dmco.save_country_posterior(dismod3.load_disease_model(data_num), model, country, sex, year, 
                            ['incidence', 'prevalence', 'remission', 'excess-mortality', 'duration', 'prevalence_x_excess-mortality'])


