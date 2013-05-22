import dismod3
reload(dismod3)
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
year = int(sys.argv[3])
country = str(sys.argv[4])
sex = str(sys.argv[5])

# run settings
# iter=500
# burn=1
# thin=2
iter=20000
burn=10000
thin=10

# load country model and add country-specific mortality estimates
model = dmco.load_new_model(data_num, country, sex, cov='average')
# dmco.add_data(model, mortality, country, sex, year)
model.keep(start_year=2005)

model = dmco.mvn(model, prior_num, 'p', country, sex, year, iter, burn, thin, rate_type='neg_binom')

# generate estimates
dmco.save_posterior(dismod3.load_disease_model(data_num), model, country, sex, year, 
                    ['prevalence'])


