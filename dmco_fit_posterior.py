'''
Run with the following arguments: data_num prior_num year param_type_list rate_type_list country sex
data_num : int, location of data
prior_num : int, 
year : int, 1990, 2005, or 2010
param_type_list : str, list of data_types separated by ' ', ex. 'p i f r' 
rate_type_list : str, list of rate_types separated by ' ', length must be equal to param_type_list or of length 1, ex. 'neg_binom binom poisson binom' 
country : str, ISO3
sex : str, 'male', 'female', or 'both'
'''

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

# assert that system arguments are correct
if len(sys.argv[5].split(' ')) != 1:
    assert len(sys.argv[5].split(' ')) == len(sys.argv[4].split(' ')), 'rate_type_list has the incorrect number of arguments--length must be 1 or match length of param_type_list'

# model id number for data
data_num = int(sys.argv[1])
prior_num = int(sys.argv[2])
year = int(sys.argv[3])
param_tye_list = sys.argv[4].split(' ')
rate_tye_list = sys.argv[5].split(' ')
country = str(sys.argv[6])
sex = str(sys.argv[7])

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

model = dmco.mvn(model, prior_num, param_type_list, country, sex, year, iter, burn, thin, rate_type=rate_type_list)

# generate estimates
dmco.save_posterior(dismod3.load_disease_model(data_num), model, country, sex, year, param_type_list)
