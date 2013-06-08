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

import time
import pylab as pl
import sys
import os

import dmco_methods as dmco
reload(dmco)

# time process
start = time.time()

# assert that system arguments are correct
if len(sys.argv[5].split(' ')) != 1:
    assert len(sys.argv[5].split(' ')) == len(sys.argv[4].split(' ')), 'rate_type_list has the incorrect number of arguments--length must be 1 or match length of param_type_list'

# model id number for data
data_num = int(sys.argv[1])
prior_num = int(sys.argv[2])
year = int(sys.argv[3])
param_type_list = sys.argv[4].split(' ')
rate_type_list = sys.argv[5].split(' ')
country = str(sys.argv[6])
sex = str(sys.argv[7])

print 'data_num', data_num
print 'prior_num', prior_num 
print 'year', year 
print 'param_type_list', param_type_list
print 'rate_type_list', rate_type_list
print 'country', country
print 'sex', sex 

# run settings
iter=102
burn=1
thin=1
# iter=20000
# burn=10000
# thin=10

# load country model and add country-specific mortality estimates
model = dmco.load_new_model(data_num, country, sex, cov='average')
mortality = pandas.read_csv('/home/j/Project/dismod/dmco_mortality.csv')
dmco.add_data(model, mortality, country, sex, year)
model.keep(start_year=year-2)
model.keep(end_year=year+2)

model = dmco.mvn(model, prior_num, param_type_list, country, sex, year, iter, burn, thin, rate_type_list)

# generate estimates
dmco.save_posterior(data_num, model, country, sex, year, param_type_list)

# record time of process
pandas.DataFrame({'job':country+str(year)+sex,'time':time.time() - start}, index=[0]).to_csv('/home/j/Project/dismod/dismod_status/prod/dm-%s/posterior/stdout/%s.csv'%(data_num,country+str(year)+sex))
