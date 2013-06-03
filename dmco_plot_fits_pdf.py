'''
Run with the following arguments: data_num prior_num year rate_type_list
data_num : int, location of data
prior_num : int, location of the priors
year : int, 1990, 2005, or 2010
rate_type_list : str, list of rate_types separated by ' ', ex. 'p i f r' 
'''

import matplotlib
matplotlib.use('AGG')

import sys

import dmco_methods as dmco
reload(dmco)

dmco.plot_fits_pdf(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split(' '))