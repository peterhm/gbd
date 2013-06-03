'''
Run with the following arguments: data_num prior_num year rate_type_list
data_num : int, location of data
prior_num : int, 
year : int, 1990, 2005, or 2010
rate_type_list : str, list of rate_types separated by ' ', ex. 'p i f r' 
'''

import pylab as pl
import pandas 
import sys
import os

# download data to j drive
os.system('/usr/local/epd-7.2-2/bin/python download_model.py %s'%(sys.argv[1]))

# load country list
country_list = pandas.read_csv('/snfs1/DATA/IHME_COUNTRY_CODES/IHME_COUNTRYCODES.CSV', index_col=None)
country_list = country_list[country_list.ix[:,'ihme_indic_country'] == 1]
country_list = list(pl.unique(country_list['iso3']))
country_list.remove('BMU')
country_list.remove('HKG')
country_list.remove('MAC')
country_list.remove('PRI')

# launch on cluster
name_list = []
for country in country_list: #['USA', 'GBR']:
    for sex in ['male', 'female']:
        name = country + str(year) + sex
        name_list.append(name)
        os.system('/usr/local/bin/SGE/bin/lx24-amd64/qsub -cwd -N ' + name + ' /homes/peterhm/gbd/dmco_fit_posterior.sh %s %s %s %s %s' %(sys.argv[1], sys.argv[2], sys.argv[3], country, sex))
        
# creating figures in .pdf
hold_str = '-hold_jid %s ' % ','.join(name_list)
os.system('/usr/local/bin/SGE/bin/lx24-amd64/qsub -cwd ' + hold_str + ' /homes/peterhm/gbd/' %(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]))
