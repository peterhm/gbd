'''
Run with the following arguments: data_num prior_num year rate_type_list
data_num : int, location of data
prior_num : int, 
year : int, 1990, 2005, or 2010
param_type_list : str, list of data_types separated by ' ', ex. 'p i f r' 
rate_type_list : str, list of rate_types separated by ' ', length must be equal to param_type_list or of length 1, ex. 'neg_binom binom poisson binom' 
'''

import pylab as pl
import pandas 
import time 
import sys
import os

print sys.argv[4]
print sys.argv[5]

# time process
start = time.time()

# assert that system arguments are correct
if len(sys.argv[5].split(' ')) != 1:
    assert len(sys.argv[5].split(' ')) == len(sys.argv[4].split(' ')), 'rate_type_list has the incorrect number of arguments--length must be 1 or match length of param_type_list'

# download data to j drive
os.system('/usr/local/epd-7.3-2/bin/python download_model.py %s'%(sys.argv[1]))

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
        name = country + str(sys.argv[3]) + sex
        name_list.append(name)
        os.system('/usr/local/bin/SGE/bin/lx24-amd64/qsub -cwd -N ' + name + ' dmco_fit_posterior.sh "%s" "%s" "%s" "%s" "%s" "%s" "%s"' %(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], country, sex))

pandas.DataFrame(name_list).to_csv('/home/j/Project/dismod/dismod_status/prod/dm-%s/posterior/stdout/name_list.csv'%(sys.argv[1]))

# creating figures in .pdf
hold_str = '-hold_jid %s ' % ','.join(name_list)
os.system('/usr/local/bin/SGE/bin/lx24-amd64/qsub -cwd ' + hold_str + ' dmco_plot_fits_pdf.sh "%s" "%s" "%s" "%s"' %(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]))

# record time of process
pandas.DataFrame({'job':'LAUNCH','time':time.time() - start}, index=[0]).to_csv('/home/j/Project/dismod/dismod_status/prod/dm-%s/posterior/stdout/timesheet.csv'%(sys.argv[1]))
