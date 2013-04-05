import pylab as pl
import pandas 
import sys
import os

# download data to j drive
os.system('python download_model.py %s'%(data_num))

# load country list
country_list = pandas.read_csv('/snfs1/DATA/IHME_COUNTRY_CODES/IHME_COUNTRYCODES.CSV', index_col=None)
country_list = country_list[country_list.ix[:,'ihme_indic_country'] == 1]
country_list = list(pl.unique(country_list['iso3']))

# launch on cluster
for country in country_list:
    for sex in ['male', 'female']:
        os.system('/usr/local/bin/SGE/bin/lx24-amd64/qsub -cwd -N ' + name + ' /home/j/Project/Models/dismodmr_rate_validation/model_comparison.sh %d %s %d' %(sys.argv[1], sys.argv[2], sys.argv[3]))

