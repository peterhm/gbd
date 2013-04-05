import pylab as pl
import pandas 
import sys
import os

# download data to j drive
os.system('python download_model.py %s'%(sys.argv[1]))

# load country list
country_list = pandas.read_csv('/snfs1/DATA/IHME_COUNTRY_CODES/IHME_COUNTRYCODES.CSV', index_col=None)
country_list = country_list[country_list.ix[:,'ihme_indic_country'] == 1]
country_list = list(pl.unique(country_list['iso3']))

# launch on cluster
for country in ['USA', 'GBR']:
    for sex in ['male', 'female']:
        os.system('/usr/local/bin/SGE/bin/lx24-amd64/qsub -cwd /homes/peterhm/gbd/dmco_fit_posterior.sh %s %s %s %s %s' %(sys.argv[1], sys.argv[2], sys.argv[3], country, sex))

