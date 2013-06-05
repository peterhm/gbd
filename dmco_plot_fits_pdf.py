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
import time 
import pandas

import dmco_methods as dmco
reload(dmco)

# time process
start = time.time()

data_num = int(sys.argv[1])
prior_num = int(sys.argv[2])
year = int(sys.argv[3])
param_type_list = sys.argv[4].split(' ')
name_list = list(pandas.read_csv('/home/j/Project/dismod/dismod_status/prod/dm-%s/posterior/stdout/name_list.csv'%(data_num), index_col=0)['0'])

# create graphics
param_type_list.append('m_all')
dmco.plot_fits_pdf(data_num, prior_num, year, param_type_list)

# compile times
timesheet = pandas.read_csv('/home/j/Project/dismod/dismod_status/prod/dm-%s/posterior/stdout/timesheet.csv'%(data_num), index_col=0)
for name in name_list:
    print name
    tmp = pandas.read_csv('/home/j/Project/dismod/dismod_status/prod/dm-%s/posterior/stdout/%s.csv'%(data_num,name), index_col=0)
    timesheet = pandas.concat([timesheet, tmp], ignore_index=True)

# record time 
timesheet = pandas.concat([timesheet, pandas.DataFrame({'job':'GRAPHICS', 'time':time.time() - start}, index=[0])], ignore_index=True)
timesheet.to_csv('/home/j/Project/dismod/dismod_status/prod/dm-%s/posterior/stdout/timesheet.csv'%(data_num))
    

