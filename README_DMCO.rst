DisMod-CO
=========

DisMod-CO aims to create country-sex-year specific estimates that are not overshadowed by regional trends (unless national data is not available). Country-sex-year specific mortality (m_all) data were not used for the GBD2010 estimates since consistency was enforced at the regional level, but are a useful addition to the country specific patterns and trends.  Isolating each country-sex-year allows the data to speak for itself and regional priors only inform estimates in the absence of data.  It also allows the compartmental model to enforce consistency between incidence, prevalence, remission, and mortality at the country level.

There are 5 essential files necessary to run DMCO with DisMod-MR.  They should be in the same folder as all of the DisMod-MR code.  These files with launch country-sex-year specific on the cluster, create estimates, and compile the estimates and the graphics into a pdf. They are the following:

dmco_launch.py
dmco_fit_posterior.py
dmco_plot_fits_pdf.py
dmco_methods.py
test_dmco.py


dmco_launch.py
---------------------

This file begins the whole process by launching country-sex-year specific models on the cluster, creating estimates, and compiling the estimates and the graphics into a pdf.  

Run with the following arguments: data_num prior_num year param_type_list rate_type_list 
data_num : int, location of data
prior_num : int, 
year : int, 1990, 2005, or 2010
param_type_list : str, list of data_types separated by ' ', ex. 'p i f r' 
rate_type_list : str, list of rate_types separated by ' ', length must be equal to param_type_list or of length 1, ex. 'neg_binom binom poisson binom' 

ex.  If I wanted to create CKD estimates in 2010 using priors from the GDB2010 study. I would run the following for a prevalence-only fit using the negative binomial distribution.

run dmco_launch.py 42165 39264 2010 'p' 'neg_binom' 

Or I could run a consistent fit using priors for prevalence, incidence, remission, and excess-mortality with a negative binomial model.

run dmco_launch.py 42165 39264 2010 'p i r f' 'neg_binom' 

Or I could run a consistent fit using priors for prevalence, incidence, remission, and excess-mortality with different distributions for each data type.

run dmco_launch.py 42165 39264 2010 'p i r f' 'neg_binom poisson normal binom' 


dmco_fit_posterior.py
---------------------

This file creates estimates for a specific country-sex-year saved in 
'/home/j/Project/dismod/dismod_status/prod/dm-%s/posterior/%s'%(data_num,folder) 
and records the length of the run in '/home/j/Project/dismod/dismod_status/prod/dm-%s/posterior/stdout/%s.csv'%(data_num,country+str(year)+sex)).  

Within this file, you may need to change the run settings (iter=102, burn=1, thin=1).  You can also change the folder name where the posterior estimates are saved, this is an option of dmco.save_posterior().  This is also where you specify if country-sex-year mortality is included and how many years of data to include in the moddel (current settings, +/- 2 years).  The years of data can be adjusted with model.keep(start_year=year-2,end_year=year+2).  NOTE: If you change the number of years to include in the model, it will also need to be changed in dmco_methods.plot_fits_pdf()!!!

Run with the following arguments: data_num prior_num year param_type_list rate_type_list country sex
data_num : int, location of data
prior_num : int, 
year : int, 1990, 2005, or 2010
param_type_list : str, list of data_types separated by ' ', ex. 'p i f r' 
rate_type_list : str, list of rate_types separated by ' ', length must be equal to param_type_list or of length 1, ex. 'neg_binom binom poisson binom' 
country : str, ISO3
sex : str, 'male', 'female', or 'both'

ex.  If I wanted to create CKD estimates for men in the Unites States in 2010 using priors from the GDB2010 study. I would run the following for a prevalence-only fit using the negative binomial distribution.

run dmco_fit_posterior.py 42165 39264 2010 'p' 'neg_binom' 'USA' 'male'

Or I could run a consistent fit using priors for prevalence, incidence, remission, and excess-mortality with a negative binomial model.

run dmco_fit_posterior.py 42165 39264 2010 'p i r f' 'neg_binom' 'USA' 'male'

Or I could run a consistent fit using priors for prevalence, incidence, remission, and excess-mortality with different distributions for each data type.

run dmco_fit_posterior.py 42165 39264 2010 'p i r f' 'neg_binom poisson normal binom' 'USA' 'male'


dmco_plot_fits_pdf.py
---------------------

This file compiles all of the graphics of country-sex-year specific estimates and the priors of the data types of interest into a pdf. Saved 
'/home/j/Project/dismod/dismod_status/prod/dm-%s/image/%s_w_prior_%s_%s.pdf'%(disease, prior, year, filename)
The estimates of all jobs are also compiled and recorded at this location '/home/j/Project/dismod/dismod_status/prod/dm-%s/posterior/stdout/timesheet.csv'%(data_num)

Run with the following arguments: data_num prior_num year param_type_list
data_num : int, location of data
prior_num : int, location of the priors
year : int, 1990, 2005, or 2010
param_type_list : str, list of rate_types separated by ' ', ex. 'p i f r' 

You may need to change the data_types that are plotted.  For example, you have a prevalence-only model and don't want to see mortality data or you only wanted priors for prevalence and incidence but also want to see remission estimates.  This can be done by changing param_type_list.append().  You can also change the filename of the pdf which is an option of dmco_methods.plot_fits_pdf().

ex. After creating country-sex-year I wanted to create CKD estimates for men in the Unites States in 2010 using priors from the GDB2010 study for prevalence, incidence, remission, and excess-mortality.

run dmco_plot_fits_pdf.py 42165 39264 2010 'p i r f'


dmco_methods.py
---------------

This file contains all of the functions and statistics for DisMod-CO.  It is well documented.


test_dmco.py
------------

This is a test suite for DisMod-CO.

run test_dmco.py

