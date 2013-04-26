import dismod3
import covariate_model
import data as dm_data
import json
import pandas 
import pymc as mc
from time import clock
import pylab as pl

full_name = {'p': 'prevalence', 
             'i': 'incidence', 
             'r': 'remission', 
             'f': 'excess-mortality', 
             'pf': 'pf', 
             'csmr': 'r_specific', 
             'm_all': 'r_all',
             'm_with': 'r_with',
             'm': 'r_other',
             'smr': 'r_standard', 
             'rr': 'relative_risk', 
             'X': 'duration'}

def disease_info(obj):
    '''find disease number from disease name or
    find a disease name from a disease number'''
    bm_path = '/snfs1/Project/GBD/dalynator/yld/best_models.csv'
    bm_csv = pandas.read_csv(bm_path,index_col=None)
    dismod_models = bm_csv.groupby('dismod_model_number').apply(lambda df: df.ix[df.index[0], 'outcome_name'])
    dismod_models = dismod_models.drop([0], axis=0)
    dismod_models = dict(pl.sort(dismod_models))
    if type(obj)==str:
        # change Series object into dictionary
        from collections import defaultdict
        reversed_dict = defaultdict(list)
        for key,value in dismod_models.iteritems():
            reversed_dict[value].append(key)
        num = reversed_dict[obj]
        if num == []:
            print 'No DisMod-MR estimates for %s'%obj
        elif len(num) > 1:
            print 'DisMod-MR has more than one model for %s'%obj
            num = [int(k) for k in num]
        else:
            num = int(num[0])
        return num
    elif type(obj)==int:
        try:
            name = dismod_models[float(obj)]
        except:
            print 'No DisMod-MR best model for %s'%obj
            name = []
        return name
    else:
        print 'Invalid entry. Please enter disease number or name'
        return []

def load_new_model(disease, country='all', sex=['total', 'male', 'female'], cov='no'):
    '''create disease model with relavtive data
    cov : str
      method to handle covariates
      default is nothing ('no')
      options include, 
        - 'drop' : drop all covartiates
        - 'zero' : missing values replaced with 0
        - 'average' : missing values replaced with average of column
    '''
    model = dismod3.data.load('/home/j/Project/dismod/output/dm-%s'%disease)
    # keep relative data
    if (type(sex)==str) & (sex != 'total'): model.keep(areas=[country], sexes=[sex, 'total'])
    else: model.keep(areas=[country], sexes=sex)
    
    if (True in pl.isnan(pl.array(model.output_template.filter(like='x_')))) | (True in pl.isnan(pl.array(model.input_data.filter(like='x_')))): 
        print 'Covariates missing, %s method used'%(cov)
        col = model.input_data.filter(like='x_').columns
        for i in col:
            if cov == 'drop': 
                model.input_data = model.input_data.drop(i,1)
                model.output_template = model.output_template.drop(i,1)
            elif cov == 'zero': 
                model.input_data[i] = model.input_data[i].fillna([0])
                model.output_template[i] = model.output_template[i].fillna([0])
            elif cov == 'average': 
                model.input_data[i] = model.input_data[i].fillna([model.input_data[i].mean()])
                model.output_template[i] = model.output_template[i].fillna(model.output_template[i].mean())
    
    return model

def geo_info(country, disease):
    '''find country region from name'''
    global_model = dm_data.ModelData()
    hierarchy = json.load(open('/home/j/Project/dismod/dismod_status/prod/dm-%s/hierarchy.json'%(disease)))
    global_model.hierarchy.add_nodes_from(hierarchy['nodes'])
    global_model.hierarchy.add_edges_from(hierarchy['edges'])
    region = global_model.hierarchy.in_edges(country)[0][0]
    return region

def get_emp(disease, data_type, country, sex, year):
    '''load posterior estimates from GBD 2010 Study for country'''
    emp = pandas.read_csv('/home/j/Project/dismod/output/dm-%s/posterior/dm-%s-%s-%s-%s-%s.csv'%(disease, disease, full_name[data_type], geo_info(country, disease), sex, year), index_col=None)

    # keep only estimates from country
    mu_rate = emp[emp['Iso3'] == country].filter(like='Draw')
    
    return mu_rate
    
def find_fnrfx(model, disease, data_type, country, sex, year):
    '''add fixed and random effects from GBD as priors to new model'''
    # create dummy model to get appropriate Model.vars fields
    dummy = load_new_model(disease, country, sex)
    dummy.vars += dismod3.ism.age_specific_rate(dummy, data_type)
    vars = dummy.vars[data_type]
    
    # save random effects
    try:
        emp_re = pandas.read_csv('/home/j/Project/dismod/output/dm-%s/posterior/re-%s-%s+%s+%s.csv'%(disease, data_type, geo_info(country,disease), sex, year), index_col=0)
        for col in emp_re.index:
            model.parameters[data_type]['random_effects'][col] = dict(dist='Constant', 
                                                                      mu=emp_re.ix[col, 'mu_coeff'], 
                                                                      sigma=emp_re.ix[col, 'sigma_coeff'])
    except:
        pass
        
    # also save empirical prior on sigma_alpha, the dispersion of the random effects
    dm = dismod3.load_disease_model(disease)
    for n in vars['sigma_alpha']:
        try:
            dm_na = dm.get_empirical_prior(full_name[data_type])['new_alpha']
            model.parameters[data_type]['random_effects'][n.__name__] = dict(dist = dm_na[n.__name__]['dist'],
                                                                             mu = dm_na[n.__name__]['mu'], 
                                                                             sigma = dm_na[n.__name__]['sigma'], 
                                                                             lower = dm_na[n.__name__]['lower'], 
                                                                             upper = dm_na[n.__name__]['upper'])
        except:
            model.parameters[data_type]['random_effects'][n.__name__] = dict(dist = 'TruncatedNormal',
                                                                             mu = .05,
                                                                             sigma = .03**-2, 
                                                                             lower = 0.01, 
                                                                             upper = 0.5)
    # save fixed effects    
    emp_fe = pandas.read_csv('/home/j/Project/dismod/output/dm-%s/posterior/fe-%s-%s+%s+%s.csv'%(disease, data_type, geo_info(country,disease), sex, year), index_col=0)
    for n, col in zip(vars['beta'], vars['X'].columns):
        model.parameters[data_type]['fixed_effects'][col] = dict(dist = 'Constant', 
                                                                 mu = emp_fe.ix[col, 'mu_coeff'], 
                                                                 sigma = emp_fe.ix[col, 'sigma_coeff'])    

def mare(model, data_type):
    try:
        pred = model.vars[data_type]['p_pred'].trace().mean(0)
    except:
        pred = 0    
    obs = model.get_data(data_type)['value']
    mare = pl.median((abs(pred - obs)/obs)*100)
    return mare

def mvn(model, disease, data_param, country, sex, year, iter, burn, thin, var_inflation=1, log_space=False):
    '''multivariate normal (variance inflation optional)
    consistent : str, 'asr' (single rate type model) or 'ism' (compartmental modeling)
    '''
    if data_param == 'consistent':
        data_types = ['f', 'i', 'p', 'r'] #, 'X']
    else:
        data_types = data_param

    # set priors
    priors = {}
    for data_type in data_types:
        # get prior for each data_type
        gbd_est = get_emp(disease, data_type, country, sex, year)
        priors[data_type] = gbd_est
        
        find_fnrfx(model, disease, data_type, country, sex, year)
    
    # add vars
    if data_param == 'consistent':
        model.vars += dismod3.ism.consistent(model, country, sex, year)
    else:
        model.vars += dismod3.ism.age_specific_rate(model, data_type, country, sex, year, mu_age_parent=None, sigma_age_parent=None)

    # add gamma priors and mc.potential
    for data_type in data_types:
        gbd_est = priors[data_type]
        if log_space == True:
            mu_rate = pl.log(pl.array(gbd_est))
            mu_rate_mean = pl.log(pl.array(gbd_est)).mean(1)
            covar = pl.cov(pl.log(pl.array(gbd_est))-pl.array([mu_rate_mean for _ in range(1000)]).T)
        else:
            mu_rate = pl.array(gbd_est)
            mu_rate_mean = pl.array(gbd_est).mean(1)
            covar = pl.cov(pl.array(gbd_est)-pl.array([mu_rate_mean for _ in range(1000)]).T)

        try:
            for gamma_k, a_k in zip(model.vars[data_type]['gamma'], model.parameters[data_type]['parameter_age_mesh']):
                if log_space == True:
                    gamma_k.value = mu_rate_mean[a_k]
                else:
                    gamma_k.value = pl.log(mu_rate_mean[a_k]+1e-10)
        except KeyError:
            pass
                
        if log_space == True:
            @mc.potential
            def parent_similarity(mu_child=model.vars[data_type]['mu_age'], mu=mu_rate_mean, C=covar*var_inflation):
                return mc.mv_normal_cov_like(pl.log(mu_child), mu, C+1.e-6*pl.eye(101))
        else:
            @mc.potential
            def parent_similarity(mu_child=model.vars[data_type]['mu_age'], mu=mu_rate_mean, C=covar*var_inflation):
                #return mc.normal_like(mu_child, mu, (pl.diag(C)+1)**-1)
                return mc.mv_normal_cov_like(mu_child, mu, C+1.e-12*pl.eye(101))
        model.vars[data_type]['parent_similarity'] = parent_similarity
    
    start = clock()
    if data_param == 'consistent':
        dismod3.fit.fit_consistent(model, iter=iter, thin=thin, burn=burn)
    else:
        dismod3.fit.fit_asr(model, data_type, iter=iter, burn=burn, thin=thin)
        
    time = clock() - start
    
    return model, priors, time, mare(model, data_type)    

def open_mortality():
    # load file
    try:
        import scikits.statsmodels.iolib as pd        
    except:
        import scikits.statsmodels.lib.io as pd
    mortality = pandas.DataFrame(pd.genfromdta('/home/j/Project/Mortality/GBD Envelopes/04. Lifetables/02. MORTMatch/cluster/results/compiled/iso3_lt_mean_uncertainty.dta'))
    # keep desired variables
    mortality = mortality.ix[:,0:7]
    mortality.columns = ['area', 'sex', 'year_start', 'age_start', 'lower_ci', 'upper_ci', 'value']
    # add input data
    mortality = mortality[((mortality['year_start']==1990)|(mortality['year_start']==2005)|(mortality['year_start']==2010))]
    return mortality

def add_data(model, mortality, country, year):
    # select desired area and year
    data = mortality[mortality['area'] == country]
    data = data[data['year_start'] == year]
    
    # add input data
    # special case for ages
    data['age_end'] = data['age_start']+5.
    data['age_end'][data['age_start'] == 0] = 1
    data['age_end'][data['age_start'] == 1] = 5

    data['age_weights'] = pl.nan
    data['data_type'] = 'm_all'
    data['effective_sample_size'] = pl.nan
    data['standard_error'] = pl.nan
    data['year_end'] = data['year_start']

    model.input_data = model.input_data.append(data, ignore_index=True)

def save_posterior(dm, model, country, sex, year, rate_type_list):
    ''' Save country level posterior in a csv file, and put the file in the 
    directory job_working_directory/posterior/country_level_posterior_dm-'id'
    dm : dismod3.load_disease_model(model_num)
    model : dataModel
    country : str
      iso country code
    sex : str
      one of 'male', 'female', or 'total'
    year : int
      one of 1990, 2005, 2010
    rate_type_list : list
      i.e. ['incidence', 'prevalence', 'remission', 'excess-mortality', 'duration', 'prevalence_x_excess-mortality']
    '''
    # job working directory
    job_wd = dismod3.settings.JOB_WORKING_DIR % dm.id

    # directory to save the file
    dir = job_wd + '/posterior/'

    # create posteriors for rate types
    for rate_type in rate_type_list:
        try:
            # make an output file
            filename = 'dm-%s-%s-%s-%s-%s.csv' % (str(dm.id), rate_type, country, sex, year)
            print('writing csv file %s' % (dir + filename))

            # set prior bounds
            t = {'incidence': 'i', 'prevalence': 'p', 'remission': 'r', 'excess-mortality': 'f',
                 'prevalence_x_excess-mortality': 'pf', 'duration': 'X'}[rate_type]
            if t in model.vars:
                if t in model.parameters and 'level_bounds' in model.parameters[t]:
                    lower=model.parameters[t]['level_bounds']['lower']
                    upper=model.parameters[t]['level_bounds']['upper']
                else:
                    lower=0
                    upper=pl.inf

                posterior = covariate_model.predict_for(model,
                                                        model.parameters[t],
                                                        country, sex, year,
                                                        country, sex, year,
                                                        True,  # population weighted averages
                                                        model.vars[t],
                                                        lower, upper).T

                # create correct number of draws
                if posterior.shape[1] < 1000: 
                    draws = posterior.shape[1]
                    print 'Output file will have fewer than 1000 draws'
                else: 
                    draws = 1000
                    posterior = posterior[:,0:draws] 
                # create file
                file = pandas.DataFrame(posterior, columns=['Draw%d'%i for i in range(draws)])
                file['Iso3'] = country
                file['Population'] = dismod3.neg_binom_model.population_by_age[(country, str(year), sex)]
                file['Rate type'] = rate_type
                file['Age'] = model.parameters['ages']
                
                # save file
                file.to_csv(dir+filename, index=False)

        except IOError, e:
            print 'WARNING: could not save country level output for %s' % rate_type
            print e

def plot_fits(disease, prior, year):
    '''Plot country fits grouped by region'''
    dir = '/home/j/Project/dismod/dismod_status/prod/'

    world = load_new_model(disease)
    for region in [region for regions in [world.hierarchy.neighbors(sr[1]) for sr in world.hierarchy.edges('all')] for region in regions]:
        country_list = world.hierarchy.neighbors(region)
        pl.figure(figsize=(24,len(country_list)*8))
        for c,country in enumerate(country_list):
            for s,sex in enumerate(['male', 'female']):
                model = load_new_model(disease, country, sex)
                for j,data_type in enumerate(['p','i','r','X','f']):
                    pl.subplot(len(country_list)*2,5,(j+1)+(s*5)+c*10)
                    dismod3.graphics.plot_data_bars(model.get_data(data_type))
                    # get estimates
                    est = pandas.read_csv(dir+'dm-%s/posterior/dm-%s-%s-%s-%s-%s.csv' % (disease, disease, full_name[data_type], country, sex, year),index_col=None)
                    est = est.filter(like='Draw')
                    gbd_est = get_emp(prior, data_type, country, sex, year)
                    
                    ymax = 0.
                    if max(est.mean(1)) > ymax: ymax = max(est.mean(1))
                    if max(gbd_est.mean(1)) > ymax: ymax = max(gbd_est.mean(1))
                    
                    # plotting
                    pl.plot(pl.array(est.mean(1)), 'k-', label='DM-CO')
                    pl.plot(pl.array(gbd_est.mean(1)), 'r-', label='GBD2010')
                    pl.plot(mc.utils.hpd(pl.array(gbd_est).T, .05), 'r:')
                    pl.plot(mc.utils.hpd(pl.array(est).T, .05), 'k:')
                    pl.title(country +' '+ data_type +' '+ sex +' '+ str(year) )
                    pl.axis([-5, 105, -ymax*.05, ymax*1.1])
                    pl.legend(loc='upper left')
                    
        pl.savefig(dir+'/dm-%s/image/%s_%s.png'%(disease, region, year))


