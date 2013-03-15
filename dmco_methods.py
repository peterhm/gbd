import dismod3
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
    dismod_models = dict(sort(dismod_models))
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

def load_new_model(disease, country='all', sex=['total', 'male', 'female']):
    '''create disease model with relavtive data'''
    model = dismod3.data.load('/home/j/Project/dismod/output/dm-%s'%disease)
    # keep relative data
    if len(sex) == 1: model.keep(areas=[country], sexes=[sex])
    else: model.keep(areas=[country], sexes=sex)
    return model

def geo_info(country, disease):
    '''find country name and region'''
    global_model = load_new_model(disease, 'all', 'total')
    region = global_model.hierarchy.in_edges(country)[0][0]
    return region

def get_emp(disease, data_type, country, sex, year):
    '''load posterior estimates from GBD 2010 Study for country'''
    emp = pandas.read_csv('/home/j/Project/dismod/output/dm-%s/posterior/dm-%s-%s-%s-%s-%s.csv'%(disease, disease, full_name[data_type], geo_info(country, disease), sex, year), index_col=None)

    # keep only estimates from country
    mu_rate = emp[emp.Iso3 == country].filter(like='Draw')
    
    return mu_rate
    
def find_fnrfx(model, disease, data_type, country, sex, year):
    '''add fixed and random effects from GBD as priors to new model'''
    # create dummy model to get appropriate Model.vars fields
    dummy = load_new_model(disease, country, sex)
    dummy.vars += dismod3.ism.age_specific_rate(dummy, data_type)
    vars = dummy.vars[data_type]
    
    # save random effects
    emp_re = pandas.read_csv('/home/j/Project/dismod/output/dm-%s/posterior/re-%s-%s+%s+%s.csv'%(disease, data_type, geo_info(country,disease), sex, year), index_col=0)
    for col in emp_re.index:
        model.parameters[data_type]['random_effects'][col] = dict(dist='Constant', 
                                                                  mu=emp_re.ix[col, 'mu_coeff'], 
                                                                  sigma=emp_re.ix[col, 'sigma_coeff'])
        
    # also save empirical prior on sigma_alpha, the dispersion of the random effects
    dm = dismod3.load_disease_model(disease)
    dm_na = dm.get_empirical_prior(full_name[data_type])['new_alpha']
    for n in vars['sigma_alpha']:
        try:
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
    pred = model.vars[data_type]['p_pred'].trace().mean(0)
    obs = model.get_data(data_type).value
    mare = pl.median((abs(pred - obs)/obs)*100)
    return mare

def data_only(model, disease, data_type, country, sex, year, iter, burn, thin):
    '''data only'''
    find_fnrfx(model, disease, data_type, country, sex, year)
    model.vars += dismod3.ism.age_specific_rate(model, data_type, country, sex, year, mu_age_parent=None, sigma_age_parent=None)
    start = clock()
    dismod3.fit.fit_asr(model, data_type, iter=iter, burn=burn, thin=thin)
    time = clock() - start
    
    pred = dismod3.covariates.predict_for(model, model.parameters[data_type], country, sex, year, country, sex, year, True, model.vars[data_type], 0, 1)
    
    return model, pred, time, mare(model, data_type)    

def gbd_prior(model, disease, data_type, country, sex, year, iter, burn, thin, var_inflation=1):
    '''empirical prior from GBD2010 (variance inflation optional)'''
    # get prior
    gbd_est = get_emp(disease, data_type, country, sex, year)
    mu_rate = pl.array(gbd_est.mean(1))
    sigma_rate = pl.array(gbd_est.std(1))
    
    find_fnrfx(model, disease, data_type, country, sex, year)
    model.vars += dismod3.ism.age_specific_rate(model, data_type, country, sex, year, mu_age_parent=mu_rate, sigma_age_parent=sigma_rate*var_inflation)
    start = clock()
    dismod3.fit.fit_asr(model, data_type, iter=iter, burn=burn, thin=thin)
    time = clock() - start
    
    pred = dismod3.covariates.predict_for(model, model.parameters[data_type], country, sex, year, country, sex, year, True, model.vars[data_type], 0, 1)
    
    return model, pred, gbd_est, time, mare(model, data_type)    

def mvn(model, disease, data_type, country, sex, year, iter, burn, thin, var_inflation=1, log_space=True):
    '''multivariate normal (variance inflation optional)'''
    # get prior
    gbd_est = get_emp(disease, data_type, country, sex, year)
    if log_space == True:
        mu_rate = pl.log(pl.array(gbd_est))
        mu_rate_mean = pl.log(pl.array(gbd_est)).mean(1)
        covar = pl.cov(pl.log(pl.array(gbd_est)))
    else:
        mu_rate = pl.array(gbd_est)
        mu_rate_mean = pl.array(gbd_est).mean(1)
        covar = pl.cov(pl.array(gbd_est))
        
    find_fnrfx(model, disease, data_type, country, sex, year)
    model.vars += dismod3.ism.age_specific_rate(model, data_type, country, sex, year, mu_age_parent=None, sigma_age_parent=None)
    for gamma_k, a_k in zip(model.vars[data_type]['gamma'], model.parameters[data_type]['parameter_age_mesh']):
        gamma_k.value = mu_rate_mean[a_k]
        
    if log_space == True:
        @mc.potential
        def parent_similarity(mu_child=model.vars[data_type]['mu_age'], mu=mu_rate_mean, C=covar*var_inflation):
            return mc.mv_normal_cov_like(pl.log(mu_child), mu, C+.01*pl.eye(101))
    else:
        @mc.potential
        def parent_similarity(mu_child=model.vars[data_type]['mu_age'], mu=mu_rate_mean, C=covar*var_inflation):
            return mc.mv_normal_cov_like(mu_child, mu, C+.01*pl.eye(101))

    model.vars[data_type]['parent_similarity'] = parent_similarity
    start = clock()
    dismod3.fit.fit_asr(model, data_type, iter=iter, burn=burn, thin=thin)
    time = clock() - start
    
    pred = dismod3.covariates.predict_for(model, model.parameters[data_type], country, sex, year, country, sex, year, True, model.vars[data_type], 0, 1)
    
    return model, pred, gbd_est, time, mare(model, data_type)    

def discrete(model, disease, data_type, country, sex, year, iter, burn, thin, var_inflation=1):
    '''discrete'''
    # get prior
    gbd_est = get_emp(disease, data_type, country, sex, year)
    mu_rate_log = pl.log(pl.array(gbd_est))
    mu_rate_log_mean = pl.log(pl.array(gbd_est)).mean(1)
    cov_log = pl.cov(pl.log(pl.array(gbd_est)))

    find_fnrfx(model, disease, data_type, country, sex, year)
    model.vars += dismod3.ism.age_specific_rate(model, data_type, country, sex, year, mu_age_parent=None, sigma_age_parent=None)
    for gamma_k, a_k in zip(model.vars[data_type]['gamma'], model.parameters[data_type]['parameter_age_mesh']):
        gamma_k.value = mu_rate_log_mean[a_k]
    
    @mc.potential
    def parent_similarity(mu_child=model.vars[data_type]['mu_age'], mu=mu_rate_log, C=cov_log*var_inflation):
        i = int(rand()*1000)
        return mc.mv_normal_cov_like(pl.log(mu_child)[::10], mu[::10,i], C[::10,::10]+.001*pl.eye(11))

    model.vars[data_type]['parent_similarity'] = parent_similarity
    start = clock()
    dismod3.fit.fit_asr(model, data_type, iter=iter, burn=burn, thin=thin)
    time = clock() - start
    
    pred = dismod3.covariates.predict_for(model, model.parameters[data_type], country, sex, year, country, sex, year, True, model.vars[data_type], 0, 1)
    
    return model, pred, gbd_est, time, mare(model, data_type)

def mvn_inflation(model, disease, data_type, country, sex, year, iter, burn, thin, log_space=True):
    '''heterogeneity inflation for multivariate normal distribution'''
    # load regional model
    if sex != 'total': sexes=[sex, 'total']
    else: sexes = sex
    dm = load_new_model(disease, geo_info(country, disease), sexes)
    
    # create heterogeneity covariate and create prior
    dm.input_data['z_age'] = .5 * (dm.input_data.age_start + dm.input_data.age_end)
    dm.vars += dismod3.ism.age_specific_rate(dm, data_type, geo_info(country,disease), sex, year, mu_age_parent=None, sigma_age_parent=None)
    start = clock()
    dismod3.fit.fit_asr(dm, data_type, iter=iter, thin=thin, burn=burn)
    prior = dismod3.covariates.predict_for(dm, dm.parameters[data_type], geo_info(country, disease), sex, year, country, sex, year, True, dm.vars[data_type], 0, 1).T
    
    # inflate variance 
    mu_rate_mean = prior.mean(1)
    sigma_rate = pl.cov(prior)   
    zeta = dm.vars[data_type]['zeta'].stats()['mean']
    for i in range(101):
        for j in range(101):
            sigma_rate[i,j] *= pl.exp(i*zeta + j*zeta)   # FIXME: is this the correct way to gross-up the uncertainty???  Hannah will do a little simulation study to test
    
    find_fnrfx(model, disease, data_type, country, sex, year)
    model.vars += dismod3.ism.age_specific_rate(model, data_type, country, sex, year, mu_age_parent=None, sigma_age_parent=None)
    for gamma_k, a_k in zip(model.vars[data_type]['gamma'], model.parameters[data_type]['parameter_age_mesh']):
        gamma_k.value = mu_rate_mean[a_k]
    
    if log_space == True:
        @mc.potential
        def parent_similarity(mu_child=model.vars[data_type]['mu_age'], mu=mu_rate_mean, C=sigma_rate):
            return mc.mv_normal_cov_like(pl.log(mu_child), mu, C+.001*pl.eye(101))
    else:
        @mc.potential
        def parent_similarity(mu_child=model.vars[data_type]['mu_age'], mu=mu_rate_mean, C=sigma_rate):
            return mc.mv_normal_cov_like(mu_child, mu, C+.001*pl.eye(101))

    model.vars[data_type]['parent_similarity'] = parent_similarity
    dismod3.fit.fit_asr(model, data_type, iter=iter, burn=burn, thin=thin)
    time = clock() - start
    
    pred = dismod3.covariates.predict_for(model, model.parameters[data_type], country, sex, year, country, sex, year, True, model.vars[data_type], 0, 1)
    
    return model, pred, prior, time, mare(model, data_type)
    
def compare(name, disease, data_type, country, sex, year, ymax, iter, burn, thin):
    # METHODS
    # data only 
    do_model = load_new_model(disease, country, sex=sex)
    do_model, do_pred, do_t, do_mare = data_only(do_model, disease, data_type, country, sex, year, iter, burn, thin)
    # GBD prior 
    p_model = load_new_model(disease, country, sex=sex)
    p_model, p_pred, p_est, p_t, p_mare = gbd_prior(p_model, disease, data_type, country, sex, year, iter, burn, thin)
    # MVN
    mvn_model = load_new_model(disease, country, sex=sex)
    mvn_model, mvn_pred, mvn_est, mvn_t, mvn_mare = mvn(mvn_model, disease, data_type, country, sex, year, iter, burn, thin, var_inflation=1, log_space=False)
    # MVN log space
    mvnlog_model = load_new_model(disease, country, sex=sex)
    mvnlog_model, mvnlog_pred, mvnlog_est, mvnlog_t, mvnlog_mare = mvn(mvnlog_model, disease, data_type, country, sex, year, 
                                                                       iter, burn, thin, var_inflation=1, log_space=True)
    # Heterogeneity
    mvnhi_model = load_new_model(disease, country, sex=sex)
    mvnhi_model, mvnhi_pred, mvnhi_est, mvnhi_t, mvnhi_mare = mvn_inflation(mvnhi_model, disease, data_type, country, sex, year, iter, burn, thin, log_space=False)
    
    # PLOTTING
    plotting = [{'model':do_model, 'pred':do_pred, 'prior':zeros((101,2)), 'time':do_t, 'mare':do_mare, 'name':'Data only'},
                {'model':p_model, 'pred':p_pred, 'prior':p_est, 'time':p_t, 'mare':p_mare, 'name':'GBD prior'},
                {'model':mvn_model, 'pred':mvn_pred, 'prior':mvn_est, 'time':mvn_t, 'mare':mvn_mare, 'name':'MVN'},
                {'model':mvnlog_model, 'pred':mvnlog_pred, 'prior':mvnlog_est, 'time':mvnlog_t, 'mare':mvnlog_mare, 'name':'MVN log space'},
                {'model':mvni_model, 'pred':mvni_pred, 'prior':mvni_est, 'time':mvni_t, 'mare':mvni_mare, 'name':'MVN heterogeneous inflation'}]

    figure(figsize=(24,10))
    for p in range(len(plotting)):
        subplot(3,6,(p/3)*3+p+4)
        model = plotting[p]['model']
        dismod3.graphics.plot_data_bars(model.get_data(data_type))

        if (p in [2, 6]): errorbar(arange(101), plotting[p]['prior'].mean(1), 1.96*pl.array(plotting[p]['prior'].std(1)*10), color='k', capsize=0, elinewidth=.5, label='Prior error', alpha=.5)
        else: errorbar(arange(101), plotting[p]['prior'].mean(1), 1.96*pl.array(plotting[p]['prior'].std(1)), color='k', capsize=0, elinewidth=.5, label='Prior error', alpha=.5)
        plot(plotting[p]['prior'].mean(1), 'k', linewidth=2, label='Prior')
        plot(plotting[p]['pred'].mean(0), 'r', linewidth=2, label=plotting[p]['name'])
        ui = mc.utils.hpd(plotting[p]['pred'], .05)
        plot(ui[:,0], 'r--', linewidth=1)
        plot(ui[:,1], 'r--', linewidth=1)

        axis([0,100,0,ymax])
        legend(loc='upper left', title='%s: \nmare=%s, \ntime=%ss'%(plotting[p]['name'], round(plotting[p]['mare'],2), plotting[p]['time']))

    subplot(1,2,1)
    for p in range(len(plotting)):
        plot(plotting[p]['pred'].mean(0), linewidth=2, label=plotting[p]['name'])   
        if p == len(plotting)-1: 
            dismod3.graphics.plot_data_bars(model.get_data(data_type))
            axis([0,100,0,.12])
            legend(loc='upper left')
            title('%s (%s), %s %s %s'%(name, disease, country, sex, year))
    