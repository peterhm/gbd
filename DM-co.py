Initialization

pwd

cd /homes/peterhm/gbd/book

import book_graphics
reload(book_graphics)

cd /homes/peterhm/gbd/

import sys
sys.path += ['../book', '../dm3-computation_only/', '../dm3-computation_only/book']
import pylab as pl
import pymc as mc
import pandas
 
import dismod3
reload(dismod3)


Functions

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

# find disease number from disease name
def disease_num(name):
    bm_path = '/snfs1/Project/GBD/dalynator/yld/best_models.csv'
    bm_csv = pandas.read_csv(bm_path,index_col=None)
    dismod_models = bm_csv.groupby('dismod_model_number').apply(lambda df: df.ix[df.index[0], 'outcome_name'])
    dismod_models = dismod_models.drop([0], axis=0)
    dismod_models = dict(sort(dismod_models))
    # change Series object into dictionary
    from collections import defaultdict
    reversed_dict = defaultdict(list)
    for key,value in dismod_models.iteritems():
        reversed_dict[value].append(key)
    num = reversed_dict[name]
    if num == []:
        print 'No DisMod-MR estimates for %s'%name
    elif len(num) > 1:
        print 'DisMod-MR has more than one model for %s'%name
        num = [int(k) for k in num]
    else:
        num = int(num[0])
    return num

# find country name

# create disease model with relavtive data
def load_new_model(disease, country, sex, year):
    model = dismod3.data.load('/home/j/Project/dismod/output/dm-%s'%disease)
    # keep relative data
    model.keep(areas=[country], sexes=[sex])
    
    # remove covariates
    model.input_data = model.input_data.drop(['x_LDI_id_Updated_7July2011', 'x_ihme_health_system_access_19jul2011'], 1)
    
    # remove increasing age pattern, to make data speak for itself
    #model.parameters['i']['increasing'] = {'age_end': 0, 'age_start': 0}
    #model.parameters['i']['smoothness']['amount'] = 'Moderately'
    return model

def get_emp(disease, country, sex, year):
    # load posterior estimates from GBD 2010 Study
    global_model = load_new_model(disease, 'all', sex, year)
    emp = pandas.read_csv('/home/j/Project/dismod/output/dm-%s/posterior/dm-%s-%s-%s-%s-%s.csv'%(disease, disease, full_name[data_type], global_model.hierarchy.in_edges(country)[0][0], sex, year), index_col=None)
 
    # remove population numbers
    del emp['Population']
 
    # keep only estimates from country
    cty_ix = (emp['Iso3'] == country)
    emp = emp[cty_ix]
    del emp['Iso3']
 
    # keep only estimates for data type
    try:
        assert pl.all(emp['Rate type'] == full_name[data_type])
    except:
        dt_ix = (emp['Rate type']) == full_name[data_type]
        emp = emp[dt_ix]
    del emp['Rate type'] 
    
    # return GBD 2010 Study posterior 
    emp.index = emp['Age']
    del emp['Age']
    return emp.mean(1), emp.std(1)*10

# set fixed effects and random effects from strength-gathering posterior mean
#model.parameters[t]['random_effects'] = ...
#model.parameters[t]['fixed_effects'] = ...
 
def find_fnrfx(t, model, vars):
    for n, col in zip(vars['alpha'], vars['U'].columns):
        stats = n.stats()
        model.parameters[t]['random_effects'][col] = dict(dist='Constant', mu=stats['mean'], sigma=stats['standard deviation'])
 
    # also save empirical prior on sigma_alpha, the dispersion of the random effects
    for n in vars['sigma_alpha']:
        stats = n.stats()
        model.parameters[t]['random_effects'][n.__name__] = dict(dist='TruncatedNormal', mu=stats['mean'], sigma=stats['standard deviation'], lower=.01, upper=.5)
        
    for n, col in zip(vars['beta'], vars['X'].columns):
        stats = n.stats()
        model.parameters[t]['fixed_effects'][col] = dict(dist='Constant', mu=stats['mean'], sigma=stats['standard deviation'])
  
Choose disease, country, age, sex, year

disease = 39098
country = 'FIN'
sex = 'male'
year = 2005
data_type = 'i'

model = load_new_model(disease, country, sex, year)

mu_prior, sigma_prior = get_emp(disease, country, sex, year)

find_fnrfx('i', fin_wp, we.vars['i'])

result = dismod3.data.ModelVars()
result['i'] = dismod3.data_model.data_model('i', fin_wp, 'i',
                                            'FIN', 'male', 2005,
                                            None, fin_mu_prior, fin_sigma_prior,
                                            rate_type='neg_binom')
fin_wp.vars += result

%time dismod3.fit.fit_asr(fin_wp, 'i', iter=iter, thin=thin, burn=burn)
