""" Tools for summarizing DM-MR results"""

import pandas as pd
population = pd.read_csv('/home/j/Project/dismod/age_weights_final_20130128.csv', index_col='age')

# make rows 0, .01, .1 into a single age group
weight = population.to_dict()['weight']
weight[0.0] = weight[0.0]+weight[0.01]+weight[0.1]
weight.pop(0.01)
weight.pop(0.1)

population = pd.DataFrame({'weight': weight})

population.index = range(len(population.index))
age_index = {i:(population['age_end']<=i).sum() for i in range(101)}
age_group = {i: population.age_start[age_index[i]] for i in range(101)}
population.index = population.age_start

def prev_draws(df):
    age_specific_cases = df.filter(like='Draw').mul(df.Population, axis='index')
    age_specific_prev = age_specific_cases.sum() / df.Population.sum()
        
    return age_specific_prev

def age_std_prev(df):
    g = df.groupby('age_group')
    age_specific_prev_draws = g.apply(prev_draws)
    age_std_draws = age_specific_prev_draws.mul(population.weight, axis='index').sum()
    return age_std_draws.describe(95).ix[['mean', '2.5%', '97.5%']].to_dict()

def calc_age_std_prevs(id):
    results = []
    print id,
    for region in ['asia_central', 'asia_east', 'asia_pacific_high_income',
                   'asia_south', 'asia_southeast', 'australasia', 'caribbean',
                   'europe_central', 'europe_eastern', 'europe_western',
                   'latin_america_andean', 'latin_america_central',
                   'latin_america_southern', 'latin_america_tropical',
                   'north_africa_middle_east', 'north_america_high_income', 'oceania',
                   'sub-saharan_africa_central', 'sub-saharan_africa_east',
                   'sub-saharan_africa_southern', 'sub-saharan_africa_west']:
        print '.',
        for year in [1990, 2005, 2010]:
            for sex in ['male', 'female']:
                df_ysr = pd.read_csv('/home/j/Project/dismod/output/dm-%d/posterior/dm-%d-prevalence-%s-%s-%s.csv'%(id, id, region, sex, year))
                df_ysr['age_group'] = df_ysr['Age'].map(age_group)
                result = age_std_prev(df_ysr)
                result.update(year=year,sex=sex,region=region)
                
                results.append(result)
    print

    results = pd.DataFrame(results).filter(['region', 'year', 'sex', 'mean', '2.5%', '97.5%'])
    results.to_csv('/home/j/Project/dismod/age_std_prev_%d.csv'%id)
