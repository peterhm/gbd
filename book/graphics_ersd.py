import sys
sys.path += ['../gbd', '../gbd/book', '../dm3-computation_only/', '../dm3-computation_only/book']
import pylab as pl
import pymc as mc
import pandas

import dismod3
reload(dismod3)

import book_graphics
reload(book_graphics)
import matplotlib as mpl

# set font
book_graphics.set_font()

def my_axis(ymax):
    pl.axis([-5,105,-ymax/10.,ymax])

def load_new_model():
    try:
        model = dismod3.data.load('/home/j/Project/dismod/output/dm-32015')
    except:
        model = dismod3.data.load('/home/j/Project/dismod/dismod_status/prod/dm-32015')
    model.keep(sexes=['male', 'total'], start_year=2005, end_year=2005)
    model.parameters['m_with'] = {}
    model.input_data = model.input_data.drop(['x_LDI_id_Updated_7July2011'],1)
    model.parameters['i']['fixed_effects']['x_sex'] = dict(dist='Normal', mu=0, sigma=.0001)
    return model

# figure ckd-data
incon_i = load_new_model()
incon_i.keep(areas=['australasia'])

pl.figure(**book_graphics.full_plus_page_params)

pl.subplot(2,2,1)
dismod3.graphics.plot_data_bars(incon_i.get_data('p'))
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .0007, .0014, .0021, .0028], [0, 0.7, 0.14, 0.21, 0.28])
my_axis(.003)
book_graphics.subtitle('(a)')


pl.subplot(2,2,2)
dismod3.graphics.plot_data_bars(incon_i.get_data('i'))
pl.xlabel('Age (years)')
pl.ylabel('Incidence \n (per 10,000 PY)'+'\n\n', ha='center')
pl.yticks([0, .0002, .0004, .0006, .0008], [0, 2, 4, 6, 8])
my_axis(.0009)
book_graphics.subtitle('(b)')


pl.subplot(2,2,3)
dismod3.graphics.plot_data_bars(incon_i.get_data('r'))
pl.xlabel('Age (years)')
pl.ylabel('Remission (per 100 PY)')
pl.yticks([0, .04, .08, .12, .16], [0, 4, 8, 12, 16])
my_axis(.19)
book_graphics.subtitle('(c)')


pl.subplot(2,2,4)
dismod3.graphics.plot_data_bars(incon_i.get_data('m_with'))
pl.xlabel('Age (years)')
pl.ylabel('With-condition mortality \n (per 100 PY)'+'\n\n', ha='center')
pl.yticks([0, .1, .2, .3, .4], [0, 10, 20, 30, 40])
my_axis(.45)
book_graphics.subtitle('(d)')


pl.subplots_adjust(top=.99, bottom=.14, wspace=.35, hspace=.25)


pl.savefig('book/graphics/ckd-data.pdf')
pl.savefig('book/graphics/ckd-data.png')

# figure ckd-incon_v_con
all_aus = load_new_model()
all_aus.keep(areas=['australasia'])

output = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_esrd.csv')

pl.figure(**book_graphics.full_plus_page_params)

pl.subplot(2,2,1)
dismod3.graphics.plot_data_bars(all_aus.get_data('p'), color='grey')
pl.plot(pl.arange(101), pl.array(output['c_p']), 'k-', linewidth=2, label='Compartmental')
pl.plot(pl.arange(101), pl.array(output['s_p']), 'k--', linewidth=2, label='Spline')
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .0007, .0014, .0021, .0028], [0, 0.7, 0.14, 0.21, 0.28])
my_axis(.003)
book_graphics.subtitle('(a)')


pl.subplot(2,2,2)
dismod3.graphics.plot_data_bars(incon_i.get_data('i'), color='grey') 
pl.plot(pl.arange(101), pl.array(output['c_i']), 'k-', linewidth=2, label='Compartmental')
pl.plot(pl.arange(101), pl.array(output['s_i']), 'k--', linewidth=2, label='Spline')
pl.xlabel('Age (years)')
pl.ylabel('Incidence \n (per 10,000 PY)'+'\n\n', ha='center')
pl.yticks([0, .0002, .0004, .0006, .0008], [0, 2, 4, 6, 8])
my_axis(.0009)
book_graphics.subtitle('(b)')


pl.subplot(2,2,3)
dismod3.graphics.plot_data_bars(incon_i.get_data('r'), color='grey')
pl.plot(pl.arange(101), pl.array(output['c_r']), 'k-', linewidth=2, label='Compartmental')
pl.plot(pl.arange(101), pl.array(output['s_r']), 'k--', linewidth=2, label='Spline')
pl.xlabel('Age (years)')
pl.ylabel('Remission (per 100 PY)')
pl.yticks([0, .04, .08, .12, .16], [0, 4, 8, 12, 16])
my_axis(.19)
book_graphics.subtitle('(c)')


pl.subplot(2,2,4)
dismod3.graphics.plot_data_bars(incon_i.get_data('m_with'), color='grey') 
pl.plot(pl.arange(101), pl.array(output['c_m']), 'k-', linewidth=2, label='Compartmental')
pl.plot(pl.arange(101), pl.array(output['s_m']), 'k--', linewidth=2, label='Spline')
pl.xlabel('Age (years)')
pl.ylabel('With-condition mortality \n (per 100 PY)'+'\n\n', ha='center')
pl.yticks([0, .1, .2, .3, .4], [0, 10, 20, 30, 40])
my_axis(.45)
book_graphics.subtitle('(d)')


pl.legend(loc='upper center', bbox_to_anchor=(-.2,-.2), fancybox=True, shadow=True, ncol=2)    
pl.subplots_adjust(top=.99, bottom=.14, wspace=.35, hspace=.25)

pl.savefig('book/graphics/ckd-incon_v_con.pdf')
pl.savefig('book/graphics/ckd-incon_v_con.png')

# figure ckd-m_with_smoothing
pl.figure(**book_graphics.full_page_params)

pl.plot(pl.array(output['c_m']), 'k-', linewidth=2, label='Compartmental')
pl.plot(pl.array(output['s_m']), 'k--', linewidth=2, label='Spline')
pl.plot(pl.array(output['s_m_smooth']), 'k:', linewidth=2, label='Spline with smoothing') 

pl.xlabel('Age (years)')
pl.ylabel('With-condition mortality (Per 100 PY)'+'\n\n', ha='center')
pl.yticks([0, .1, .2, .3, .4], [0, 10, 20, 30, 40])
my_axis(.45)
pl.legend(loc='upper left', fancybox=True, shadow=True)

pl.subplots_adjust(hspace=.35)
pl.subplots_adjust(wspace=.35)

pl.savefig('book/graphics/ckd-m_with_smoothing.pdf')
pl.savefig('book/graphics/ckd-m_with_smoothing.png')

# figure ckd-asp_scatter
scatter = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_esrd_regions.csv')

pl.figure(**book_graphics.full_page_params)

pl.plot(pl.array(scatter['compartmental']), pl.array(scatter['spline']), 'ks', label = 'GBD 2010 Region')
pl.plot([-1,1], [-1,1], 'k-')
pl.xlabel('Compartmental prevalence \nestimates (%)')
pl.ylabel('Spline prevalence \nestimates (%)')
pl.yticks([0, .0004, .0006, .0008, .001], [0, 0.04, 0.06, 0.08, 0.10])
pl.xticks([0, .0004, .0006, .0008, .001], [0, 0.04, 0.06, 0.08, 0.10])
pl.axis([.0003,.0011,.0003,.0011])

pl.savefig('book/graphics/ckd-asp_scatter.pdf')
pl.savefig('book/graphics/ckd-asp_scatter.png')

pl.show()