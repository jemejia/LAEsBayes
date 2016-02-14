from matplotlib import *
from numpy import *
use('Agg')
import corner
import matplotlib

import matplotlib.pyplot as P
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True #Guarantee type 1 plots to 
rcParams['text.usetex'] = True #Very important to force python to recognize Latex
rcParams['legend.numpoints']=1 
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}
rc('font', **font)


parspace=genfromtxt('../output/LAEBayes_mcmc_dM_final.out')
parspace[:,1]+=parspace[:,0]

fig = corner.corner(parspace[200:],bins=10, labels=["$M_{min}$", "$M_{max}$"],quantiles=[0.16, 0.5, 0.84])
fig.savefig('../output/likelyplot.png', format='png', dpi=600,bbox_inches='tight')
