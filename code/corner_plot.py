from numpy import *
import corner
from matplotlib import *
use('Agg')
#import matplotlib

from matplotlib.pyplot import *
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True #Guarantee type 1 plots to 
rcParams['text.usetex'] = True #Very important to force python to recognize Latex
rcParams['legend.numpoints']=1 
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}
rc('font', **font)


parspace=genfromtxt('../output/LAEBayes_mcmc_dM_1deg.out')
parprob=genfromtxt('../output/LAEBayes_mcmc_dM_chi_1deg.out')
parall=genfromtxt('../output/focc_1deg.txt')
parspace[:,1]+=parspace[:,0]

arg_ini=0
fig = corner.corner(parspace[arg_ini:],weights=exp(parprob[arg_ini:]),bins=20, smooth=0.6, labels=[r"$\log M_{\rm min}\left[\rm{M_{\odot}h^{-1}}\right]$", r"$\log M_{\rm max}\left[\rm{M_{\odot}h^{-1}}\right]$"],quantiles=[0.16, 0.5, 0.84])
fig.savefig('../output/likelyplot_1deg.png',  format='png', dpi=150,bbox_inches='tight')
fig = corner.corner(transpose([parall[:,0][arg_ini:],parall[:,2][arg_ini:]]),weight=exp(parprob)[ :len( parall[:,0] ) ][arg_ini:],bins=20, smooth=0.6, labels=[r"$\log M_{\rm min}\left[\rm{M_{\odot}h^{-1}}\right]$", r"$\log f_{\rm occ}$"],quantiles=[0.16, 0.5, 0.84])
fig.savefig('../output/likelyplot_1degf.png', format='png', dpi=150,bbox_inches='tight')
fig = corner.corner(transpose([parall[:,0][arg_ini:], (parall[:,1]-parall[:,0])[arg_ini:] ]), weight=exp(parprob)[ :len( parall[:,0] ) ][arg_ini:], bins=20, smooth=0.6, labels=[r"$\log M_{\rm min}\left[\rm{M_{\odot}h^{-1}}\right]$", r"$\Delta \log M\left[\rm{M_{\odot}h^{-1}}\right]$"],quantiles=[0.16, 0.5, 0.84])
fig.savefig('../output/likelyplot_1degdM.png', format='png', dpi=150,bbox_inches='tight')

fig = corner.corner(transpose([ (parall[:,1]-parall[:,0])[arg_ini:], parall[:,2][arg_ini:] ]), weight=exp(parprob)[ :len( parall[:,0] ) ][arg_ini:], bins=20,  smooth=0.6,labels=[r"$\Delta \log M\left[\rm{M_{\odot}h^{-1}}\right]$", r"$\log f_{\rm occ}$"],quantiles=[0.16, 0.5, 0.84])
fig.savefig('../output/likelyplot_1degdMf.png', format='png', dpi=150,bbox_inches='tight')


font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}
rc('font', **font)



fig = corner.corner(transpose([parall[:,0][arg_ini:],parall[:,1][arg_ini:],parall[:,1][arg_ini:]-parall[:,0][arg_ini:],parall[:,2][arg_ini:]]),weight=exp(parprob)[ :len( parall[:,0] ) ][arg_ini:],bins=20,  smooth=0.6,labels=[r"$\log M_{\rm min}\left[\rm{M_{\odot}h^{-1}}\right]$", r"$\log M_{\rm max}\left[\rm{M_{\odot}h^{-1}}\right]$", r"$\Delta \log M\left[\rm{M_{\odot}h^{-1}}\right]$",r"$\log f_{\rm occ}$"],quantiles=[0.16, 0.5, 0.84])
fig.savefig('../output/likelyplot_1deg_total.png', format='png', dpi=150,bbox_inches='tight')




font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 24}
rc('font', **font)

rcParams['figure.figsize']=8,8
figure()
fig1=scatter(parall[400:,0],parall[400:,4]-parall[400:,3],c=(parall[400:,1]-parall[400:,0]), cmap=cm.winter_r,vmin=0.1,vmax=1.5,s=15,marker='D')
fig2=scatter(parall[400:,0],parall[400:,6]-parall[400:,5],c=(parall[400:,1]-parall[400:,0]), cmap=cm.autumn_r,vmin=0.1,vmax=1.5,s=15,marker='o')
colorbar(fig1)
colorbar(fig2)
xlabel(r'$M_{\rm min}\left[\rm{M_{\odot}h^{-1}}\right]$')
ylabel(r'$\Delta \log \left( N_{\rm halos}/\tilde{N}_{\rm halos}\right)$')
savefig('../output/mmin_dfocc1.png', format='png', dpi=150,bbox_inches='tight')


figure()
scatter(parall[400:,0],parall[400:,6]-parall[400:,5],color='k')
xlabel(r'$\log M_{\rm min}\left[\rm{M_{\odot}h^{-1}}\right]$')
ylabel(r'$\log\left(f_{\rm occ}\left(97.5\%\right)/f_{\rm occ}\left(2.5\%\right)\right)$')
savefig('../output/mmin_dfocc2.png', format='png', dpi=150,bbox_inches='tight')



font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 24}
rc('font', **font)


rcParams['figure.figsize']=14,8
figure()
#fig2=scatter(parall[:,0],parall[:,11],c=(parall[:,1]), cmap=cm.autumn_r,vmin=0.1,vmax=1.5,s=15)
#fig1=scatter(parall[:,0],parall[:,9],c=(parall[:,1]), cmap=cm.winter_r,vmin=0.1,vmax=1.5,s=15)
#fig=scatter(parall[:,0],parall[:,7],c=(parall[:,1]), cmap=cm.gray_r,vmin=0.1,vmax=1.5,s=15)


fig2=scatter(parall[:,0],parall[:,11],c=(parall[:,1]-parall[:,0]), cmap=cm.autumn_r,vmin=0.1,vmax=1.5,s=20,marker='s')
fig1=scatter(parall[:,0],parall[:,9],c=(parall[:,1]-parall[:,0]), cmap=cm.winter_r,vmin=0.1,vmax=1.5,s=15,marker='D')
fig=scatter(parall[:,0],parall[:,7],c=(parall[:,1]-parall[:,0]), cmap=cm.gray_r,vmin=0.1,vmax=1.5,s=15,marker='o')


colorbar(fig2)
colorbar(fig1)
colorbar(fig)
plot([9,13],[9,13],color='g')
xlim(9.2,12.5)
ylim(9.2,12.5)
xlabel(r'$\log M_{\rm min}\left[\rm{M_{\odot}h^{-1}}\right]$')
ylabel(r'$M_{\rm \%}\left[\rm{M_{\odot}h^{-1}}\right]$')
savefig('../output/mmin_mmed_colordm.png', format='png', dpi=150,bbox_inches='tight')

th,co,dco=np.genfromtxt('../obs/ACF/Bielby2015.txt',unpack=1)
th,c1,dcl1,dcu1=np.genfromtxt('../output/corr105_d025.txt',unpack=1)
th,c2,dcl2,dcu2=np.genfromtxt('../output/corr105_d05.txt',unpack=1)
th,c3,dcl3,dcu3=np.genfromtxt('../output/corr105_d10.txt',unpack=1)
th,c4,dcl4,dcu4=np.genfromtxt('../output/corr105_d20.txt',unpack=1)
th,c5,dcl5,dcu5=np.genfromtxt('../output/corr105_d30.txt',unpack=1)

dcl1[dcl1>=c1] = c1[dcl1>=c1]*.999999
dcl2[dcl2>=c2] = c2[dcl2>=c2]*.999999
dcl3[dcl3>=c3] = c3[dcl3>=c3]*.999999
dcl4[dcl4>=c4] = c4[dcl4>=c4]*.999999
dcl5[dcl5>=c5] = c5[dcl5>=c5]*.999999

rcParams['figure.figsize']=12,6
figure();errorbar(th,c1,np.array([dcl1,dcu1]),label=r'$M_{\rm min}=10.5$'+' '+r'$\Delta M=0.25$')
errorbar(th,c2,np.array([dcl2,dcu2]),label=r'$\Delta M=0.5$')
errorbar(th,c3,np.array([dcl3,dcu3]),label=r'$\Delta M=1.0$')
errorbar(th,c4,np.array([dcl4,dcu4]),label=r'$\Delta M=2.0$')
errorbar(th,c5,np.array([dcl5,dcu5]),label=r'$\Delta M=3.0$')
errorbar(th,co,dco,label='ACF obs. Bielby+2016')
xlabel(r'$\theta [\rm arcmin]$')
ylabel(r'$\omega(\theta)$')
xlim(0.5,4.0)
ylim(0.006,30.0)
legend(fontsize=12)
xscale('log')
yscale('log')
savefig('../output/corr.png', format='png', dpi=150,bbox_inches='tight')
show()




