#bin/#!/usrt/env python
import numpy as np 
#import matplotlib.pyplot as P
#import correlation_lib as corr
import sys
import itertools
import emcee
import time
from pylab import *
tini=time.time()



load_snap_file=1 #load the snapshot positions and massess. Default 0.
indexes_subcat=1 #find the indexes for the subcatalogs. Default 1
plot_ndist=0  # plot the number halo distrubition in the subacatalog. Default 0

niter=100 #number of times to iter emcee
ndim=3 #number of parameters, 3 in this time: mmin,mma,focc.
nwalkers=16 # number of independent walks
ncores=16 #to optimize it should be set to nwalkers

#The center and widths of the seeds for EMCEE 
Mo=11.0 # 
dmmin=1.0
dmmax=1.0
fo=0.5
dfo=0.5

#Below are the seeds to be passed to EMCEE, Here I generete them randomly using a gaussian distribution.

mmin=Mo+dmmin*randn(nwalkers) 
mmin=np.abs(mmin)
mmax=Mo+0.2+dmmax*randn(nwalkers)
mmax=np.abs(mmax)
Focc=fo+dfo*randn(nwalkers)
Focc=np.abs(Focc)    

#initial seeds for the emcee walkers
pos = [[mmin[i],mmax[i],Focc[i]] for i in range(nwalkers)]

simwidth=250.0 # physical width of the simulation in Mpc/h
z=3.1
Nlae_mean=199 # The (mean) number of LAEs in the observational  catalog(s)
Dc=6403.7 #Mpc 
Da=1561.9 #Mpc 
Dl=26255.0 #Mpc
h=0.70 

Dc=Dc*h #convertion to Mpc/h
Da=Da*h
Dl=Dl*h

omegam=0.30711
omegal=0.69289

scale=7.855 #kpc/".
scale_arcmin=scale*60/1000 # Mpc/(')
xwidth=34 #width in arcmin of the observational catalog
ywidth=27 #arcmin
xwidth*=scale_arcmin*(1+z)*h #converting the size to the simulation units:  Mpc/h
ywidth*=scale_arcmin*(1+z)*h #converting the size to Mpc/h
zwidth=58.5*h #converting the size to Mpc/h. 58.5 is taken from the difference in comoving distance between the edges of the LAE filter


#The number of subcatalogs within the volume of the simulation is is nx*ny*nz
nx=np.int(np.floor(simwidth/xwidth))
ny=np.int(np.floor(simwidth/ywidth))
nz=np.int(np.floor(simwidth/zwidth))




def DD_histogram(X,Y,distance,th_min,th_max,theta_bins,logtheta=False):
    '''
    DD histogram
    Computes the angular two point distribution  in the observational catalog
    '''

    n_points=len(X)
    d_arr=[]
    #print "number of LAEs=" + str(n_points)
    for i in range(n_points):
        
        for j in range(i+1,n_points):
            d=(X[i] - X[j])*(X[i] - X[j]) + (Y[i] - Y[j])*(Y[i] - Y[j])
            d=np.sqrt(d)/distance
            d_arr=np.append(d_arr,d)
            
    
    theta=206265.0*d_arr
    #print d_arr, theta, "darr-theta"
    if logtheta: theta=np.log10(theta)
    
    DD,bins=np.histogram(theta,bins=theta_bins, range=(th_min,th_max))
    
    N=1.0*n_points*(1.0*n_points-1.0)/2.0
    DD=DD/N
    
    return DD,bins

'''
RR histogram
Computes the angular two point distribution  in the random catalog
'''


def RR_histogram(Xr,Yr,distance,th_min,th_max,theta_bins,cat_number=3,logtheta=False):
    
    n_points=len(Xr)
    #Xmin=1.0*np.floor( np.amin(X) )
    #Xmax=1.0*np.ceil( np.amax(X) )
    #Ymin=1.0*np.floor( np.amin(Y) )
    #Ymax=1.0*np.ceil( np.amax(Y) )
    

    #Xr= Xmin +  ( Xmax - Xmin )*np.random.random_sample(n_points*cat_number)
    #Yr=Ymin +   ( Ymax - Ymin )*np.random.random_sample(n_points*cat_number)
    d_arr=[]

    for m in range(cat_number):
        
        for i in range(n_points):
            for j in range(i+1,n_points):
                
                d=(Xr[i + n_points*m ] - Xr[j + n_points*m ])*(Xr[i + n_points*m ] - Xr[j + n_points*m ]) + (Yr[i + n_points*m ] - Yr[j + n_points*m ])*(Yr[i + n_points*m ] - Yr[j + n_points*m ])
                d=np.sqrt(d)/distance
                d_arr=np.append(d_arr,d)
                
    theta=206265*d_arr 
    if logtheta: theta=np.log10(theta)
    RR,bins=np.histogram(theta,bins=theta_bins, range=(th_min,th_max))
    
    N=1.0*n_points*(1.0*n_points-1.0)/2.0
    N=N*cat_number
    RR=RR/N            
    
    return RR,bins



'''
DR histogram
Computes the angular two point distribution  betweeen the observational and the random catalog
'''


def DR_histogram(X,Y,Xr,Yr,distance,th_min,th_max,theta_bins,cat_number=3,logtheta=False):
    n_points=len(X)
    #Xmin=1.0*np.floor( np.amin(X) )
    #Xmax=1.0*np.ceil( np.amax(X) )
    #Ymin=1.0*np.floor( np.amin(Y) )
    #Ymax=1.0*np.ceil( np.amax(Y) )
    

    #Xr= Xmin +  ( Xmax - Xmin )*np.random.random_sample(n_points*cat_number)
    #Yr=Ymin +   ( Ymax - Ymin )*np.random.random_sample(n_points*cat_number)
    d_arr=[]
    
    for m in range(cat_number): 
        for i in range(n_points):
            for j in range(n_points):
                d=( X[i] - Xr[j + n_points*m] )*( X[i] - Xr[j + n_points*m] ) + ( Y[i] - Yr[j + n_points*m] )*( Y[i] - Yr[j + n_points*m] )
                d=np.sqrt(d)/distance
                d_arr=np.append(d_arr,d)
    theta=206265*d_arr 
    if logtheta: theta=np.log10(theta)
    
    DR,bins=np.histogram(theta,bins=theta_bins, range=(th_min,th_max))
    
    N=1.0*n_points*(n_points)
    N=(N*cat_number)            
    DR=DR/N
    
    return DR,bins

# Different Angular correlatio function approaches
def landy_correlation(DD,RR,DR):
    CORR= (DD - 2.0*DR + RR)/RR
    return CORR
def peebles_correlation(DD,DR):
    CORR=DD/DR - 1.0
    return CORR
def standard_correlation(DD,RR):
    CORR=DD/RR - 1.0
    return CORR 



def randsample(arraySize,sampleSize):
#taken from stackoverflow forum: how do i create a LIST of unique random numbers
# This function takes a random subsample of length SampleSize from an sample of lenght arraySize 
    answer = set()
    answerSize = 0

    while answerSize < sampleSize:
        r = np.random.randint(0,arraySize)
        if r not in answer:
            answerSize += 1
            answer.add(r)
    answer=np.array([int(x) for x in answer])
    return answer





# Here I am reading the position an mass of the halos in the snapshot
#I have previosly taken from the main snapshot file the relevant columns into a single file
# by means of the command line using: awk '{print $1}' snapfile > posx.txt .Here $1 represent the position of the desired column in the snapshot file 

if load_snap_file:
    filename="../data/FOFhalos/posx.txt"
    x=np.genfromtxt(filename)
    filename="../data/FOFhalos/posy.txt"
    y=np.genfromtxt(filename)
    
    filename="../data/FOFhalos/posz.txt"
    z=np.genfromtxt(filename)
    
    filename="../data/FOFhalos/mass.txt"
    mass=np.genfromtxt(filename)
    
    cat_path="../data/FOFhalos/subcatalogs/"

def subcatalogs(x,y,z,mass,savepath="../data/FOFhalos/subcatalogs/"):
    """ This function returns a list of the arrays of the indexes (index) and 
    an array of the lenghts of the subarrays
    with the size of the observational catalog (n_array)"""
    first=1
    n_array=[]
    print "creating sub-catalogs"
    for i in range(nx):
        xmin=i*xwidth
        xmax=(i+1)*xwidth
        wxmin=x>xmin
        wxmax=x<xmax
        for j in range(ny):
            ymin=j*ywidth
            ymax=(j+1)*ywidth
            wymin=y>ymin
            wymax=y<ymax
            for k in range(nz):
                zmin=k*ywidth
                zmax=(k+1)*ywidth
                wzmin=z>zmin
                wzmax=z<zmax
                where=wxmin*wxmax*wymin*wymax*wzmin*wzmax
                n_array.append( len(x[where]) )
                if first==1:index=[np.array(range(len(x)))[where]];first=0
                else:index=index+[np.array(range(len(x)))[where]]
    print "finishing sub-catalog creation"
    return index, n_array
                
if indexes_subcat:
    index,n_array=subcatalogs(x,y,z,mass,savepath=cat_path)

cat_path="../data/FOFhalos/subcatalogs/"

# This is to plot the distrubion of the number of halos in the subcatalog in different mass bins
# plot_ndist is set to 0 by default
if plot_ndist:
    
    mmin=np.arange(10.0,12.0,0.5)
    mmax=np.arange(10.5,12.5,0.5)

    figure()
    n, bins, patches = plt.hist(np.log10(n_array), np.sqrt(len(n_array))/2.0, normed=True,histtype='step', cumulative=False,label='all')
    file=open("mminmax_density.txt","w")
    file.write("")
    file.close()
    file=open("../output/mminmax_density.txt","a")

    for Mmin,Mmax in itertools.product(mmin,mmax):
        if Mmin>=Mmax-0.00001: continue
        Nmass=[]
        #print Mmin,Mmax
        for i,j,k in itertools.product(range(nx),range(ny),range(nz)):
            filename=cat_path+"catalog_"+str(i)+"_"+str(j)+"_"+str(k)+".txt"
            x,y,z,m=np.genfromtxt(filename,unpack=1)
            m=np.log10(m)
            wmin=m>Mmin; wmax=m<Mmax
            wm=wmin*wmax
            xm=x[wm]; ym=y[wm]; zm=z[wm]; mm=m[wm]
            Nmass.append(len(xm))
            del(xm);del(ym);del(zm);del(mm)
            print i,j,k
            n, bins, patches = plt.hist(np.log10(Nmass), np.sqrt(len(Nmass))/2.0, normed=True,histtype='step', cumulative=False,label="Mmin="+str(Mmin)+" Mmax="+str(Mmax))
            Nmass=np.append(np.array([mmin,mmax]), Nmass)
            np.savetxt(file,Nmass)
    plt.legend()
    file.close()
    del(file)
    te=time.time()-tini-t1-t2



def corelationlike(Mmin,Mmax,focc,x,y,z,mass,index,n_array,xwidth,ywidth,zwidth,Dc=Dc,Da=Da,Dl=Dl,cat_number=1,logtheta=1,convert_to_arcsec=1,estimator='landy',acf_file='../data/obs/ACF/Bielby2015.txt',Nlae=Nlae_mean):
    """ This function computes the likelyhood of a model with Mmin, Mmax, focc to macth the observational ACF and mean number density of galaxies.
    It is fundamental for the EMCEE code.""" 
    print "Mmin=",Mmin," Mmax=",Mmax," focc=",focc
    theta,corrobs,dcorr=np.genfromtxt(acf_file,unpack=1)
    if convert_to_arcsec: theta=theta*60.0
    if logtheta: theta=np.log10(theta)
    
    theta_bins=theta-(theta[1]-theta[0])/2.0
    theta_bins=np.append(theta_bins,theta[-1]+(theta[1]-theta[0])/2.0)
    th_min=theta_bins[0];th_max=theta_bins[-1]
    corr=np.zeros( len(theta)  )
    m=np.log10(mass)
    if Mmin>=Mmax-0.00001 or focc>=1 or focc<0 or Mmin<9.5 or Mmax<9.6 or  Mmax>13.0 or Mmin>12.5: chi2=np.inf;print "chi2=inf"
    else:
        print th_max, th_min, "thmax-min"
        count=1
        #for i,j,k in itertools.product(range(nx),range(ny),range(nz)):
        
        x1=x[index[count]];y1=y[index[count]];m1=m[index[count]]
        wmin=m1>Mmin; wmax=m1<Mmax
        wm=wmin*wmax
        print x1, y1, wm
        #raw_input("Press Enter to continue...")
        xm=x1[wm]; ym=y1[wm]; n_points1=len(xm); n_points=len(xm)
        Nhalo=n_points*mean(n_array)/n_array[0] # a crude stimation of the mean number of halos with Mmin<m<Mmax in a box
        focch=Nlae/Nhalo
        Nlae_teo=focc*Nhalo
        dNlae=np.sqrt(Nlae)
        if  (1.0*Nlae/Nlae_teo)>6.3 or (1.0*Nlae_teo/Nlae)>6.3:
            chi2=np.inf
        else:
            count=0
            corr=np.zeros( len(theta)  )
            narray=[]
            for i,j,k in itertools.product(range(nx),range(ny),range(nz)):
                x1=x[index[count]];y1=y[index[count]];m1=m[index[count]]
                wmin=m1>Mmin; wmax=m1<Mmax
                wm=wmin*wmax
                xm=x1[wm]; ym=y1[wm]; n_points1=len(xm);n_points=len(xm)
                if Nlae<=n_points:
                    rpos=randsample(n_points,Nlae)
                    xm=xm[rpos];ym=ym[rpos]
                    n_points=Nlae
                
                narray=np.append(narray,n_points1)
                DD,bins=DD_histogram(xm,ym,Dc,th_min,th_max,theta_bins,logtheta=logtheta)
                xr= xwidth*(i+np.random.random_sample(n_points*cat_number))
                yr= ywidth*(j+np.random.random_sample(n_points*cat_number))
                RR,bins=RR_histogram(xr,yr,Dc,th_min,th_max,theta_bins,cat_number=cat_number,logtheta=logtheta)
                DR,bins=DR_histogram(xm,ym,xr,yr,Dc,th_min,th_max,theta_bins,cat_number=cat_number,logtheta=logtheta)
                #print "DD=",DD," RR=",RR," DR=",DR
                if estimator=='landy':corr=corr+landy_correlation(DD,RR,DR)
                elif estimator=='peebles':corr=corr+peebles_correlation(DD,DR)
                elif estimator=='standard':  corr=corr+standard_correlation(DD,RR)
                else:corr=corr+landy_correlation(DD,RR,DR)
                count=count+1
            corr=np.ma.masked_invalid(corr/count)
            Nhalo=np.mean(n_array)
            focch=Nlae/Nhalo
            Nlae_teo=focc*Nhalo
            dNlae=np.sqrt(Nlae)
        chi2=((Nlae_teo-Nlae)*(Nlae_teo-Nlae)/(dNlae*dNlae))+np.sum( (corr-corrobs)*(corr-corrobs)/(dcorr*dcorr) )/(np.ma.count(corr)-1) 
    print "chi2",chi2
    print "corr=",corr,(np.ma.count(corr))
    print "corrobs=",corrobs
    #raw_input("Press Enter to continue...")
    return -1*chi2




def lnlike(params,x,y,z,mass,index,n_array,xwidth,ywidth,zwidth):
    """This is the way the likelihood function should be redefined to be  passed to emcee"""
    Mmin,Mmax,focc=params
    return corelationlike(Mmin,Mmax,focc,x,y,z,mass,index,n_array,xwidth,ywidth,zwidth)


#Setting emcee 6
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike, args=(x,y,z,mass,index,n_array,xwidth,ywidth,zwidth),threads=ncores)

# Runing emcee for a few times 
pos0, prob, state = sampler.run_mcmc(pos, 4)

sampler.reset()

fn = "../output/LAEBayes_mcmc_new.out"
f = open(fn, "w")
f.close()

# Restarting  emcee over a large loop of niter times
for pos, prob, rstate in sampler.sample(pos0, prob, state, iterations=niter):
    # Write the current position to a file, one line per walker                                                                                                                     
    f = open(fn, "a")
    f.write("\n".join(["\t".join([str(q) for q in p]) for p in pos]))
    f.write("\n")
    f.close()
