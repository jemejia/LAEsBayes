#bin/#!/usrt/env python
import numpy as np 
import os
#import matplotlib.pyplot as P
#import correlation_lib as corr
import sys
import itertools
import emcee
import time
import matplotlib
matplotlib.use('Agg')
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}
matplotlib.rc('font', **font)

from pylab import *
ioff()
try:
    os.mkdir('../output/')
except: print '../output already exist. Not a problem'
#tini=time.time()




load_snap_file=1 #load the snapshot positions and massess. Default 0.
indexes_subcat=1 #find the indexes for the subcatalogs. Default 1
plot_ndist=0  # plot the number halo distrubition in the subacatalog. Default 0
run_emcee=1 
compute_RR=1
ncatRR=1000


sufix='_1deg'
niter=400 #number of times to iter emcee
ndim=2 #number of parameters, 3 in this time: mmin,mma,focc.
nwalkers=24 # number of independent walks
ncores=24 #to optimize it should be set to nwalkers

#The center and widths of the seeds for EMCEE 
Mo=10.0 # 
dmmin=1.0
dmmax=1.0
fo=0.5
dfo=0.5

#Below are the seeds to be passed to EMCEE, Here I generete them randomly using a gaussian distribution.

mmin=Mo+dmmin*randn(nwalkers) 
mmin=np.abs(mmin)
dM=1.0*randn(nwalkers)
dM=np.abs(dM)
Focc=fo+dfo*randn(nwalkers)
Focc=np.abs(Focc)    

#initial seeds for the emcee walkers
if ndim==3:
    pos = [[mmin[i],dM[i],Focc[i]] for i in range(nwalkers)]
if ndim==2:
    pos = [[mmin[i],dM[i]] for i in range(nwalkers)]
simwidth=250.0 # physical width of the simulation in Mpc/h
redshift=3.1
Nlae_mean=643 # The (mean) number of LAEs in the observational  catalog(s)
Dc=6403.7 #Mpc 
Da=1561.9 #Mpc 
Dl=26255.0 #Mpc
h=0.70 

Dc=Dc*h #convertion to Mpc/h
Da=Da*h
Dl=Dl*h

omegam=0.30711
omegal=0.69289

scale=7.572 #kpc/".
scale_arcmin=scale*60/1000 # Mpc/(')
xwidth=62#34 #width in arcmin of the observational catalog
ywidth=62#27 #arcmin
xwidth*=scale_arcmin*(1+redshift)*h #converting the size to the simulation units:  Mpc/h
ywidth*=scale_arcmin*(1+redshift)*h #converting the size to Mpc/h
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
    to=time.time()
    n_points=len(X)
    #print "npoints=", n_points
    
    #print "number of LAEs=" + str(n_points)
    dx=np.array([(a-b)*(a-b) for a,b in itertools.product(X,X)])
    dy=np.array([(a-b)*(a-b) for a,b in itertools.product(Y,Y)])
    wx=dx!=0
    wy=dy!=0
    if len(dy[wy])<len(dx[wx]):
        dy=dy[wx]
        dx=dx[wx]
    else:
        dy=dy[wy]
        dx=dx[wy]

    d_arr=np.sqrt(dx+dy)/distance
    """
    d_arr=[]
    for i in range(n_points):
        for j in range(i+1,n_points):
            d=(X[i] - X[j])*(X[i] - X[j]) + (Y[i] - Y[j])*(Y[i] - Y[j])
            d=np.sqrt(d)/distance
            d_arr=np.append(d_arr,d)
    """
    t1=time.time()
    #print "tiempo armando d array=", t1-to
    try:
        theta=206265.0*d_arr
    except:
        theta=1000*th_max*np.ones(3)

    ##print d_arr, theta, "darr-theta"
    if logtheta: theta=np.log10(theta)

    t2=time.time()
    #print "tiempo def theta=", t2-t1

    
    DD,bins=np.histogram(theta,bins=theta_bins, range=(th_min,th_max))

    
    t3=time.time()
    #print "tiempo armando histo=", t3-t2

    
    #N=1.0*n_points*(1.0*n_points-1.0)/2.0
    N=1.0*n_points*(1.0*n_points-1.0)
    DD=DD/N
    
    return DD,bins

'''
RR histogram
Computes the angular two point distribution  in the random catalog
'''


def RR_histogram(Xr,Yr,distance,th_min,th_max,theta_bins,cat_number=3,logtheta=False):

    #to=time.time()
    n_points=len(Xr)
    #Xmin=1.0*np.floor( np.amin(X) )
    #Xmax=1.0*np.ceil( np.amax(X) )
    #Ymin=1.0*np.floor( np.amin(Y) )
    #Ymax=1.0*np.ceil( np.amax(Y) )
    

    #Xr= Xmin +  ( Xmax - Xmin )*np.random.random_sample(n_points*cat_number)
    #Yr=Ymin +   ( Ymax - Ymin )*np.random.random_sample(n_points*cat_number)
    dx=np.array([(a-b)*(a-b) for a,b in itertools.product(Xr,Xr)])
    dy=np.array([(a-b)*(a-b) for a,b in itertools.product(Yr,Yr)])
    wx=dx!=0
    wy=dy!=0
    if len(dy[wy])<len(dx[wx]):
        dy=dy[wx]
        dx=dx[wx]
    else:
        dy=dy[wy]
        dx=dx[wy]


    d_arr=np.sqrt(dx+dy)/distance
    

    """
    d_arr=[]
    for i in range(n_points):
        for j in range(i+1,n_points):
            
            d=(Xr[i] - Xr[j])*(Xr[i] - Xr[j]) + (Yr[i] - Yr[j])*(Yr[i] - Yr[j])
            d=np.sqrt(d)/distance
            d_arr=np.append(d_arr,d)
    """
    #t1=time.time()
    #print "tiempo armando d array=", t1-to
    
    
    try:
        theta=206265.0*d_arr
    except:
        theta=1000*th_max*np.ones(3)
    

    
    if logtheta: theta=np.log10(theta)

    
    #t2=time.time()
    #print "tiempo def theta=", t2-t1

    
    
    RR,bins=np.histogram(theta,bins=theta_bins, range=(th_min,th_max))
    
    #t3=time.time()
    #print "tiempo armando histo=", t3-t2

    

    
    
    #N=1.0*n_points*(1.0*n_points-1.0)/2.0
    N=1.0*n_points*(1.0*n_points-1.0)
    N=N*cat_number
    RR=RR/N            
    
    return RR,bins



'''
DR histogram
Computes the angular two point distribution  betweeen the observational and the random catalog
'''


def DR_histogram(X,Y,Xr,Yr,distance,th_min,th_max,theta_bins,cat_number=1,logtheta=False):
    #to=time.time()
    n_points=len(X)
    #Xmin=1.0*np.floor( np.amin(X) )
    #Xmax=1.0*np.ceil( np.amax(X) )
    #Ymin=1.0*np.floor( np.amin(Y) )
    #Ymax=1.0*np.ceil( np.amax(Y) )
    

    #Xr= Xmin +  ( Xmax - Xmin )*np.random.random_sample(n_points*cat_number)
    #Yr=Ymin +   ( Ymax - Ymin )*np.random.random_sample(n_points*cat_number)
    dx=np.array([(a-b)*(a-b) for a,b in itertools.product(X,Xr)])
    dy=np.array([(a-b)*(a-b) for a,b in itertools.product(Y,Yr)])
    wx=dx!=0
    wy=dy!=0
    if len(dy[wy])<len(dx[wx]):
        dy=dy[wx]
        dx=dx[wx]
    else:
        dy=dy[wy]
        dx=dx[wy]


    d_arr=np.sqrt(dx+dy)/distance
    
    """
    d_arr=[]
    
    for i in range(n_points):
        for j in range(n_points):
            d=( X[i] - Xr[j] )*( X[i] - Xr[j] ) + ( Y[i] - Yr[j] )*( Y[i] - Yr[j] )
            d=np.sqrt(d)/distance
            d_arr=np.append(d_arr,d)
    """
    #t1=time.time()
    #print "tiempo armando d array=", t1-to
        
    
    try:
        theta=206265.0*d_arr
    except:
        theta=1000*th_max*np.ones(3)

                


    if logtheta: theta=np.log10(theta)

    
    #t2=time.time()
    #print "tiempo def theta=", t2-t1

    
    
    DR,bins=np.histogram(theta,bins=theta_bins, range=(th_min,th_max))
    
    #t3=time.time()
    #print "tiempo armando histo=", t3-t2


    

    
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
    filename="../FOFhalos/posx.txt"
    x=np.genfromtxt(filename)
    filename="../FOFhalos/posy.txt"
    y=np.genfromtxt(filename)
    
    filename="../FOFhalos/posz.txt"
    z=np.genfromtxt(filename)
    
    filename="../FOFhalos/mass.txt"
    mass=np.genfromtxt(filename)
    
    cat_path="../FOFhalos/subcatalogs/"

try:
    os.mkdir(cat_path)
except:
    pass

def subcatalogs(x,y,z,mass,savepath="../FOFhalos/subcatalogs/"):
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
                zmin=k*zwidth
                zmax=(k+1)*zwidth
                wzmin=z>zmin
                wzmax=z<zmax
                where=wxmin*wxmax*wymin*wymax*wzmin*wzmax
                n_array.append( len(x[where]) )
                #print first
                if first==1:index=[np.array(range(len(x)))[where]];first=0
                else:index=index+[np.array(range(len(x)))[where]]
    print "finishing sub-catalog creation"
    return index, n_array
                
if indexes_subcat:
    index,n_array=subcatalogs(x,y,z,mass,savepath=cat_path)

cat_path="../FOFhalos/subcatalogs/"

# This is to plot the distrubion of the number of halos in the subcatalog in different mass bins
# plot_ndist is set to 0 by default

if plot_ndist:
    
    mmin=np.arange(9.5,12.0,0.5)
    mmax=np.arange(10.0,13.0,0.5)

    figure()
    #n, bins, patches = plt.hist(np.log10(n_array)-np.log10(np.percentile(n_array,50)), np.sqrt(len(n_array))/2.0, normed=False,histtype='step', cumulative=False,label='all')
    file=open("../output/mminmax_density"+sufix+".txt","w")
    file.write("")
    file.close()
    file=open("../output/mminmax_density"+sufix+".txt","a")
    file1=open("../output/mmedian"+sufix+".txt","w")
    file1.write("")
    file1.close()
    file1=open("../output/mmedian"+sufix+".txt","a")
    mmedian=[]
    low=[]
    up=[]
    low2=[]
    up2=[]
    MMIN=[]
    MMAX=[]

    
    for Mmin,Mmax in itertools.product(mmin,mmax):
        if Mmin>=Mmax-0.00001: continue
        Nmass=[]
        #print Mmin,Mmax
        count=0
        wmin=log10(mass)>Mmin; wmax=log10(mass)<Mmax; wm=wmin*wmax; mm=log10(mass[wm]);mmedian=np.append(mmedian,percentile(mm,50))
        up=np.append(up,percentile(mm,84))
        low=np.append(low,percentile(mm,16))

        up2=np.append(up2,percentile(mm,97.5))
        low2=np.append(low2,percentile(mm,2.5))

        MMIN=np.append(MMIN,Mmin)
        MMAX=np.append(MMAX,Mmax)
        strmmed="%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n"%(Mmin,Mmax,percentile(mm,50),percentile(mm,16),percentile(mm,84),percentile(mm,2.5),percentile(mm,97.5))
        file1.write(strmmed)
        for i,j,k in itertools.product(range(nx),range(ny),range(nz)):
            filename=cat_path+"catalog_"+str(i)+"_"+str(j)+"_"+str(k)+".txt"
            #x,y,z,m=np.genfromtxt(filename,unpack=1)
            
            x1=x[index[count]];y1=y[index[count]];z1=z[index[count]];m1=mass[index[count]]
            
        

            m1=np.log10(m1)
            wmin=m1>Mmin; wmax=m1<Mmax
            wm=wmin*wmax
            xm=x1[wm]; ym=y1[wm]; zm=z1[wm]; mm=m1[wm]
            Nmass.append(len(xm))
            del(xm);del(ym);del(zm);del(mm)
            count=count+1
            #print i,j,k
        n, bins, patches = plt.hist(np.log10(Nmass)-np.log10(np.percentile(Nmass,50)), np.sqrt(len(Nmass))/2.0, normed=False,histtype='step', cumulative=False,label=r"$M_{\rm min}=$"+'%.1f'%Mmin+" "+r"$M_{\rm max}=$"+'%.1f'%Mmax)
        #Nmass=np.append(np.array([mmin,mmax]), Nmass)
        stringNmass="%10.2f %10.2f "%(Mmin,Mmax)
        
        for ns in Nmass: stringNmass+="% 10.2f"%ns
        stringNmass+="\n"
        file.write(stringNmass)
        #plt.legend()
    xlabel(r"$\log\left(N_{\rm halos}/\tilde{N}_{\rm halos}\right)$")
    ylabel(r"$N_{\rm fields}$")
    xlim(-0.25,0.50)
    legend(fontsize=9)
    savefig('../output/ndist.png',format='png', dpi=150,bbox_inches='tight')
    file.close()
    file1.close()
    del(file)
    del(file1)

acf_file='../obs/ACF/Bielby2015.txt'
theta,corrobs,dcorrobs=np.genfromtxt(acf_file,unpack=1)
theta=theta*60.0
theta=np.log10(theta)
    
theta_bins=theta-(theta[1]-theta[0])/2.0
theta_bins=np.append(theta_bins,theta[-1]+(theta[1]-theta[0])/2.0)
th_min=theta_bins[0];th_max=theta_bins[-1]
if compute_RR:
    RR=np.zeros_like(theta)
    for i in range(ncatRR): xr= xwidth*(np.random.random_sample(Nlae_mean)); yr= ywidth*(np.random.random_sample(Nlae_mean));RR1,bins=RR_histogram(xr,yr,Dc,th_min,th_max,theta_bins,cat_number=1,logtheta=1);RR=RR+RR1
    
    RR=RR/(1.0*ncatRR)
    np.savetxt('../output/RR.txt',np.transpose([theta,RR]))
else:
    theta,RR=np.genfromtxt('../output/RR.txt',unpack=1)
def corelationlike(Mmin,dM,focc,x,y,z,mass,index,n_array,xwidth,ywidth,zwidth,Dc=Dc,Da=Da,Dl=Dl,cat_number=1,logtheta=1,convert_to_arcsec=1,estimator='landy',acf_file='../obs/ACF/Bielby2015.txt',Nlae=Nlae_mean,RR1=RR,return_corr=0):
    """ This function computes the likelyhood of a model with Mmin, Mmax, focc to macth the observational ACF and mean number density of galaxies.
    It is fundamental for the EMCEE code.""" 
    Mmax=Mmin+dM
    print "Mmin=",Mmin," Mmax=",Mmax," focc=",focc
    theta,corrobs,dcorrobs=np.genfromtxt(acf_file,unpack=1)
    if convert_to_arcsec: theta=theta*60.0
    if logtheta: theta=np.log10(theta)
    
    theta_bins=theta-(theta[1]-theta[0])/2.0
    theta_bins=np.append(theta_bins,theta[-1]+(theta[1]-theta[0])/2.0)
    th_min=theta_bins[0];th_max=theta_bins[-1]
    corr=np.zeros( len(theta)  )
    m=np.log10(mass)
    if Mmin>=Mmax-0.00001 or focc>1 or focc<0 or Mmin<9.19 or Mmax<9.24 or  Mmax>13.4 or Mmin>12.5: chi2=np.inf; print chi2; return -1.0*np.inf
    else:
        
        #print th_max, th_min, "thmax-min"
        count=1
        #for i,j,k in itertools.product(range(nx),range(ny),range(nz)):
        
        x1=x[index[count]];y1=y[index[count]];m1=m[index[count]]
        wmin=m1>Mmin; wmax=m1<Mmax
        wm=wmin*wmax
        #print x1, y1, wm
        #raw_input("Press Enter to continue...")
        xm=x1[wm]; ym=y1[wm]; n_points1=len(xm); n_points=len(xm)
        Nhalo=n_points*mean(n_array)/n_array[0] # a crude stimation of the mean number of halos with Mmin<m<Mmax in a box
        focch=Nlae/Nhalo
        Nlae_teo=focc*Nhalo
        dNlae=np.sqrt(Nlae)
        
        if focc==1:
            Nlae_teo=Nlae
            if (Nhalo<Nlae/3):
                chi2=np.inf
                print chi2
                return -1*chi2
        if  (1.0*Nlae/Nlae_teo)>6.3 or (1.0*Nlae_teo/Nlae)>6.3:
            chi2=np.inf
            print chi2
            return -1*chi2
        else:
            count=0
            #print count
            corr=np.zeros( len(theta)  )
            corr_arr=np.zeros( (len(theta),nx*ny*nz)  )
            narray=[]
            for i,j,k in itertools.product(range(nx),range(ny),range(nz)):
                
                x1=x[index[count]];y1=y[index[count]];m1=m[index[count]]
                wmin=m1>Mmin; wmax=m1<Mmax
                wm=wmin*wmax
                xm=x1[wm]; ym=y1[wm]; n_points1=len(xm);n_points=len(xm)
                

                if Nlae<=n_points:
                    
                    rpos=randsample(n_points1,Nlae)
                    xm1=xm[rpos];ym1=ym[rpos]
                    n_points=Nlae
                    DD,bins=DD_histogram(xm1,ym1,Dc,th_min,th_max,theta_bins,logtheta=logtheta)
                    
                     
                    xr= xwidth*(i+np.random.random_sample(n_points*cat_number))
                    yr= ywidth*(j+np.random.random_sample(n_points*cat_number))
                    if len(RR1)!=0:
                        RR=RR1
                    else:
                        RR,bins=RR_histogram(xr,yr,Dc,th_min,th_max,theta_bins,cat_number=cat_number,logtheta=logtheta)
                    DR,bins=DR_histogram(xm1,ym1,xr,yr,Dc,th_min,th_max,theta_bins,cat_number=cat_number,logtheta=logtheta)

                    
                   
                else:
                    DD,bins=DD_histogram(xm,ym,Dc,th_min,th_max,theta_bins,logtheta=logtheta)
                    xr= xwidth*(i+np.random.random_sample(n_points*cat_number))
                    yr= ywidth*(j+np.random.random_sample(n_points*cat_number))
                    if RR1!=0:
                        RR=RR1
                    else:
                        RR,bins=RR_histogram(xr,yr,Dc,th_min,th_max,theta_bins,cat_number=cat_number,logtheta=logtheta)
                        
                    DR,bins=DR_histogram(xm,ym,xr,yr,Dc,th_min,th_max,theta_bins,cat_number=cat_number,logtheta=logtheta)
                    
                DD[DD<=0]=1.0/(1.0*Nlae*Nlae)
                RR[RR<=0]=1.0/(1.0*Nlae*Nlae)
                DR[DR<=0]=1.0/(1.0*Nlae*Nlae)
                    
                narray=np.append(narray,n_points1)
                #print "DD=",DD," RR=",RR," DR=",DR
                if estimator=='landy':corr1=landy_correlation(DD,RR,DR);corr=corr+corr1;corr_arr[:,count]=corr1
                elif estimator=='peebles':corr1=peebles_correlation(DD,DR);corr=corr+corr1;corr_arr[:,count]=corr1
                elif estimator=='standard':corr1=standard_correlation(DD,RR);corr=corr+corr1;corr_arr[:,count]=corr1
                else:corr1=landy_correlation(DD,RR,DR);corr=corr+corr1;corr_arr[:,count]=np.ma.masked_invalid(corr1)
                    
                count=count+1
                #print count, " end"
            corr=np.ma.masked_invalid(corr/count)
            
            std=np.array([np.std(corr_arr[i,:]) for i in range(len(theta))])
            low=np.array([np.percentile(corr_arr[i,:],16) for i in range(len(theta))])
            up=np.array([np.percentile(corr_arr[i,:],84) for i in range(len(theta))])
            median=np.array([np.percentile(corr_arr[i,:],50) for i in range(len(theta))])
            err_l=median-low
            err_u=up-median
            #print "std,err_l,err_u",std,err_l,err_u
            Nhalo=np.mean(n_array)
            focch=Nlae/Nhalo
            Nlae_teo=focc*Nhalo
            dNlae=np.sqrt(Nlae)
            dcorr=dcorrobs+0.000001
            for i in range(len(theta)):
                if (corr-corrobs)[i]<0:
                    dcorr[i]=dcorrobs[i]+err_u[i]
                else:
                    dcorr[i]=dcorrobs[i]+err_l[i]
        if focc==1:
            chi2=np.sum( (corr-corrobs)*(corr-corrobs)/(dcorr*dcorr) )/(np.ma.count(corr)-1)
        else:
            dcorr=dcorrobs
            chi2=((Nlae_teo-Nlae)*(Nlae_teo-Nlae)/(dNlae*dNlae))+np.sum( (corr-corrobs)*(corr-corrobs)/(dcorr*dcorr) )/(np.ma.count(corr)-1) 
        print "chi2",chi2
    #print "corr=",corr,(np.ma.count(corr))
    #print "corrobs=",corrobs
    #raw_input("Press Enter to continue...")
    if return_corr==1:
        return theta,corr,err_u,err_l
    else:
        return -1*chi2


def lnlike(params,x,y,z,mass,index,n_array,xwidth,ywidth,zwidth):
    """This is the way the likelihood function should be redefined to be  passed to emcee"""
    
    Mmin,dM,focc=params
    return corelationlike(Mmin,dM,focc,x,y,z,mass,index,n_array,xwidth,ywidth,zwidth)


def lnlike1(params,x,y,z,mass,index,n_array,xwidth,ywidth,zwidth):
    Mmin,dM=params
    return corelationlike(Mmin,dM,1,x,y,z,mass,index,n_array,xwidth,ywidth,zwidth)

if run_emcee:
    #Setting emcee 6
    if ndim==3:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike, args=(x,y,z,mass,index,n_array,xwidth,ywidth,zwidth),threads=ncores)
    
    if ndim==2:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike1, args=(x,y,z,mass,index,n_array,xwidth,ywidth,zwidth),threads=ncores)
    # Runing emcee for a few times 
    to=time.time()
    pos0, prob, state = sampler.run_mcmc(pos, 1)
    t1=time.time()
    print "time for 1 iter=", (t1-to)/(60*60)
    print "estimated time =", niter*(t1-to)/(60*60)
    sampler.reset()

    fn = "../output/LAEBayes_mcmc_dM"+sufix+".out"


    f = open(fn, "w")
    f.close()


    fn1 = "../output/LAEBayes_mcmc_dM_chi"+sufix+".out"

    f1 = open(fn, "w")
    f1.close()


    # Restarting  emcee over a large loop of niter times
    for pos, prob, rstate in sampler.sample(pos0, prob, state, iterations=niter):
        # Write the current position to a file, one line per walker                                                                                                                     
        f = open(fn, "a")
        f1 = open(fn1, "a")
        f.write("\n".join(["\t".join([str(q) for q in p]) for p in pos]))
        #f1.write("\n".join(["\t".join([str(q) for q in p]) for p in pos]))
        f1.write("\n".join([str(q) for q in prob]))
        #print "rstate=  ", rstate
        #f1.write("\t".join([str(q) for q in prob]))
        f.write("\n")
        f1.write("\n")
        f.close()
        f1.close()





