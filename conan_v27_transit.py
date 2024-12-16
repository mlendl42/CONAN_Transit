import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#import MCcubed as mc3   ---> https://github.com/pcubillos/mc3
import mc3

from basecoeff_v27 import *
import fitfunc_v27 as ff
from plots_v27 import *
from ecc_om_par import *
from outputs_v27 import *


file=open('input.dat')
dump=file.readline()
dump=file.readline()
fpath=file.readline()         # the path where the files are
fpath=fpath.rstrip()
dump=file.readline()
dump=file.readline()
# now read in the lightcurves one after the other
names=[]                    # array where the LC filenames are supposed to go
filters=[]                  # array where the filter names are supposed to go
lamdas=[]
bases=[]                    # array where the baseline exponents are supposed to go
groups=[]                   # array where the group indices are supposed to go
grbases=[]
dump=file.readline()        # read the next line
dump=file.readline()        # read the next line
while dump[0] != '#':           # if it is not starting with # then
    adump=dump.split()          # split it
    names.append(adump[0])      # append the first field to the name array
    filters.append(adump[1])    # append the second field to the filters array
    lamdas.append(adump[2])     # append the second field to the filters array
    strbase=adump[3:10]         # string array of the baseline function exponents
    base = [int(i) for i in strbase]
    bases.append(base)
    group = int(adump[10])
    groups.append(group)
    grbase=int(adump[9])
    grbases.append(grbase)
    dump=file.readline()    # and read the following line - either a lc or something starting with "#"

nphot=len(names)             # the number of photometry input files
njumpphot=np.zeros(nphot)
filnames=np.array(list(sorted(set(filters),key=filters.index))) 
ulamdas=np.array(list(sorted(set(lamdas),key=lamdas.index))) 
grnames=np.array(list(sorted(set(groups))))
nfilt=len(filnames)
ngroup=len(grnames)
dump=file.readline()

RVbases=[]
gammas=[]
gamsteps=[]
gampri=[]
gamprilo=[]
gamprihi=[]

dump=file.readline()
dump=file.readline()
adump=dump.split()
adump[3] = (0. if adump[1] == 'n' else adump[3])
adump[7] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[7])
adump[8] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[8])
adump[9] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[9])
rprs_in=[float(adump[2]),float(adump[3]),float(adump[4]),float(adump[5]),float(adump[7]),float(adump[8]),float(adump[9])] 
rprs0=np.copy(rprs_in[0])
if rprs_in[1] != 0.:
    njumpphot=njumpphot+1

erprs0=np.copy(0)    
dump=file.readline()
adump=dump.split()
adump[3] = (0. if adump[1] == 'n' else adump[3])
adump[7] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[7])
adump[8] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[8])
adump[9] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[9])
inc_in=[float(adump[2]),float(adump[3]),float(adump[4]),float(adump[5]),float(adump[7]),float(adump[8]),float(adump[9])]

if inc_in[1] != 0.:
    njumpphot=njumpphot+1
    
dump=file.readline()
adump=dump.split()
adump[3] = (0. if adump[1] == 'n' else adump[3])
adump[7] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[7])
adump[8] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[8])
adump[9] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[9])
aRs_in=[float(adump[2]),float(adump[3]),float(adump[4]),float(adump[5]),float(adump[7]),float(adump[8]),float(adump[9])]
if aRs_in[1] != 0.:
    njumpphot=njumpphot+1

dump=file.readline()
adump=dump.split()
adump[3] = (0. if adump[1] == 'n' else adump[3])
adump[7] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[7])
adump[8] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[8])
adump[9] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[9])
T0_in=[float(adump[2]),float(adump[3]),float(adump[4]),float(adump[5]),float(adump[7]),float(adump[8]),float(adump[9])]
if T0_in[1] != 0.:
    #njumpRV=njumpRV+1
    njumpphot=njumpphot+1
    
dump=file.readline()
adump=dump.split()
adump[3] = (0. if adump[1] == 'n' else adump[3])
adump[7] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[7])
adump[8] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[8])
adump[9] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[9])
per_in=[float(adump[2]),float(adump[3]),float(adump[4]),float(adump[5]),float(adump[7]),float(adump[8]),float(adump[9])]
if per_in[1] != 0.:
    njumpRV=njumpRV+1
    njumpphot=njumpphot+1
    
dump=file.readline()
adump=dump.split()
adump[3] = (0. if adump[1] == 'n' else adump[3])
adump[7] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[7])
adump[8] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[8])
adump[9] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[9])
eccpri=np.copy(adump[6])
ecc_in=[float(adump[2]),float(adump[3]),float(adump[4]),float(adump[5]),float(adump[7]),float(adump[8]),float(adump[9])]
if ecc_in[1] != 0.:
    njumpRV=njumpRV+1
    
dump=file.readline()
adump=dump.split()
adump[3] = (0. if adump[1] == 'n' else adump[3])
adump[7] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[7])
adump[8] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[8])
adump[9] = (0. if (adump[6] == 'n' or adump[1] == 'n' or float(adump[3]) == 0.) else adump[9])
opri=np.copy(adump[6])
omega_in=[float(adump[2]),float(adump[3]),float(adump[4]),float(adump[5]),float(adump[7]),float(adump[8]),float(adump[9])]
omega_in=np.multiply(omega_in,np.pi)/180.
if omega_in[1] != 0.:
    njumpRV=njumpRV+1
    
K_in=[0.,0.,0.,0.,0.,0.,0.]

if ((eccpri == 'y' and opri == 'n') or (eccpri == 'n' and opri == 'y')):
    print('priors on eccentricity and omega: either both on or both off')
    print(nothing)
    
eos_in,eoc_in = ecc_om_par(ecc_in, omega_in)

ddfYN = 'n'
divwhite = 'n'
dump=file.readline()

grprs=np.zeros(ngroup)   # the group rprs values
egrprs=np.zeros(ngroup)  # the uncertainties of the group rprs values
dwfiles = []             # the filenames of the white residuals

dwCNMarr=np.array([])      # initializing array with all the dwCNM values
dwCNMind=[]                # initializing array with the indices of each group's dwCNM values
dwind=np.array([])    

if (ddfYN=='n' and np.max(grbases)>0):
    print('no dDFs but groups? Not a good idea!')
    print(base)
    print(nothing)

dump=file.readline()

# setup the arrays for the LD coefficients
c1_in=np.zeros((nfilt,7))
c2_in=np.zeros((nfilt,7))
c3_in=np.zeros((nfilt,7))
c4_in=np.zeros((nfilt,7))

for i in range(nfilt):
    dump=file.readline()
    adump=dump.split()
    j=np.where(filnames == adump[0])               # make sure the sequence in this array is the same as in the "filnames" array
    k=np.where(np.array(filters) == adump[0])

    c1_in[j,:]=[float(adump[2]),float(adump[3]),-3.,3.,float(adump[2]),float(adump[4]),float(adump[5])]  # the limits are -3 and 3 => very safe
    c1_in[j,5] = (0. if (adump[1] == 'n' or float(adump[3]) == 0.) else c1_in[j,5])
    c1_in[j,6] = (0. if (adump[1] == 'n' or float(adump[3]) == 0.) else c1_in[j,6])
    if c1_in[j,1] != 0.:
        njumpphot[k]=njumpphot[k]+1
        
    c2_in[j,:]=[float(adump[6]),float(adump[7]),-3.,3.,float(adump[6]),float(adump[8]),float(adump[9])]
    c2_in[j,5] = (0. if (adump[1] == 'n' or float(adump[7]) == 0.) else c2_in[j,5])
    c2_in[j,6] = (0. if (adump[1] == 'n' or float(adump[7]) == 0.) else c2_in[j,6])
    if c2_in[j,1] != 0.:
        njumpphot[k]=njumpphot[k]+1

    c3_in[j,:]=[float(adump[10]),float(adump[11]),-3.,3.,float(adump[10]),float(adump[12]),float(adump[13])]
    c3_in[j,5] = (0. if (adump[1] == 'n' or float(adump[11]) == 0.) else c3_in[j,5])
    c3_in[j,6] = (0. if (adump[1] == 'n' or float(adump[11]) == 0.) else c3_in[j,6])
    if c3_in[j,1] != 0.:
        njumpphot[k]=njumpphot[k]+1

    c4_in[j,:]=[float(adump[14]),float(adump[15]),-3.,3.,float(adump[14]),float(adump[16]),float(adump[17])]
    c4_in[j,5] = (0. if (adump[1] == 'n' or float(adump[15]) == 0.) else c4_in[j,5])
    c4_in[j,6] = (0. if (adump[1] == 'n' or float(adump[15]) == 0.) else c4_in[j,6])
    if c4_in[j,1] != 0.:
        njumpphot[k]=njumpphot[k]+1

    # convert the input u1, u2 to c1, c2 if the limb-darkening law is quadratic
    if (c3_in[j,0] == 0. and c4_in[j,0]==0 and c3_in[j,1] == 0. and c4_in[j,1] == 0.):
        print('Limb-darkening law: quadratic')
        v1=2.*c1_in[j,0]+c2_in[j,0]
        v2=c1_in[j,0]-c2_in[j,0]
        ev1=np.sqrt(4.*c1_in[j,1]**2+c2_in[j,1]**2)
        ev2=np.sqrt(c1_in[j,1]**2+c2_in[j,1]**2)
        lov1=np.sqrt(4.*c1_in[j,5]**2+c2_in[j,5]**2)
        lov2=np.sqrt(c1_in[j,5]**2+c2_in[j,5]**2)
        hiv1=np.sqrt(4.*c1_in[j,6]**2+c2_in[j,6]**2)
        hiv2=np.sqrt(c1_in[j,6]**2+c2_in[j,6]**2) 
        c1_in[j,0]=np.copy(v1)
        c2_in[j,0]=np.copy(v2)
        c1_in[j,4]=np.copy(v1)
        c2_in[j,4]=np.copy(v2)
        c1_in[j,1]=np.copy(ev1)
        c2_in[j,1]=np.copy(ev2)
        if (adump[1] == 'y'):    # prior on LDs
            c1_in[j,5]=np.copy(lov1)
            c1_in[j,6]=np.copy(hiv1)
            c2_in[j,5]=np.copy(lov2)
            c2_in[j,6]=np.copy(hiv2)


# setup the arrays for the contamination factors
cont=np.zeros((nfilt,2))   # contains for each filter [value, error]

# read the stellar input properties
dump=file.readline()
dump=file.readline()
dump=file.readline()
adump=dump.split()
Rs_in = float(adump[1])
sRs_lo = float(adump[2])
sRs_hi = float(adump[3])
dump=file.readline()
adump=dump.split()
Ms_in = float(adump[1])
sMs_lo = float(adump[2])
sMs_hi = float(adump[3])
dump=file.readline()
adump=dump.split()
howstellar = adump[1]

# read the MCMC setup 
dump=file.readline()
dump=file.readline()
adump=dump.split() 
nsamples=int(adump[1])   # total number of integrations
dump=file.readline()
adump=dump.split()
nchains=int(adump[1])  #  number of chains
dump=file.readline()
adump=dump.split()
nproc=int(adump[1])   #  number of processes
dump=file.readline()
adump=dump.split()
burnin=int(adump[1])    # Length of bun-in
dump=file.readline()
adump=dump.split()
walk=adump[1]            # Differential Evolution?          
dump=file.readline()
adump=dump.split()
grtest = True if adump[1] == 'y' else False  # GRtest done?
dump=file.readline()
adump=dump.split()
plots = True if adump[1] == 'y' else False  # Make plots done
dump=file.readline()
adump=dump.split()
leastsq = True if adump[1] == 'y' else False  # Do least-square?
dump=file.readline()
adump=dump.split()
savefile = adump[1]   # Filename of save file
dump=file.readline()
adump=dump.split()
savemodel = adump[1]   # Filename of model save file
adaptBL = 'y'
paraCNM = 'n'
baseLSQ = 'y'
                                   
inc_in=np.multiply(inc_in,np.pi)/180.
                                   
tarr=np.array([]) # initializing array with all timestamps
farr=np.array([]) # initializing array with all flux values
earr=np.array([]) # initializing array with all error values
xarr=np.array([]) # initializing array with all x_shift values
yarr=np.array([]) # initializing array with all y_shift values
aarr=np.array([]) # initializing array with all airmass values
warr=np.array([]) # initializing array with all fwhm values
sarr=np.array([]) # initializing array with all sky values
lind=np.array([]) # initializing array with the lightcurve indices
barr=np.array([]) # initializing array with all bisector values
carr=np.array([]) # initializing array with all contrast values

indlist = []    # the list of the array indices
bvars    = []   # a list that will contain lists of [0, 1] for each of the baseline parameters, for each of the LCs. 0 means it's fixed. 1 means it's variable

if ddfYN == 'y':   # if ddFs are fit: set the Rp/Rs to the value specified at the jump parameters, and fix it.
    rprs_in=[rprs_in[0],0,0,1,0,0,0]
    nddf=nfilt
else:
    nddf=0

# set up the parameters
params   = np.array([T0_in[0], rprs_in[0], inc_in[0], aRs_in[0], per_in[0], eos_in[0], eoc_in[0], K_in[0]])  # initial guess params
stepsize = np.array([T0_in[1], rprs_in[1], inc_in[1], aRs_in[1], per_in[1], eos_in[1], eoc_in[1], K_in[1]])  # stepsizes
pmin     = np.array([T0_in[2], rprs_in[2], inc_in[2], aRs_in[2], per_in[2], eos_in[2], eoc_in[2], K_in[2]])  # Boundaries (min)
pmax     = np.array([T0_in[3], rprs_in[3], inc_in[3], aRs_in[3], per_in[3], eos_in[3], eoc_in[3], K_in[3]])  # Boundaries (max)
prior    = np.array([T0_in[4], rprs_in[4], inc_in[4], aRs_in[4], per_in[4], eos_in[4], eoc_in[4], K_in[4]])  # Prior centers
priorlow = np.array([T0_in[5], rprs_in[5], inc_in[5], aRs_in[5], per_in[5], eos_in[5], eoc_in[5], K_in[5]])  # Prior sigma low side
priorup  = np.array([T0_in[6], rprs_in[6], inc_in[6], aRs_in[6], per_in[6], eos_in[6], eoc_in[6], K_in[6]])  # Prior sigma high side
pnames   = np.array(['T_0', 'RpRs', 'inc_[d]', 'aRs', 'Period_[d]', 'esin(w)', 'ecos(w)', 'K']) # Parameter names

if (divwhite=='y'):           # do we do a divide-white? If yes, then fix all the transit shape parameters
    stepsize[0:6] = 0
    prior[0:6] = 0

for i in range(nfilt):  # add the LD coefficients for the filters to the parameters
    params=np.concatenate((params, [c1_in[i,0], c2_in[i,0], c3_in[i,0], c4_in[i,0]]))
    stepsize=np.concatenate((stepsize, [c1_in[i,1], c2_in[i,1], c3_in[i,1], c4_in[i,1]]))
    pmin=np.concatenate((pmin, [c1_in[i,2], c2_in[i,2], c3_in[i,2], c4_in[i,2]]))
    pmax=np.concatenate((pmax, [c1_in[i,3], c2_in[i,3], c3_in[i,3], c4_in[i,3]]))
    prior=np.concatenate((prior, [c1_in[i,4], c2_in[i,4], c3_in[i,4], c4_in[i,4]]))
    priorlow=np.concatenate((priorlow, [c1_in[i,5], c2_in[i,5], c3_in[i,5], c4_in[i,5]]))
    priorup=np.concatenate((priorup, [c1_in[i,6], c2_in[i,6], c3_in[i,6], c4_in[i,6]]))
    pnames=np.concatenate((pnames, [filnames[i]+'_c1',filnames[i]+'_c2',filnames[i]+'_c3',filnames[i]+'_c4']))
    
nbc_tot = np.copy(0)  # total number of baseline coefficients let to vary (leastsq OR jumping)

for i in range(nphot):
    t, flux, err, xshift, yshift, airm, fwhm, sky, eti = np.loadtxt(fpath+names[i], usecols=(0,1,2,3,4,5,6,7,8), unpack = True)  # reading in the data
    if (divwhite=='y'): # if the divide - white is activated, divide the lcs by the white noise model before proceeding
        dwCNM = np.copy(dwCNMarr[dwCNMind[groups[i]-1]])
        flux=np.copy(flux/dwCNM)
    
    sky=sky-np.mean(sky)
    tarr=np.concatenate((tarr,t), axis=0)
    farr=np.concatenate((farr,flux), axis=0)
    earr=np.concatenate((earr,err), axis=0)
    xarr=np.concatenate((xarr,xshift), axis=0)
    yarr=np.concatenate((yarr,yshift), axis=0)
    aarr=np.concatenate((aarr,airm), axis=0)
    warr=np.concatenate((warr,fwhm), axis=0)
    sarr=np.concatenate((sarr,sky), axis=0)
    barr=np.concatenate((barr,np.zeros(len(t),dtype=np.int)), axis=0)   # bisector array: filled with 0s
    carr=np.concatenate((carr,np.zeros(len(t),dtype=np.int)), axis=0)   # contrast array: filled with 0s
    lind=np.concatenate((lind,np.zeros(len(t),dtype=np.int)+i), axis=0)
    indices=np.where(lind==i)
    indlist.append(indices)
    
    A_in,B_in,C_in,D_in,E_in,G_in,H_in,nbc = basecoeff(bases[i])  # the baseline coefficients for this lightcurve; each is a 2D array
    nbc_tot = nbc_tot+nbc # add up the number of jumping baseline coeff
    # if the least-square fitting for the baseline is turned on (baseLSQ = 'y'), then set the stepsize of the jump parameter to 0
    if (baseLSQ == "y"):
        abvar=np.concatenate(([A_in[1,:],B_in[1,:],C_in[1,:],D_in[1,:],E_in[1,:],G_in[1,:],H_in[1,:]]))
        abind=np.where(abvar!=0.)
        bvars.append(abind)
        A_in[1,:]=B_in[1,:]=C_in[1,:]=D_in[1,:]=E_in[1,:]=G_in[1,:]=H_in[1,:]=0                             # the step sizes are set to 0 so that they are not interpreted as MCMC JUMP parameters

    # append these to the respective mcmc input arrays
    params=np.concatenate((params,A_in[0,:],B_in[0,:],C_in[0,:],D_in[0,:],E_in[0,:],G_in[0,:],H_in[0,:]))
    stepsize=np.concatenate((stepsize,A_in[1,:],B_in[1,:],C_in[1,:],D_in[1,:],E_in[1,:],G_in[1,:],H_in[1,:]))
    pmin=np.concatenate((pmin,A_in[2,:],B_in[2,:],C_in[2,:],D_in[2,:],E_in[2,:],G_in[2,:],H_in[2,:]))
    pmax=np.concatenate((pmax,A_in[3,:],B_in[3,:],C_in[3,:],D_in[3,:],E_in[3,:],G_in[3,:],H_in[3,:]))
    prior=np.concatenate((prior, np.zeros(len(A_in[0,:])+len(B_in[0,:])+len(C_in[0,:])+len(D_in[0,:])+len(E_in[0,:])+len(G_in[0,:])+len(H_in[0,:]))))
    priorlow=np.concatenate((priorlow, np.zeros(len(A_in[0,:])+len(B_in[0,:])+len(C_in[0,:])+len(D_in[0,:])+len(E_in[0,:])+len(G_in[0,:])+len(H_in[0,:]))))
    priorup=np.concatenate((priorup, np.zeros(len(A_in[0,:])+len(B_in[0,:])+len(C_in[0,:])+len(D_in[0,:])+len(E_in[0,:])+len(G_in[0,:])+len(H_in[0,:]))))
    pnames=np.concatenate((pnames, [names[i]+'_A0', names[i]+'_A1',names[i]+'_A2',names[i]+'_A3',names[i]+'_A4',names[i]+'_B1',names[i]+'_B2',names[i]+'_C1', names[i]+'_C2',names[i]+'_C3', names[i]+'_C4',names[i]+'_C5',names[i]+'_D1',names[i]+'_D2',names[i]+'_E1',names[i]+'_E2',names[i]+'_G1',names[i]+'_G2',names[i]+'_G3',names[i]+'_H1',names[i]+'_H2']))

    
inmcmc='n'
indparams = [tarr,farr,xarr,yarr,warr,aarr,sarr,barr,carr, nphot, indlist, filters, nfilt, filnames,nddf, inmcmc, paraCNM, baseLSQ, bvars, cont,names,earr,divwhite,dwCNMarr,dwCNMind,Rs_in,sRs_lo,sRs_hi,Ms_in,sMs_lo,sMs_hi]

func=ff.fitfunc

# ========= plot the initial lightcurves to see where we start=================
yval0=func(params, tarr, farr,xarr,yarr,warr,aarr,sarr,barr,carr, nphot, indlist, filters, nfilt, filnames,nddf,inmcmc, paraCNM, baseLSQ,bvars, cont,names,earr,divwhite,dwCNMarr,dwCNMind,Rs_in,sRs_lo,sRs_hi,Ms_in,sMs_lo,sMs_hi)
mcmc_plots(yval0,tarr,farr,earr,xarr,yarr,warr,aarr,sarr,barr,carr,lind, nphot, indlist, filters,names, 'init_',params)

inparams=np.copy(params)
insteps=np.copy(stepsize)

print('start least-square fit')

output_opt = mc3.fit(farr, earr, func, params, indparams=indparams,
                    pmin=pmin, pmax=pmax, pstep=stepsize, leastsq='lm')

chibp = output_opt['bestp']
chim = output_opt['best_model']# the model of the chi2 fit
    
# plot the chi2 result to check of it's looking OK
mcmc_plots(chim,tarr,farr,earr,xarr,yarr,warr,aarr,sarr,barr,carr,lind, nphot, indlist, filters,names,'lssq_', params)

rarr=farr-chim  # the full residuals

newparams=np.copy(chibp)   

inmcmc='y'
indparams = [tarr,farr,xarr,yarr,warr,aarr,sarr,barr,carr, nphot, indlist, filters, nfilt, filnames,nddf, inmcmc, paraCNM, baseLSQ, bvars, cont, names, earr,divwhite,dwCNMarr,dwCNMind,Rs_in,sRs_lo,sRs_hi,Ms_in,sMs_lo,sMs_hi]

mc3_output = mc3.sample(data=farr, uncert=earr, func=func, params=newparams, indparams=indparams,
                    pstep=stepsize, prior=prior, pmin=pmin, pmax=pmax, pnames=pnames,
                    priorlow=priorlow, priorup=priorup, chisqscale=False, leastsq=None,
                    nsamples=nsamples, nchains=nchains, sampler=walk, grtest=grtest, ncpu=nproc,
                    burnin=burnin, plots=plots, savefile=savefile, savemodel=None, log='MCMC.log')   

ijnames = np.where(stepsize != 0.)
jnames = pnames[[ijnames][0]]  # jnames are the names of the jump parameters
nijnames = np.where(stepsize == 0.)
njnames = pnames[[nijnames][0]]  # njnames are the names of the fixed parameters

bp = mc3_output['bestp']
jbp = bp[ijnames]
posterior = mc3_output['posterior'][mc3_output['zmask']]
s1bp = [mc3_output['CRlo'] ,mc3_output['CRhi']]

inmcmc='n'
yval=func(bp,      tarr, farr,xarr,yarr,warr,aarr,sarr, barr, carr,nphot, indlist, filters, nfilt, filnames,nddf,inmcmc, paraCNM, baseLSQ, bvars, cont, names, earr,divwhite,dwCNMarr,dwCNMind,Rs_in,sRs_lo,sRs_hi,Ms_in,sMs_lo,sMs_hi)

mcmc_plots(yval,tarr,farr,earr,xarr,yarr,warr,aarr,sarr,barr,carr,lind, nphot, indlist, filters, names, 'fit_', bp)

dim=posterior.shape
# calculate PDFs of the stellar parameters given 
Rs_PDF = get_PDF_Gauss(Rs_in,sRs_lo,sRs_hi,dim)
Ms_PDF = get_PDF_Gauss(Ms_in,sMs_lo,sMs_hi,dim)

medvals,maxvals=mcmc_outputs(posterior,jnames, ijnames, njnames, nijnames, bp, s1bp, ulamdas, Rs_in, Ms_in, Rs_PDF, Ms_PDF, nfilt, filnames, howstellar)

npar=len(jnames)
if (baseLSQ == "y"):
    npar = npar + nbc_tot   # add the baseline coefficients if they are done by leastsq

ndat = len(tarr)

bic=get_BIC(npar,ndat)

medp=inparams
medp[[ijnames][0]]=medvals
maxp=inparams
maxp[[ijnames][0]]=maxvals

inmcmc='n'
mval=func(medp, tarr, farr,xarr,yarr,warr,aarr,sarr,barr, carr, nphot, indlist, filters, nfilt, filnames,nddf,inmcmc, paraCNM, baseLSQ, bvars, cont, names, earr,divwhite,dwCNMarr,dwCNMind,Rs_in,sRs_lo,sRs_hi,Ms_in,sMs_lo,sMs_hi)

mcmc_plots(mval,tarr,farr,earr,xarr,yarr,warr,aarr,sarr,barr,carr,lind, nphot, indlist, filters, names, 'med_',medp)
   
mval2=func(maxp, tarr, farr,xarr,yarr,warr,aarr,sarr,barr, carr, nphot, indlist, filters, nfilt, filnames,nddf,inmcmc, paraCNM, baseLSQ, bvars, cont, names, earr,divwhite,dwCNMarr,dwCNMind,Rs_in,sRs_lo,sRs_hi,Ms_in,sMs_lo,sMs_hi)

mcmc_plots(mval2,tarr,farr,earr,xarr,yarr,warr,aarr,sarr,barr,carr,lind, nphot, indlist, filters, names, 'max_', maxp)


