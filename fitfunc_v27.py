import sys
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy
import scipy.constants as cn
import batman

from outputs_v27 import *

def fitfunc(params, tarr, farr,cxarr,cyarr,warr,aarr,sarr,barr,carr, nphot, indlist, filters, nfilt, filnames,nddf,inmcmc, paraCNM, baseLSQ, bvars, cont, names, earr,divwhite,dwCNMarr,dwCNMind,Rs_in,sRs_lo,sRs_hi,Ms_in,sMs_lo,sMs_hi):
    zarr=np.array([])
    ypos=np.array([])

    Rsolar = 6.957e8   # IAU value
    Rjup = 7.1492e7    # IAU value
    Rearth = 6.3781e6  # IAU value
    
    GMsolar = 1.3271244e20 # IAU value
    GMjup = 1.2668653e17   # IAU value
    GMearth = 3.986004e14  # IAU value
    
    Msolar = GMsolar / cn.G # IAU suggestion
    Mjup = GMjup / cn.G     # IAU suggestion
    Mearth = GMearth / cn.G # IAU suggestion
    
    rhoSolar = 3. / (4. * np.pi) * Msolar / Rsolar**3 
    rhoJup = 3. / (4. * np.pi) * Mjup / Rjup**3 
    rhoEarth = 3. / (4. * np.pi) * Mearth / Rearth**3 

    Rs = get_PDF_Gauss(Rs_in,sRs_lo,sRs_hi,[1])[0]
    Ms = get_PDF_Gauss(Ms_in,sMs_lo,sMs_hi,[1])[0]

    vcont = np.copy(cont[:,0])

    nantest=np.isnan(np.min(params))
   
    ecc = params[5]**2+params[6]**2
    if (ecc >= 0.99):
        ecc = 0.99
        if (params[6]/np.sqrt(ecc) < 1.):
            ome2 = np.arccos(params[6]/np.sqrt(ecc))
        else:
            ome2 = 0   
        params[5] = np.sqrt(ecc)*np.sin(ome2)
    
    if (ecc>0.00001):
        if (np.abs(params[5]<0.00000001)):
            ome = np.arctan(np.abs(params[5]/params[6]))
        else:
            ome = np.abs(np.arcsin(params[5]/np.sqrt(ecc)))

        if (params[5]<0):
            if (params[6]<0):
                ome = ome + np.pi
            else:
                ome = 2.*np.pi - ome
        else:
            if (params[6]<0):
                ome = np.pi - ome            

    else:
        ome=0.
        ecc=0.
    
    efac1=np.sqrt(1.-ecc**2)/(1.+ecc*np.sin(ome))
    efac2=params[2]*(1.-ecc**2)/(1.+ecc*np.sin(ome))
    
    aRs = np.copy(params[3])   
    inc = np.copy(params[2])   
           
    # calculate the true -> eccentric -> mean anomaly at transit -> perihelion time
    TA_tra = np.pi/2. - ome
    TA_tra = np.mod(TA_tra,2.*np.pi)
    EA_tra = 2.*np.arctan( np.tan(TA_tra/2.) * np.sqrt((1.-ecc)/(1.+ecc)) )
    EA_tra = np.mod(EA_tra,2.*np.pi)
    MA_tra = EA_tra - ecc * np.sin(EA_tra)
    MA_tra = np.mod(MA_tra,2.*np.pi)
    mmotio = 2.*np.pi/params[4]   # the mean motion, i.e. angular velocity [rad/day] if we had a circular orbit
    T_peri = params[0] - MA_tra/mmotio
    
    # =========== Transit model calculation =======================
    # now, for all lightcurves, calculate the z values for the timestamps
    for j in range(nphot):
        tt = np.copy(tarr[indlist[j][0]]) # time values of lightcurve j
        MA_lc = (tt-T_peri)*mmotio
        MA_lc = np.mod(MA_lc,2*np.pi)
        # source of the below equation: http://alpheratz.net/Maple/KeplerSolve/KeplerSolve.pdf
        EA_lc = MA_lc + np.sin(MA_lc)*ecc + 1./2.*np.sin(2.*MA_lc)*ecc**2 + \
            (3./8.*np.sin(3.*MA_lc) - 1./8.*np.sin(MA_lc))*ecc**3 + \
                (1./3.*np.sin(4.*MA_lc) - 1./6.*np.sin(2*MA_lc))*ecc**4 + \
                    (1./192*np.sin(MA_lc)-27./128.*np.sin(3.*MA_lc)+125./384.*np.sin(5*MA_lc))*ecc**5 + \
                        (1./48.*np.sin(2.*MA_lc)+27./80.*np.sin(6.*MA_lc)-4./15.*np.sin(4.*MA_lc))*ecc**6
        EA_lc = np.mod(EA_lc,2*np.pi)
        TA_lc = 2.*np.arctan(np.tan(EA_lc/2.) * np.sqrt((1.+ecc)/(1.-ecc)) )
        TA_lc = np.mod(TA_lc,2*np.pi)
        
        # calculate the z array using the same formalism as Batman (Kreidberg 2015)
        z_ma = aRs * (1.0 - ecc**2)/(1.0 + ecc*np.cos(TA_lc)) * np.sqrt(1.0 - np.sin(ome + TA_lc)**2 * np.sin(inc)**2)
        
       # calculate the z array using the formalism of M. Gillon 
        R_lc = aRs*(1.-ecc*np.cos(EA_lc))  #normalized (to Rs) planet-star separation
        b_lc = params[2]*(1.-ecc*np.cos(EA_lc))
        y_lc = np.sqrt(R_lc**2 - b_lc**2)*np.cos(TA_lc + ome - np.pi/2.)

        
        zarr = np.concatenate((zarr,z_ma), axis=0)
        ypos = np.concatenate((ypos,y_lc), axis=0)

    cnmarr=[]  
    
    # calculate the mid-transit points for each lightcurve
    marr=np.array([])
    for j in range(nphot):
        # number of periods elapsed between given T0 and lightcurve start
        tt = np.copy(tarr[indlist[j][0]]) # time values of lightcurve j
        ft = np.copy(farr[indlist[j][0]]) # flux values of lightcurve j
        et = np.copy(earr[indlist[j][0]]) # error values of lightcurve j
        z  = np.copy(zarr[indlist[j][0]]) # z values of lightcurve j
        y  = np.copy(ypos[indlist[j][0]]) # y values of lightcurve j
        
        #  the supplementary parameter: am, cx, cy, fwhm, sky
        am=np.copy(aarr[indlist[j][0]])
        cx=np.copy(cxarr[indlist[j][0]])
        cy=np.copy(cyarr[indlist[j][0]])
        fwhm=np.copy(warr[indlist[j][0]])
        sky=np.copy(sarr[indlist[j][0]])
        npo=len(z)                # number of lc points
        
        # normalize the timestamps to the center of the transit
        delta = np.round(np.divide((tt[0]-params[0]),params[4]))
        T0_lc=params[0]+delta*params[4]  
        ts=tt-T0_lc
        
        # identify the filter index of this LC
        k = np.where(filnames == filters[j])  # k is the index of the LC in the filnames array
        k = np.asscalar(k[0])
        u1ind = 8+nddf+4*k  # index in params of the first LD coeff of this filter
        u2ind = 9+nddf+4*k  # index in params of the second LD coeff of this filter
       
        RR=np.copy(params[1])
        
        # convert the LD coefficients to u1 and u2
        u1=(params[u1ind] + params[u2ind])/3.
        u2=(params[u1ind] - 2.*params[u2ind])/3.
                
        batpars = batman.TransitParams()
        batpars.t0 = np.copy(params[0])                       #time of inferior conjunction
        batpars.per = np.copy(params[4])         #orbital period
        batpars.rp = np.copy(params[1])          #planet radius (in units of stellar radii)
        batpars.a =  np.copy(params[3])          #semi-major axis (in units of stellar radii)
        batpars.inc = np.copy(params[2]*180/np.pi)         #orbital inclination (in degrees)
        batpars.ecc = np.copy(ecc)               #eccentricity
        batpars.w = np.copy(ome*180/np.pi)                 #longitude of periastron (in degrees)
        batpars.u = [np.copy(u1), np.copy(u2)]   #limb darkening coefficients [u1, u2]
        batpars.limb_dark = "quadratic"          #limb darkening model
        t = np.copy(tt)
        m = batman.TransitModel(batpars, t)      #initializes model
        mm0 = m.light_curve(batpars)             #calculates light curve        
        
        mm0[y<0]=1.
        # correct here for the contamination
        mm = (mm0-1)/(vcont[k]+1) + 1.
 
        bfstart= 8+nddf+nfilt*4   # the first index in the param array that refers to a baseline function
        incoeff = list(range(bfstart+j*19,bfstart+j*19+19))  # the indices for the coefficients for the base function        
        
        if (baseLSQ == 'y'):
            mres=ft/mm
            ivars = np.copy(bvars[j][0])
            incoeff=np.array(incoeff)
            coeffstart = np.copy(params[incoeff[ivars]])   # RANDOM NOTE: you can use lists as indices to np.arrays but not np.arrays to lists or lists to lists
            icoeff,dump = scipy.optimize.leastsq(para_minfunc, coeffstart, args=(ivars, mm, ft, ts, am, cx, cy, fwhm, sky))
            coeff = np.copy(params[incoeff])   # the full coefficients -- set to those defined in params (in case any are fixed non-zero)
            coeff[ivars] = np.copy(icoeff)     # and the variable ones are set to the result from the minimization
            
        else:        
            coeff = np.copy(params[incoeff])   # the coefficients for the base function
    
        bfunc=basefunc_noCNM(coeff, ts, am, cx, cy, fwhm, sky)
        mod=mm*bfunc

        fco=ft/bfunc
        
        if (inmcmc == 'n'):
            outfile=names[j][:-4]+'_out.dat'
            of=open(outfile,'w')
        
            for k in range(len(tt)):
                of.write('%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n' % (tt[k], ft[k], et[k], mod[k],bfunc[k],mm[k],fco[k])) 
        
            of.close()
        
        if (min(mm0)>0.99999):
            mod = mod + 10000      
        
        marr=np.append(marr,mod)
    
    chisq = np.sum((marr-farr)**2/earr**2)    
    
    return marr

def basefunc_noCNM(coeff, ts, am, cx, cy, fwhm, sky):
    # the full baseline function calculated with the coefficients given; of which some are not jumping and set to 0
    bfunc=coeff[0]+coeff[1]*ts+coeff[2]*np.power(ts,2)+coeff[3]*np.power(ts,3)+coeff[4]*np.power(ts,4)+coeff[5]*am+coeff[6]*np.power(am,2)+coeff[7]*cx+coeff[8]*np.power(cx,2)+coeff[9]*cy+coeff[10]*np.power(cy,2)+coeff[11]*cy*cx+coeff[12]*fwhm+coeff[13]*np.power(fwhm,2)+coeff[14]*sky+coeff[15]*np.power(sky,2)+coeff[16]*np.sin(ts*coeff[17]+coeff[18])
    
    return bfunc

def para_minfunc(icoeff, ivars, mm, ft, ts, am, cx, cy, fwhm, sky):
    icoeff_full = np.zeros(19)
    icoeff_full[ivars] = np.copy(icoeff)
    return (ft - mm * basefunc_noCNM(icoeff_full, ts, am, cx, cy, fwhm, sky))

