def mcmc_plots(yval,tarr,farr,earr,xarr,yarr,warr,aarr,sarr,barr,carr, lind, nphot, indlist, filters,names,prefix,params):
    
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt    
    import numpy as np
    from matplotlib.ticker import FormatStrFormatter    
    from scipy.stats import binned_statistic

    for j in range(nphot):
        mod=yval[indlist[j][0]]
        time=tarr[indlist[j][0]]
        flux=farr[indlist[j][0]]
        err=earr[indlist[j][0]]
        am=aarr[indlist[j][0]]
        cx=xarr[indlist[j][0]]
        cy=yarr[indlist[j][0]]
        fwhm=warr[indlist[j][0]]
        sky=sarr[indlist[j][0]]
        
        infile=names[j][:-4]+'_out.dat'
        tt, ft, et, mt, bfunc, mm, fco = np.loadtxt(infile, usecols=(0,1,2,3,4,5,6), unpack = True)  # reading in the lightcurve data
        
        # bin the lightcurve data
        binsize=30./(24.*60.)   # 60 minute bin size in days
        bintime = binsize *24.*60.
        nbin = int((np.max(tt)-np.min(tt))/binsize)  # number of bins
        binlims=np.zeros(nbin+1)
        tbin=np.zeros(nbin)
        binnps=np.zeros(nbin)  #number of points per bin
        binlims[0]=min(tt)
        binind=[]
        for k in range(1,nbin+1):
            binlims[k]=binlims[k-1]+binsize
            for k in range(nbin):
                tbin[k]=binlims[k]+0.5*binsize
                binnps[k]=len(tt[(tt>binlims[k]) & (tt<binlims[k+1])])
        
        ftbin, dump, dump2 = binned_statistic(tt,ft,statistic='mean',bins=binlims)
        mtbin, dump, dump2 = binned_statistic(tt,mt,statistic='mean',bins=binlims)
        fcobin, dump, dump2 = binned_statistic(tt,fco,statistic='mean',bins=binlims)
        mmbin, dump, dump2 = binned_statistic(tt,mm,statistic='mean',bins=binlims)
        etbin, dump, dump2 = binned_statistic(tt,et,statistic='mean',bins=binlims)        
        etbin = etbin/np.sqrt(binnps)
        
        rmsfit=np.std(ftbin-mtbin)
        tit="%.0f"% bintime+ ' min RMS to fit '+"%.6f" % rmsfit 
        outname=prefix+names[j][:-4]+'_fit.png'
        plt.figure(10)
        plt.clf()
        plt.errorbar(tt-7000, ft, yerr=et, fmt=".g", label='data')
        plt.plot(tt-7000, mt, "-r", label='MCMC best fit')
        plt.errorbar(tbin-7000, ftbin, yerr=etbin, fmt="ob", label='binned data')
        plt.title(tit)
        plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%5.2f'))
        plt.legend(loc="best")
        plt.xlabel("HJD - 2457000")
        plt.ylabel("Relative Flux")
        plt.savefig(outname)
        
        res=flux-mod
        outname=prefix+names[j][:-4]+'_res.png'
        plt.figure(10)
        plt.clf()
        plt.errorbar(tt-7000, ft-mt, yerr=et, fmt=".g", label='residuals')
        plt.plot(tt-7000, mt-mt, "-r")
        plt.errorbar(tbin-7000, ftbin-mtbin, yerr=etbin, fmt="ob", label='binned residuals')
        plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%5.2f'))
        plt.title(tit)
        plt.legend(loc="best")
        plt.xlabel("HJD - 2457000")
        plt.ylabel("Relative Flux")
        plt.savefig(outname)
        
        outname=prefix+names[j][:-4]+'_basecor.png'
        plt.figure(10)
        plt.clf()
        plt.errorbar(tt-7000, fco, yerr=et, fmt=".g", label='baseline corrected data')
        plt.plot(tt-7000, mm, "-r")
        plt.errorbar(tbin-7000, fcobin, yerr=etbin, fmt="ob")
        plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%5.2f'))
        plt.title(tit)
        plt.legend(loc="best")
        plt.xlabel("HJD - 2457000")
        plt.ylabel("Relative Flux")
        plt.savefig(outname)

        outname=prefix+names[j][:-4]+'_sum.png'
        plt.figure(figsize=(10,15))
        plt.clf()
        
        plt.errorbar(tt-7000, ft, yerr=et, fmt=".g", label='data')
        plt.plot(tt-7000, bfunc, "-y")
        plt.plot(tt-7000, mt, "-r", label='MCMC best fit')
        plt.errorbar(tbin-7000, ftbin, yerr=etbin, fmt="ob", label='binned data')
        l = 1-min(mt)

        plt.errorbar(tt-7000, fco-l, yerr=et, fmt=".g", label='baseline corrected data')
        plt.plot(tt-7000, mm-l, "-r")
        plt.errorbar(tbin-7000, fcobin-l, yerr=etbin, fmt="ob")
        m = 1.5*l-min(mm)
        plt.errorbar(tt-7000, ft-mt-m, yerr=et, fmt=".g", label='residuals')
        plt.plot(tt-7000, mt-mt-m, "-r")
        plt.errorbar(tbin-7000, ftbin-mtbin-m, yerr=etbin, fmt="ob", label='binned residuals')
       
        plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%5.2f'))
        plt.title(tit)
#        plt.legend(loc="best")
        plt.xlabel("HJD - 2457000")
        plt.ylabel("Relative Flux")
        plt.savefig(outname)


def param_hist(vals,pname,mv,s1v,s3v,mav,s1m,s3m):
    
    # this needs to be written. Just make a histogram plot of the parameter (values vals), label it (pname), and indicate the 1- and 3- sigma limits (s1, s3)

    import numpy as np
    import matplotlib 
    matplotlib.use('Agg')
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt

    fig = plt.figure()

    matplotlib.rcParams.update({'font.size': 10})
    
    num_bins = 50
    # the histogram of the data
    n, bins, patches = plt.hist(vals, num_bins, facecolor='green', alpha=0.5)
    plt.xlabel(pname)
    plt.ylabel('N samples')
    l=plt.axvline(x=mv, color='r')    
    l=plt.axvline(x=mv+s1v[0], color='b')
    l=plt.axvline(x=mv+s1v[1], color='b')
    l=plt.axvline(x=mv+s3v[0], color='y')
    l=plt.axvline(x=mv+s3v[1], color='y')

    l=plt.axvline(x=mav, color='r', linestyle='dashed')    
    l=plt.axvline(x=mav+s1m[0], color='b', linestyle='dashed')
    l=plt.axvline(x=mav+s1m[1], color='b', linestyle='dashed')
    l=plt.axvline(x=mav+s3m[0], color='y', linestyle='dashed')
    l=plt.axvline(x=mav+s3m[1], color='y', linestyle='dashed')
    
    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)
    outname="hist_"+pname+".png"
    plt.savefig(outname)

    outfile="posterior_"+pname+".dat"
    of=open(outfile,'w')
    for ii in range(len(vals)):
        of.write('%14.8f\n' % (vals[ii]) ) 

    of.close()


def param_histbp(vals,pname,mv,s1v,s3v,mav,s1m,s3m,bpm,s1bpm):
    
    import numpy as np
    import matplotlib 
    matplotlib.use('Agg')
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt

    fig = plt.figure()

    matplotlib.rcParams.update({'font.size': 10})
    
    num_bins = 50
    # the histogram of the data
    n, bins, patches = plt.hist(vals, num_bins, facecolor='green', alpha=0.5)
    plt.xlabel(pname)
    plt.ylabel('N samples')
    l=plt.axvline(x=mv, color='r', linestyle='dotted')    
    l=plt.axvline(x=mv+s1v[0], color='b', linestyle='dotted')
    l=plt.axvline(x=mv+s1v[1], color='b', linestyle='dotted')
    l=plt.axvline(x=mv+s3v[0], color='y', linestyle='dotted')
    l=plt.axvline(x=mv+s3v[1], color='y', linestyle='dotted')

    l=plt.axvline(x=mav, color='r', linestyle='dashed')    
    l=plt.axvline(x=mav+s1m[0], color='b', linestyle='dashed')
    l=plt.axvline(x=mav+s1m[1], color='b', linestyle='dashed')
    l=plt.axvline(x=mav+s3m[0], color='y', linestyle='dashed')
    l=plt.axvline(x=mav+s3m[1], color='y', linestyle='dashed')

    l=plt.axvline(x=bpm, color='r')    
    
    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)
    outname="hist_"+pname+".png"
    plt.savefig(outname)

    outfile="posterior_"+pname+".dat"
    of=open(outfile,'w')
    for ii in range(len(vals)):
        of.write('%14.8f\n' % (vals[ii]) ) 

    of.close()

