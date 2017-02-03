import pylab as pl
import numpy as np
from stellpops import specTools as st
from stellpops import indexTools as it
import re
try:
    import astropy.io.fits as pf
except:
    try:
        import pyfits as pf
    except:
        raise ImportError("astropy.io.fits nor pyfits found.")
import pdb
import copy as cp
import glob as g
import atpy as at
import matplotlib as mpl
import scipy as si
import os
from stellpops.indexTools import drawGrid

fs=15
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['xtick.labelsize']=fs
mpl.rcParams['ytick.labelsize']=fs
mpl.rcParams['axes.titlesize']=fs
mpl.rcParams['axes.labelsize']=fs
mpl.rcParams['axes.labelweight']='normal'
mpl.rcParams['font.size']=fs
mpl.rcParams['font.weight']='normal'
mpl.rcParams['legend.fontsize']=int(fs*0.75)
mpl.rcParams['errorbar.capsize']=0


V15path='/home/houghton/z/data/stellar_pops/Vaz15'


def loadV15ssps(V15path=V15path, ageRange=None, ZRange=None, exactIMFs=None, exactAlphas=None, isochrone="BaSTI", IMFtype="un", verbose=True):
    """
    Ryan Houghton
    4/10/16

    
    Function to read in V15 SSP SEDs and return as spectrum class (modified V03tools.py).

    Inputs:
    - filepath: Path and filename string of FITS file for V03 SED
    - ageRange\_ Give the min and max values for age/Z
    - ZRange  /
    - exactIMFs   : give a list of IMFs to load: 'un:1.30'
    - exactAlphas : give a list of alphas to load: 'Ep0.00', 'Ep0.40', 'baseFe' 

    Outputs:
    - spectrum: A spectrum class for the given SSP SED.
                Initialised with units (lambda=A, flux=erg/s/cm2/A) @ D=10pc


    e.g.
    from stellpops import V15tools as v15
    # load all available models on disk
    vall = v15.loadV15ssps()
    # load a subset of the models:
    vsub = v15.loadV15ssps(ZRange=[0.0,0.3], exactAlphas=['Ep0.00','Ep0.40'], ageRange=[0.1,3.0], verbose=True)

    """

    if exactAlphas is not None:
        assert isinstance(exactAlphas, list), "exactAlphas is not a list"
    if exactIMFs is not None:
        assert isinstance(exactIMFs, list), "exactIMFs is not a list"

    pathToModels=V15path+"/"+isochrone+"/"+IMFtype+"/"

    files = g.glob(pathToModels+"M*.fits")
    nfiles = len(files)
    assert nfiles > 0, "Failed to find any files"
    if verbose: print "Found "+str(nfiles)+" files in V15 library"

    imftypeTAGs=[]; imfslopeTAGs=[]; imfTAGs=[]; alphatypeTAGs=[]; alphaTAGs=[]; ZTAGs=[];
    imftypes=[]; imfslopes=[]; imf=[]; Zs=[]; ages=[]; alphas=[]; alphatypes=[]

    fs=[]
    for flong in files:
        f=flong.split("/")[-1]
        # split the name up into model params
        splits=re.split('M|Z|T|_iTp|_|.fits',f)[1:-1]
        imftypes.append(splits[0][0:2])
        imfslopes.append(splits[0][2:])
        imf.append(imftypes[-1]+":"+str(imfslopes[-1]))
        pm=1.0
        if (splits[1][0]=="m"): pm=pm*-1.
        Z=float(splits[1][1:])*pm
        Zs.append(Z)
        ages.append(float(splits[2]))
        alphas.append(splits[3])
        alphatypes.append(splits[4])

        fs.append(f)

        ## single TAGs for IMF, alpha and Z
        #imfTAGs.append('I='+imf[-1])
        #alphaTAGs.append('A='+alphatypes[-1])
        #ZTAGs.append('Z='+str(Zs[-1]))

    # make unique lists, avoiding sets (which change order)
    uimfs=[]; [uimfs.append(i) for i in imf if not uimfs.count(i)]
    ualphas=[]; [ualphas.append(i) for i in alphatypes if not ualphas.count(i)]
    uZs=[]; [uZs.append(i) for i in Zs if not uZs.count(i)]
    uages=[]; [uages.append(i) for i in ages if not uages.count(i)]

    # sort IMFs, Zs
    uZs.sort()
    ualphas.sort()
    uimfs.sort()

    if verbose:
        print "Found IMFs : ", uimfs
        print "Found Z values:", uZs
        print "Found alpha types: ", ualphas
        print "Found ages ", uages

    # clip the parameter ranges if asked
    uages, ualphas, uZs, uimfs = clipParamRanges(uages, ualphas, uZs, uimfs, \
                                          ageRange=ageRange, ZRange=ZRange, exactAlphas=exactAlphas, exactIMFs=exactIMFs, verbose=verbose)

    nimfs = len(uimfs)
    nZs   = len(uZs)
    nalphas=len(ualphas)
    nages = len(uages)

    uages = np.array(uages)

    # load the files
    # init lists/arrays - 3D: alpha, Z, IMF + 1D of age => 4D
    aflams= [[[[] for z in range(nalphas)] for y in range(nZs)] for x in range(nimfs)]
    aages = [[[[] for z in range(nalphas)] for y in range(nZs)] for x in range(nimfs)]
    aZs   = [[[[] for z in range(nalphas)] for y in range(nZs)] for x in range(nimfs)]
    aalphas=[[[[] for z in range(nalphas)] for y in range(nZs)] for x in range(nimfs)]
    aIMFs = [[[[] for z in range(nalphas)] for y in range(nZs)] for x in range(nimfs)]
    aMass = [[[[] for z in range(nalphas)] for y in range(nZs)] for x in range(nimfs)]
    for ii, imf in zip(range(nimfs), uimfs):
        for Ai, A in zip(range(nalphas), ualphas): 
            # load the mass file
            if A[0]=="E":
                atag = A[1:]
            else:
                atag = A
            massfile = g.glob(V15path+"/masses/hssp_mass_"+isochrone+"*_"+imf[0:2].upper()+"_"+atag+"_v10.0.txt")[0]
            masses = at.Table(massfile, type='ascii')
            
            for Zi, Z in zip(range(nZs), uZs): 

            
                # load specs for all ages, fixed alpha, Z, IMF
                spec = loadV15specs(pathToModels, imf=imf, Z=Z, A=A, ageRange=ageRange, verbose=verbose)

                aflams[ii][Zi][Ai].append(spec.flam)
                aages[ii][Zi][Ai].append(spec.age)
                aalphas[ii][Zi][Ai].append(np.repeat(A,nages))
                aZs[ii][Zi][Ai].append(np.repeat(Z,nages))
                aIMFs[ii][Zi][Ai].append(np.repeat(imf,nages))
                # add masses from mass file
                aMass[ii][Zi][Ai].append(masses['M(*+remn)'][np.where((np.round(masses['[M/H]'],2)==Z) \
                                                                      & (masses.slope==float(imf[3:]))\
                                                                      & (masses.Age >= uages.min()) \
                                                                      & (masses.Age <= uages.max()))])
    
                
    specs = st.spectrum(lam=spec.lam, lamspec=np.array(aflams), \
                        age=np.array(aages), alpha=np.array(aalphas), mass=np.squeeze(np.array(aMass)), \
                        Z=np.array(aZs), IMF=np.array(aIMFs), resolution=[None,2.56], wavesyst="air", debug=False)
    # return the dict
    return specs

def clipParamRanges(uages, ualphas, uZs, uimfs, \
                    ageRange=None, ZRange=None, exactAlphas=None, exactIMFs=None, verbose=False):
    """
    Clip the parameter ranges if asked, to avoid always loading in the entire library.
    """

    agetest = ageRange is not None
    Ztest = ZRange is not None
    atest = exactAlphas is not None
    itest = exactIMFs is not None
    if agetest:
        uages = np.array(uages)
        aR=np.array(ageRange)
        ageloc = np.where((uages >= aR.min()) & (uages <=aR.max()))[0]
        uages = uages[ageloc]
    if Ztest:
        uZs = np.array(uZs)
        ZR=np.array(ZRange)
        Zloc = np.where((uZs >= ZR.min()) & (uZs <=ZR.max()))[0]
        uZs = uZs[Zloc]
    if atest:
        aloc=[]
        ualphas=np.array(ualphas)
        for aval in exactAlphas:
            aloc.append(np.where(ualphas==aval)[0][0])
        aloc=np.array(aloc)
        ualphas=ualphas[aloc]

    if itest:
        iloc=[]
        uimfs=np.array(uimfs)
        for ival in exactIMFs:
            iloc.append(np.where(uimfs==ival)[0][0])
        iloc=np.array(iloc)
        uimfs=uimfs[iloc]

    if (agetest or Ztest or atest or itest) and verbose:
        print "Loading IMFs : ", uimfs
        print "Loading Z values:", uZs
        print "Loading alpha types: ", ualphas
        print "Loading ages: ", uages

    return uages, ualphas, uZs, uimfs

def loadV15specs(V15path, imf='um:1.3', Z=0.0, A=0.0, ageRange=None, verbose=True):
    """
    Load the V15 specs for a given alpha, Z, IMF (i.e. ALL ages)
    """
    # make p or m for Z
    if Z>0.0:
        pm="p"
    elif Z<0.0:
        pm="m"
    elif Z==0.0:
        pm="p"
    else:
        raise ValueError("This should never happen.")

    # find the files for this IMF, Z, A and ALL ages
    itype, ival = imf.split(":")
    glob = "M"+itype+str(ival)+"*Z"+pm+"*"+str(np.fabs(Z))+"*T*"+"i?p*"+str(A)+"*.fits"

    files, nfiles = st.findfiles(V15path, glob=glob)

    assert nfiles>0, "No files found for this IMF, alpha and Z"
    if verbose: print "Found "+str(nfiles)+" files (ages) for IMF="+imf+", Z="+str(Z)+" and alpha="+str(A)

    # load first spectrum 
    flong=files[0]
    h=pf.open(flong, mode='readonly', memmap=False)
    naxis = h[0].header['NAXIS1']
    crval = h[0].header['CRVAL1']
    cdelt = h[0].header['CDELT1']
    crpix = h[0].header['CRPIX1']
    h.close()

    fluxes=[] # init
    theseages=[]
    for flong in files:
        # split filename
        f=flong.split("/")[-1]
        splits=re.split('M|Z|T|_iTp|_|.fits',f)[1:-1]
        # get info from filename
        imftype=splits[0][0:2]
        imfslope=splits[0][2:]
        pm=1.0
        if (splits[1][0]=="m"): pm=pm*-1.
        Z=float(splits[1][1:])*pm
        thisage=float(splits[2])

        # check if ageRange specified
        if ageRange is not None:
            aR = np.array(ageRange)
            # if ageRange given, only add spectra within this age range
            if (thisage >= aR.min()) & (thisage <= aR.max()):
                theseages.append(thisage)
                alpha=splits[3]
                alphatype=splits[4]
                # load file into memory
                h=pf.open(flong, mode='readonly', memmap=False)
                #sanity check
                assert (crval==h[0].header['CRVAL1'] and cdelt==h[0].header['CDELT1'] and \
                        crpix==h[0].header['CRPIX1'] and naxis==h[0].header['NAXIS1']), \
                        "This spectrum ("+f+") doesn't match the format of the first"
                # make flux array
                factor= st.L_sun/(4.*np.pi*(10.0*st.pc*100.0)**2.0)
                flux = h[0].data * factor
                fluxes.append(flux) # save flux
                h.close() # close file
        else:
           theseages.append(thisage)
           alpha=splits[3]
           alphatype=splits[4]
           # load file into memory
           h=pf.open(flong, mode='readonly', memmap=False)
           #sanity check
           assert (crval==h[0].header['CRVAL1'] and cdelt==h[0].header['CDELT1'] and \
                   crpix==h[0].header['CRPIX1'] and naxis==h[0].header['NAXIS1']), \
                   "This spectrum ("+f+") doesn't match the format of the first"
           # make flux array
           factor= st.L_sun/(4.*np.pi*(10.0*st.pc*100.0)**2.0)
           flux = h[0].data * factor
           fluxes.append(flux) # save flux
           h.close() # close file 

    # make wavelength array
    lam = crval + np.arange(naxis) * cdelt
    # make spectrum class for all these ages
    spec = st.spectrum(lamspec=fluxes, lam=lam, age=theseages, Z=Z, model="V15", IMF=imftype+str(imfslope), wavesyst="air")

    return spec

def calcV15M2L(filter='r_prime', imf="un", ageRange=None, ZRange=None, exactIMFs=None, exactAlphas=['baseFe'], isochrone="Padova"):
    """
    Calc the M2L for the V15 models, using either Padova or BaSTI isochrones

    NOTE: masses and M/L change a lot depending on isochrones
    
    """
    assert len(exactAlphas)==1, "More than one alpha value not allowed"
    quickFile = V15path+"/masses/hssp_m2l_"+filter+"_"+imf.upper()+"_"+isochrone+exactAlphas[0]+".txt"
    filePresent = os.path.exists(quickFile)
    if filePresent:
        tab = at.Table(quickFile, type='ascii')
    else:
        specs = loadV15ssps(V15path=V15path, ageRange=ageRange, ZRange=ZRange, exactIMFs=None, exactAlphas=exactAlphas, \
                            isochrone=isochrone, IMFtype=imf, verbose=True)

        filts = st.loadMyFilters()
        m2l = specs.calcM2L(filts[filter])
        tab = at.Table()
        tab.add_column("IMF", np.repeat(imf, np.prod(np.array(specs.dims))))
        tab.add_column("slope", np.array([x[3:] for x in specs.IMF.ravel()]), dtype=np.float)
        tab.add_column('(M/L)r', m2l.ravel())
        tab.add_column('[M/H]', specs.Z.ravel())
        tab.add_column('Age', specs.age.ravel())
        tab.write(quickFile, type="ascii")
    
    return tab

def loadV12M2L(imf="un", filterset="SDSS", set={"SDSS":"sdss", "JC":"phot"}, \
                    massfile="ssp_", version="9.2", massdir="masses", calcMass=False):
    """
    Load the M/L values for the MILES models

    """

    filename = V15path+"/"+massdir+"/"+massfile+set[filterset]+"_miuscat_"+imf.upper()+"_v"+str(version)

    tab = at.Table(filename, type="ascii")

    if calcMass:
        if filterset=="JC":
            Usol, Bsol, Vsol, Rsol, Isol = [5.60, 5.43, 4.82, 4.46, 4.15]
            MassU = tab['(M/L)U'] / 10.0**(-0.4*Usol) * 10.0**(-0.4*tab['U'])
            MassB = tab['(M/L)B'] / 10.0**(-0.4*Bsol) * 10.0**(-0.4*tab['B'])
            MassV = tab['(M/L)V'] / 10.0**(-0.4*Vsol) * 10.0**(-0.4*tab['V'])
            MassR = tab['(M/L)R'] / 10.0**(-0.4*Rsol) * 10.0**(-0.4*tab['R'])
            MassI = tab['(M/L)I'] / 10.0**(-0.4*Isol) * 10.0**(-0.4*tab['I'])
            Mass = np.mean(np.array([MassU, MassB, MassV, MassR, MassI]),axis=0)
        elif filterset=="SDSS":
            usol, gsol, rsol, isol = [6.55, 5.12, 4.68, 4.57]
            Massu = tab['(M/L)u'] / 10.0**(-0.4*usol) * 10.0**(-0.4*tab['u'])
            Massg = tab['(M/L)g'] / 10.0**(-0.4*gsol) * 10.0**(-0.4*tab['g'])
            Massr = tab['(M/L)r'] / 10.0**(-0.4*rsol) * 10.0**(-0.4*tab['r'])
            Massi = tab['(M/L)i'] / 10.0**(-0.4*isol) * 10.0**(-0.4*tab['i'])
            Mass = np.mean(np.array([Massu, Massg, Massr, Massi]),axis=0)
        else:
            raise ValueError("Filterset not understood")
            
        tab.add_column("mass", Mass)

        #tab2 = at.Table()
        #tab2.add_column('IMF', tab['IMF'])
        #tab2.add_column('slope', tab['slope'])
        #tab2.add_column('[M/H]', tab['[M/H]'])
        #tab2.add_column('Age', tab.Age)
        #tab2.add_column('Mass', Mass)
        #tab2.write('Masses_'+imf.upper()+".txt", type='ascii')
        #tab = tab2
        
    return tab

def loadJohnsonCousinsM2L(imf="un", massfile="ssp_phot_miuscat", version="9.2", massdir="masses"):
    """
    Load the JohnsonCousins M/L values for the MILES models

    """

    filename = V15path+"/"+massdir+"/"+massfile+"_"+imf.upper()+"_v"+str(version)

    tab = at.Table(filename, type="ascii")

    #Usol, Bsol, Vsol, Rsol, Isol = [5.60, 5.43, 4.82, 4.46, 4.15]
    #MassU = tab['(M/L)U'] / 10.0**(-0.4*Usol) * 10.0**(-0.4*tab['U'])
    #MassB = tab['(M/L)B'] / 10.0**(-0.4*Bsol) * 10.0**(-0.4*tab['B'])
    #MassV = tab['(M/L)V'] / 10.0**(-0.4*Vsol) * 10.0**(-0.4*tab['V'])
    #MassR = tab['(M/L)R'] / 10.0**(-0.4*Rsol) * 10.0**(-0.4*tab['R'])
    #MassI = tab['(M/L)I'] / 10.0**(-0.4*Isol) * 10.0**(-0.4*tab['I'])
    #
    #Mass = np.mean(np.array([MassU, MassB, MassV, MassR, MassI]),axis=0)
    #
    #tab2 = at.Table()
    #tab2.add_column('IMF', tab['IMF'])
    #tab2.add_column('slope', tab['slope'])
    #tab2.add_column('[M/H]', tab['[M/H]'])
    #tab2.add_column('Age', tab.Age)
    #tab2.add_column('Mass', Mass)
    #tab2.write('Masses_'+imf.upper()+".txt", type='ascii')
    return tab


def loadSDSSM2L(imf="un", massfile="ssp_sdss_miuscat", version="9.2", massdir="masses"):
    """
    Load the SDSS M/L values for the MILES models

    """

    filename = V15path+"/"+massdir+"/"+massfile+"_"+imf.upper()+"_v"+str(version)

    tab = at.Table(filename, type="ascii")

    #usol, gsol, rsol, isol = [6.55,5.12,4.68,4.57]
    #Massu = tab['(M/L)u'] / 10.0**(-0.4*usol) * 10.0**(-0.4*tab['u'])
    #Massg = tab['(M/L)g'] / 10.0**(-0.4*gsol) * 10.0**(-0.4*tab['g'])
    #Massr = tab['(M/L)r'] / 10.0**(-0.4*rsol) * 10.0**(-0.4*tab['r'])
    #Massi = tab['(M/L)i'] / 10.0**(-0.4*isol) * 10.0**(-0.4*tab['i'])
    #
    #Mass = np.mean(np.array([Massu, Massg, Massr, Massi]),axis=0)
    #
    #tab2 = at.Table()
    #tab2.add_column('IMF', tab['IMF'])
    #tab2.add_column('slope', tab['slope'])
    #tab2.add_column('[M/H]', tab['[M/H]'])
    #tab2.add_column('Age', tab.Age)
    #tab2.add_column('Mass', Mass)
    #tab2.write('Masses_'+imf.upper()+".txt", type='ascii')
    return tab
    

def investigateAZI(lamRange=[4000,6000], sigma=200.0):
    """
    Investigate how spectra change for alpha, Z and IMF changes at fixed age
    """

    vs_un  = loadV15ssps(ZRange=[0.0,0.3], exactAlphas=['Ep0.00','Ep0.40'], ageRange=[13.5,13.5], IMFtype="un", verbose=True)
    vs_un2 = vs_un.copy()
    vs_bi  = loadV15ssps(ZRange=[0.0,0.3], exactAlphas=['Ep0.00','Ep0.40'], ageRange=[13.5,13.5], IMFtype="bi", verbose=True)
    vs_bi2 = vs_bi.copy()
    
    from stellpops import indexTools as it
    libair=it.loadLickIndicesAir()
    vs_un.clipSpectralRange(lamRange[0],lamRange[1])
    vs_un.gaussVelConvolve(0.0,sigma)
    vs_un.normaliseSpec(polyOrder=5,indLib=libair)
    vs_un2.clipSpectralRange(lamRange[0],lamRange[1])
    vs_un2.gaussVelConvolve(0.0,sigma)

    vs_bi.clipSpectralRange(lamRange[0],lamRange[1])
    vs_bi.gaussVelConvolve(0.0,sigma)
    vs_bi.normaliseSpec(polyOrder=5,indLib=libair)
    vs_bi2.clipSpectralRange(lamRange[0],lamRange[1])
    vs_bi2.gaussVelConvolve(0.0,sigma)

    # find location of Salpeter and Chabrier
    uiloc = np.where(vs_un.IMF[:,0,0]=='un:1.30')[0][0]
    biloc = np.where(vs_bi.IMF[:,0,0]=='bi:1.30')[0][0]
    
    pl.figure()
    pl.subplot(311)
    su = vs_un2.flam[uiloc,0,1]/vs_un2.flam[uiloc,0,0]
    sb = vs_bi2.flam[biloc,0,1]/vs_bi2.flam[biloc,0,0]
    pl.plot(vs_un2.lam,su/np.median(su))
    pl.plot(vs_bi2.lam,sb/np.median(sb))
    pl.title("Alpha changes - Ratio")
    pl.subplot(312)
    su=vs_un2.flam[uiloc,2,0]/vs_un2.flam[uiloc,0,0]
    sb=vs_un2.flam[biloc,2,0]/vs_bi2.flam[biloc,0,0]
    pl.plot(vs_un2.lam,su/np.median(su))
    pl.plot(vs_bi2.lam,sb/np.median(sb))
    pl.title("Z changes - Ratio")
    pl.subplot(313)
    su=vs_un2.flam[-1,0,0]/vs_un2.flam[uiloc,0,0]
    sb=vs_un2.flam[-1,0,0]/vs_bi2.flam[biloc,0,0]
    sub=vs_un2.flam[uiloc,0,0]/vs_bi2.flam[biloc,0,0]
    hsub=vs_un2.flam[-1,0,0]/vs_bi2.flam[-1,0,0]
    pl.plot(vs_un2.lam,su/np.median(su))
    pl.plot(vs_bi2.lam,sb/np.median(sb))
    pl.plot(vs_un2.lam, sub/np.median(sub), "k-")
    pl.plot(vs_un2.lam, hsub/np.median(hsub), "k:")
    
    pl.title("IMF changes - Ratio")

    pdb.set_trace()
    
    pl.figure()
    pl.subplot(311)
    pl.plot(vs_un.lam,vs_un.flam[6,0,0])
    pl.plot(vs_un.lam,vs_un.flam[6,0,1])
    pl.title("Alpha changes - normalised")
    pl.subplot(312)
    pl.plot(vs_un.lam,vs_un.flam[6,0,0])
    pl.plot(vs_un.lam,vs_un.flam[6,2,0])
    pl.title("Z changes - normalised")
    pl.subplot(313)
    pl.plot(vs_un.lam,vs_un.flam[6,0,0])
    pl.plot(vs_un.lam,vs_un.flam[11,0,0])
    pl.title("IMF changes - normalised")
    

    pl.figure()
    pl.subplot(311)
    pl.plot(vs_un2.lam,vs_un2.flam[6,0,0]/np.median(vs_un2.flam[6,0,0]))
    pl.plot(vs_un2.lam,vs_un2.flam[6,0,1]/np.median(vs_un2.flam[6,0,1]))
    pl.title('Alpha changes - absolute')
    pl.subplot(312)
    pl.plot(vs_un2.lam,vs_un2.flam[6,0,0]/np.median(vs_un2.flam[6,0,0]))
    pl.plot(vs_un2.lam,vs_un2.flam[6,2,0]/np.median(vs_un2.flam[6,2,0]))
    pl.title('Z changes - absolute')
    pl.subplot(313)
    pl.plot(vs_un2.lam,vs_un2.flam[6,0,0]/np.median(vs_un2.flam[6,0,0]))
    pl.plot(vs_un2.lam,vs_un2.flam[11,0,0]/np.median(vs_un2.flam[11,0,0]))
    pl.title('IMF changes - absolute')


def compareMass2Light(band="r", age=12.5893, Z=0.0, minAge=3.0, maxAge=14.0, minZ=-0.5, showUN=True, showBI=True, trackMin=False):
    """
    Compare M/L plots for different IMFs
    
    """


    tabUN = loadV12M2L(imf="un", filterset="SDSS")
    tabBI = loadV12M2L(imf="bi", filterset="SDSS")

    pl.figure()
    # helper function to do plot
    def plotForAge(age, Z, label=None, showUN=True, showBI=True):
        locUN = np.where((tabUN['Age']==age) & (tabUN['[M/H]']==Z))[0]
        locBI = np.where((tabBI['Age']==age) & (tabBI['[M/H]']==Z))[0]

        if label is None:
            labelUN=None
            labelBI=None
        else:
            labelUN="UN "+label
            labelBI="BI "+label
        if showUN: pl.plot(tabUN['slope'][locUN], tabUN['(M/L)'+band][locUN], label=labelUN)
        if showBI: pl.plot(tabBI['slope'][locBI], tabBI['(M/L)'+band][locBI], ":", label=labelBI)

        # calc extreme values - for most top heavy, salpeter/chab and most bottom heavy
        UNv13 = np.where(tabUN['slope'][locUN]==1.3)[0][0]
        BIv13 = np.where(tabBI['slope'][locBI]==1.3)[0][0]

        # left val, right val, min val, salp/chab val
        UNex = [tabUN['(M/L)'+band][locUN][0], tabUN['(M/L)'+band][locUN][-1], \
                [tabUN['(M/L)'+band][locUN].min(), \
                 tabUN['slope'][locUN[tabUN['(M/L)'+band][locUN].argmin()]]], \
                tabUN['(M/L)'+band][locUN][UNv13]]
        BIex = [tabBI['(M/L)'+band][locBI][0], tabBI['(M/L)'+band][locBI][-1], \
                [tabBI['(M/L)'+band][locBI].min(), \
                 tabBI['slope'][locBI[tabBI['(M/L)'+band][locBI].argmin()]]], \
                 tabBI['(M/L)'+band][locBI][BIv13]]
        
        return UNex, BIex


    # if single age and z, do single plot
    if (age is not None) and (Z is not None):
        plotForAge(age, Z, showUN=showUN, showBI=showBI)
    # otherwise do range of plots, taking into account min and max
    else:
        ages = np.array(list(set(list(tabUN.Age))))
        ages.sort()
        if minAge is not None:
            cut = np.where(ages>=minAge)[0]
            ages = ages[cut]
        if maxAge is not None:
            cut = np.where(ages<=maxAge)[0]
            ages = ages[cut]
        zs = np.array(list(set(list(tabUN['[M/H]']))))
        if minZ is not None:
            cut = np.where(zs>=minZ)[0]
            zs=zs[cut]
        ages = [ages.min(), ages.max()]
        zs = [zs.min(), zs.max()]
        UNexs=[]
        BIexs=[]
        for a in ages:
            for z in zs:
                Ux, Bx = plotForAge(a, z, label=str(a)+" Gyrs, Z="+str(z), showUN=showUN, showBI=showBI)
                UNexs.append(Ux)
                BIexs.append(Bx)

    if trackMin:
        UNmin = []
        BImin = []
        for ux, bx in zip(UNexs, BIexs):
            UNmin.append(ux[2])
            BImin.append(bx[2])

        UNmin = np.array(UNmin)
        BImin = np.array(BImin)
        if showUN: pl.plot(UNmin[:,1], UNmin[:,0], "k-")
        if showBI: pl.plot(BImin[:,1], BImin[:,0], "k:")
    pl.xlabel('IMF Slope')
    pl.ylabel(r'M/L$_'+band+'$')
    pl.legend(loc=0)

    return UNexs, BIexs

def showMinMLs(minAge=1.):
    """
    Wrapper for below - used to create paper plot
    """

    pl.figure(figsize=(12,6))
    pl.subplot(121)
    pl.title("Unimodal")
    showMinML("un", minAge=minAge)
    pl.legend(loc=0)
    pl.xlabel(r'$\alpha_u$')
    pl.ylabel('M/L (r-band)')
    pl.subplot(122)
    pl.title("Bimodal")
    showMinML("bi", minAge=minAge)
    pl.xlabel(r'$\alpha_b$')
    

    pl.tight_layout()
    pl.savefig("MinML_IMFforms.png")
    
def showMinML(imf, minAge=None, maxAge=None, minZ=None, maxZ=None, alpha=0.4):
    """
    Plot the min ML, across all ages and Zs.
    
    """
    
    tabUN = loadM2L(imf=imf, filterset="SDSS")
    

    ages = np.array(list(set(tabUN.Age))) ; ages.sort()
    Zs = np.array(list(set(tabUN['[M/H]']))); Zs.sort()

    if minAge is not None: ages=ages[np.where(ages>=minAge)]
    if maxAge is not None: ages=ages[np.where(ages<=maxAge)]
    if minZ is not None: Zs=Zs[np.where(Zs>=minZ)]
    if maxZ is not None: Zs=Zs[np.where(Zs<=maxZ)]
    

    a0s=[]
    x=np.linspace(0.3,3.3,101)
    cols=['b', 'g', 'r', 'c', 'm', 'y', 'k']
    for nzi, z in enumerate(Zs):
        mins=[]
        for a in ages:
            UNloc = np.where( (tabUN.Age==a) & (tabUN['[M/H]']==z) )[0]
            
            spl = si.interpolate.UnivariateSpline(tabUN['slope'][UNloc], tabUN['(M/L)r'][UNloc], k=3, s=0.0)
            dspl = spl.derivative()
            xmin = x[np.fabs(dspl(x)-0.0).argmin()]

            #if z==Zs[-1]:
            if a==ages[-1]:
                #pl.plot(tabUN['slope'][UNloc], tabUN['(M/L)r'][UNloc], "k:", alpha=alpha)
                pl.plot(x,spl(x), cols[nzi]+'-', alpha=alpha)
                print a, z
                
            #rdspl = dspl.roots() # get minima 
            #UNminloc = tabUN['(M/L)r'][UNloc].argmin()
            #mins.append([tabUN['slope'][UNloc][UNminloc], tabUN['(M/L)r'][UNloc][UNminloc]])
            mins.append([xmin, spl(xmin)])
            print xmin
            
        mins=np.array(mins)
        loc = np.where(mins[1:,0]-mins[:-1,0]>0.)[0][0]
        a0 = np.mean(np.array([ages[loc],ages[loc+1]])) # age at which min ML changes
        a0s.append(a0)
        # can't get splines to work here
        #spl = si.interpolate.UnivariateSpline(mins[:,0], mins[:,1], k=3, s=10)
        #xx = np.linspace(mins[:,0].min(), mins[:,0].max(), 101)
        #pl.plot(xx, spl(xx), label='[M/H]='+"{:.2f}".format(z))

        pl.plot(mins[:,0], mins[:,1], label='[M/H]='+"{:.2f}".format(z))
        
    pl.yscale('log')
    
    return a0s # ages for which min ML transitions from one IMF to another


def createAgeZAlphaIndexGrid(indices=None, ageRange=None, ZRange=None, exactIMFs=None, exactAlphas=None, isochrone='BaSTI', IMFtype='un', \
                             sigma=None, atLickRes=False, calcMgFe=True, Fe1='Fe5270', Fe2='Fe5335', \
                             verbose=True, indexMethod=0, lib=None):
    """
    Create an index grid using the V15 models.
    
    """

    if sigma is not None:
        assert atLickRes is False, "Can't specify velocity broadening AND Lick resolution."

    if lib is None:
        # load models
        lib = loadV15ssps(V15path=V15path, ageRange=ageRange, ZRange=ZRange, exactIMFs=exactIMFs, \
                          exactAlphas=exactAlphas, isochrone=isochrone, IMFtype=IMFtype, verbose=verbose)

    # load AIR indices
    airInds = it.loadLickIndicesAir(atLickRes=atLickRes)
    if indices is not None:
        airInds = airInds.subset(indices)
        
    if sigma is not None:
        it.addSigma2IndLib(airInds, sigma=sigma)

    # calc indices
    vals={}
    for ind in airInds.names:
        vals[ind]=lib.calcIndex(airInds[ind], disp=None, method=indexMethod, verbose=False)

    if calcMgFe:
        # calc meanFe and MgFe index
        MgFe, avFe = it.calcMgFe(vals['Mg_b'], vals[Fe1], vals[Fe2], returnMeanFe=True)
        vals['<MgFe>']=MgFe
        vals['<Fe>']=avFe

    return vals, lib, airInds

def showZgrid(vals, lib, xindex, yindex, minage=0.5, maxage=4.0, minZ=-0.5, maxZ=0.5, alpha=0.0, color="k", \
              ageLabels=True, Zlabels=True, labelsize=10, ZforAgeLabel='min', ageForZlabel='max', \
              aha='right', ava='top', zha='right', zva='top', showAlpha=False, \
              ageZforAlphaLabel=['max','max'], Fe1='Fe5270', Fe2='Fe5335', tab=None, \
              alphalabs={'A=0.0':r'[$\alpha$/Fe]=0.0', 'A=0.3':r'[$\alpha$/Fe]=0.3', 'A=0.5':r'[$\alpha$/Fe]=0.5'}):
    """
    Use the previously calculated index grid to make an index-index plot
    """

    

    alphaKey = 'Ep'+'{0:.2f}'.format(alpha)

    # get original dims
    nZ, nAlpha, nAge = lib.dims

    # which samples to use along each dim
    alphaloc = np.where(lib.alpha[0,:,0]==alphaKey)[0]
    Zloc = np.where( (lib.Z[:,0,0]>minZ) & (lib.Z[:,0,0]<maxZ) )[0]
    ageloc = np.where( (lib.age[0,0,:]>minage) & (lib.age[0,0,:]<maxage) )[0]

    # new size
    nAlphaNew = len(alphaloc)
    nZnew = len(Zloc)
    nAgeNew = len(ageloc)

    # make address arrays
    alphai = np.tile(alphaloc,(nZnew,1,nAgeNew))
    Zi     = np.tile(np.tile(Zloc, (1,nAlphaNew,1)).T, (1,1,nAgeNew))
    agei   = np.tile(ageloc, (nZnew, nAlphaNew, 1))

    # make x and y arrays of indices.
    # If lib has variance (for bespoke models accounting for variance weighted indices),
    # then we ignore the errors calculated for the indices
    if lib.eflam is None:
        xarray = np.squeeze(vals[xindex][Zi, alphai, agei])
        yarray = np.squeeze(vals[yindex][Zi, alphai, agei])
    else:
        xarray = np.squeeze(vals[xindex][0, Zi, alphai, agei]) # 0 = index, 1 = error for index
        yarray = np.squeeze(vals[yindex][0, Zi, alphai, agei])


    ages = np.squeeze(lib.age[Zi,alphai,agei][0,0,:])
    Zs   = np.squeeze(lib.Z[Zi,alphai,agei][:,0,0])
    
    drawGrid(xindex, yindex, xarray, yarray, ages, Zs, alpha, color=color,
             ageLabels=ageLabels, Zlabels=Zlabels, labelsize=labelsize, \
             ZforAgeLabel=ZforAgeLabel, ageForZlabel=ageForZlabel, \
             aha=aha, ava=ava, zha=zha, zva=zva, showAlpha=showAlpha, \
             ageZforAlphaLabel=ageZforAlphaLabel, ageSuffix="Gyr", Zprefix='[M/H]', alphaPrefix=r'[$\alpha$/Fe]')
