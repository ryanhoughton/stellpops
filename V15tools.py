import pylab as pl
import numpy as np
from stellpops import specTools as st
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

V15path='/home/houghton/z/data/stellar_pops/Vaz15/'


def loadV15ssps(V15path=V15path, verbose=True):
    """
    Ryan Houghton
    4/10/16

    
    Function to read in V15 SSP SEDs and return as spectrum class (modified V03tools.py).

    Inputs:
    - filepath: Path and filename string of FITS file for V03 SED

    Outputs:
    - spectrum: A spectrum class for the given SSP SED.
                Initialised with units (lambda=A, flux=erg/s/cm2/A) @ D=10pc

    """

    files, nfiles = st.findfiles(V15path, glob='M*.fits')
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
        print "Found IMF : ", uimfs
        print "Found Z values:", uZs
        print "Found alpha types: ", ualphas

    nimfs = len(uimfs)
    nZs   = len(uZs)
    nalphas=len(ualphas)
    nages = len(uages)

    # load the files
    # init lists/arrays - 3D: alpha, Z, IMF + 1D of age => 4D
    aflams= [[[[] for z in range(nalphas)] for y in range(nZs)] for x in range(nimfs)]
    aages = [[[[] for z in range(nalphas)] for y in range(nZs)] for x in range(nimfs)]
    aZs   = [[[[] for z in range(nalphas)] for y in range(nZs)] for x in range(nimfs)]
    aalphas=[[[[] for z in range(nalphas)] for y in range(nZs)] for x in range(nimfs)]
    aIMFs = [[[[] for z in range(nalphas)] for y in range(nZs)] for x in range(nimfs)]
    for ii, imf in zip(range(nimfs), uimfs): 
        for Zi, Z in zip(range(nZs), uZs): 
            for Ai, A in zip(range(nalphas), ualphas): 
                # load specs for all ages, fixed alpha, Z, IMF
                spec = loadV15specs(V15path, imf=imf, Z=Z, A=A, verbose=verbose)

                aflams[ii][Zi][Ai].append(spec.flam)
                aages[ii][Zi][Ai].append(spec.age)
                aalphas[ii][Zi][Ai].append(np.repeat(A,nages))
                aZs[ii][Zi][Ai].append(np.repeat(Z,nages))
                aIMFs[ii][Zi][Ai].append(np.repeat(imf,nages))
                
    specs = st.spectrum(lam=spec.lam, lamspec=np.array(aflams), \
                        age=np.array(aages), alpha=np.array(aalphas), \
                        Z=np.array(aZs), IMF=np.array(aIMFs), resolution=[[None],[2.56]], wavesyst="air")
    
    # return the dict
    return specs

def loadV15specs(V15path, imf='um:1.3', Z=0.0, A=0.0, verbose=True):
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
    glob = "M"+itype+str(ival)+"*Z"+pm+"*"+str(np.fabs(Z))+"*T*"+"iTp*"+str(A)+"*.fits"

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
"""

Look at which bits of spectrum vary more for Z than for A:

pl.plot(vs['I=un:1.30']['Z=0.4']['A=Ep0.00'].flam[40]/vs['I=un:1.30']['Z=0.06']['A=Ep0.00'].flam[40], vs['I=un:1.30']['Z=0.06']['A=Ep0.00'].flam[40]/vs['I=un:1.30']['Z=0.06']['A=Ep0.40'].flam[40], ".", alpha=0.2)

pl.plot(x,x*1.35-0.07,"k-")


loc=np.where(vs['I=un:1.30']['Z=0.06']['A=Ep0.00'].flam[40]/vs['I=un:1.30']['Z=0.06']['A=Ep0.40'].flam[40] > (vs['I=un:1.30']['Z=0.4']['A=Ep0.00'].flam[40]/vs['I=un:1.30']['Z=0.06']['A=Ep0.00'].flam[40]*1.35-0.07))[0]

pl.figure()
pl.plot(vs['I=un:1.30']['Z=0.06']['A=Ep0.00'].lam,vs['I=un:1.30']['Z=0.06']['A=Ep0.00'].flam[40])
pl.plot(vs['I=un:1.30']['Z=0.06']['A=Ep0.00'].lam[loc],vs['I=un:1.30']['Z=0.06']['A=Ep0.00'].flam[40][loc], "ro")



"""
