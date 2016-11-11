import numpy as np
import pylab as pl
import string as s
import re
import specTools as t
import pdb
from os.path import expanduser

L_sun = 3.826E33 # the L_sun defined by BC03 in erg/s

def loadM05spec(fname, massfile=None, agecol=0, zcol=1, lamcol=2, fluxcol=3, \
                angstscale=1.0, fluxscale=1.0, skip=0, minAge=None, maxAge=None, \
                Zdict={"10m4":0.00025, "0001":0.001, "001":0.01, \
                "002":0.02, "004":0.04, "007":0.07}, \
                ZHdict={"10m4":-2.25, "0001":-1.35, "001":-0.33, "002":0.0, \
                        "004":0.35, "007":0.67}, \
                resolution=[None,{'vis':(5,10), 'nir':(20,100)}], ):
    """
    Author: Ryan Houghton (14/4/11)

    Updated with ZHdict to load masses for M11 models (11/11/2016)

    Purpose: To load in a stellar pop spectrum from M05, keeping
       - Age (Gyr)
       - Metallicity Z
       - Wavelength (AA)
       - Flux density (erg/s/AA/cm^2)

    Inputs:
       fname - the filename to load data from
       Zdict - a dictionary to work out model metallicity from filename
    """

    # get metallicity
    Ztag = fname.split(".")[1].split("z")[-1]
    Z = Zdict[Ztag]
    
    # read in the raw data
    alldata = np.real(np.loadtxt(fname, usecols=[agecol, zcol, lamcol, fluxcol], \
                    dtype="D", unpack=True)) # dtype="F" isn't good enough

    # get the array size and spec size 
    nx, ny = alldata.shape

    loc = (np.where(alldata[0,:]==alldata[0,0]))[0]
    nlam=len(loc)

    # check if file has whole number of spectra
    assert (float(ny)/float(nlam) % 1) == 0.0, "NLINE != N*NLAM ?!"

    nspec = ny/nlam

    # init
    specs = np.zeros((nspec,nlam))
    ages  = np.zeros(nspec)
    lam   = np.zeros(nlam)

    specs[0,:] = alldata[3, loc]
    lam[:]     = alldata[2, loc]
    ages[0]    = alldata[0, loc[0]] # in Gyrs
    
    # sort the other data
    for ix in range(1,nspec):
        loc = (np.where(alldata[0,:]==alldata[0,ix*nlam]))[0]
        specs[ix,:] = alldata[3, loc]
        ages[ix] = alldata[0,loc[0]]

    # convert: erg/s/AA => * 1.0 /(4.*pi*D[cm]**2) => erg/s/cm**2/AA (@10pc)
    factor= (1.0/(10.0*t.pc*100)**2.0) / (4.0*np.pi) 
    specs = specs * factor

    
    # load the stellar masses
    if massfile!=None:
        massdata = np.real(np.loadtxt(massfile, unpack=True, dtype="D"))
        zloc = np.where(massdata[0] == ZHdict[Ztag])[0] #alldata[1][0])[0]
        mass = massdata[2,zloc]
    else:
        mass=None

    # cut ages if given
    if (minAge is not None) & (maxAge is not None):
        aloc = np.where((ages>=minAge) & (ages<=maxAge))[0]
    elif (minAge is not None) & (maxAge is None):
        aloc = np.where(ages>=minAge)[0]
    elif (minAge is None) & (maxAge is not None):
        aloc = np.where(ages<=maxAge)[0]
    else:
        aloc = np.arange(len(ages))
    ages = ages[aloc]
    mass = mass[aloc]
    specs = specs[aloc,:]
    
    spectra = t.spectrum(lam=lam, lamspec=specs, age=ages, Z=Z, mass=mass, model="M05", \
                         resolution=resolution, wavesyst="air")

    return spectra
    

def loadM05ssps(sedpath="~/z/data/stellar_pops/M05", \
                imf="salpeter", glob="00[1,2,4]", morph="rhb", \
                IMFdict={"salpeter":"ss", "kroupa":"kr"}, \
                Zdict={"10m4":0.00025, "0001":0.001, "001":0.01, \
                "002":0.02, "004":0.04, "007":0.07}, multiD=True, verbose=False):
    """
    Author: Ryan Houghton (20/4/11)

    Purpose: Load the Maraston 2005 SEDs.

    Inputs:
       sedpath - location of the models
       imf     - [ salpeter | Krouper ] Type of IMF used: 2 choices only
       glob    - the REGEXP to select one/all the required metallicities/files 
       IMFdict - a dictionary to convert IMF choice to file name
       Zdict   - a dictionary to convert filenames to metallicities
       multiD  - defaults to True but for backward compatibility, set to false for list of spectra
    
    """

    # get the file name
    files, nfiles = t.findfiles(expanduser(sedpath)+"/"+"sed."+IMFdict[imf]+"z"+glob+"."+morph)

    # load the first to get the array size
    if verbose:
        print "Reading data from ..."
        for fname in files: print fname

    massfile=expanduser(sedpath)+"/"+"stellarmass."+imf

    specs = []
    specs.append(loadM05spec(files[0], massfile=massfile))
    if verbose: print "Read "+files[0]+" (Z="+str(specs[-1].Z)+")"

    # load in all the others
    for fname in files[1:]:
        specs.append(loadM05spec(fname, massfile=massfile))
        if verbose: print "Read "+fname+" (Z="+str(specs[-1].Z)+")"

    # put specs into multiD format
    if multiD:
        flams = []
        ages = []
        Zs = []
        masses=[]
        for s in specs:
            flams.append(s.flam)
            ages.append(s.age)
            masses.append(s.mass)
            Zs.append(np.tile(np.array(s.Z),s.nspec)) # Z is scalar, not array, so tile up
        # make multi-D spec
        specs = t.spectrum(lam=s.lam,lamspec=flams, age=ages, Z=Zs, mass=masses, model="M05", IMF=imf, \
                           resolution=s.resolution, wavesyst=s.wavesyst)

    return specs
    

    
