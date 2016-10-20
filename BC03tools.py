import numpy as np
import pylab as pl
import string as s
import re
import pdb
import os
from stellpops import specTools as t
from os.path import expanduser

L_sun = 3.826E33 # the L_sun defined by BC03 in erg/s

def readsec(txt, startline):

    """
    Author: Ryan Houghton (16/4/11)

    Purpose: readlines() gives a massive string array; use this function to start from *startline*.
             This code is setup for the BC03 files, where the first number on the line indicates
             the number of entries to be read

    Returns:
       nval - number of values it was meant to read
       vals - all the values read (excluding above)
       n    - the line at which it stopped reading (nvals reached)

    """

    # get the first value on the line: number of elems to read
    line=s.split(txt[startline])
    nval=int(s.atof(line[0]))
    vals=[]
    for elem in line[1:]:
        vals.append(elem)

    # read upto nval values
    n=startline
    while len(vals) < nval:
        n+=1
        line=s.split(txt[n])
        for elem in line: vals.append(elem)

    vals=vals[0:nval]
    
    # convert to float
    nvals = []
    for val in vals: nvals.append(s.atof(val)) 

    return nval, nvals, n
    

def loadBC03spec(fname, IMF=None, loadMassFile=True, resolution=[None,{'3200AA-9500AA':(3.1,3.4)}]):

    """
    Author: Ryan Houghton (16/4/11)

    Purose: Read in a stellar pop spectrum from BC03, keeping
       - Age (Gyr)
       - Metallicity Z (one per file, given as arg)
       - Wavelength (AA)
       - Flux density (output is erg/s/AA even though input is L_sun/AA)

    Returns:
       nspec  - number of spectra it was meant to read
       ages   - the ages of the spectra
       nlam   - number of wavelength intervals for each spectrum 
       lams   - the wavelength of each spectrum pixel
       alldata- the data array containing spectra (nlam x nspec)
       
    """

    nlinelam = 636
    nlinedat = 1483

    # load spectrum file
    f=open(fname, "r")
    txt = f.readlines()
    f.close()
    nspec, ages, stopline = readsec(txt, 0)
    # advance to lines starting "^Padova"
    linenum=stopline
    while ( (re.search("^Padova", txt[linenum]) is None) & (linenum<len(txt)-1) ): linenum+=1
    specsStartLine=linenum+3
    Zline=linenum+2
    nlam, lams, stopline = readsec(txt, specsStartLine)
    # load complementary mass file
    if loadMassFile:
        subnames = s.split(fname,".")
        mname = s.join(subnames[:-1],sep="")+".4color"
        mages,mbol,bmag,vmag,mlb,mlv,mstar,mgas,mgal,sfr = np.loadtxt(mname, unpack=True)
        mass = mstar # in M_sun units
        mass = np.insert(mass, 0, mass[0]) # add extra t=0 element
    else:
        mass=None

    Z = s.atof(txt[Zline].strip().split("=")[-1])

    ages=np.array(ages) / 1e9 # put into Gyrs
    lams=np.array(lams)

    # note the order here: nlam seems to be on less efficient dimention... but it's not, numpy is odd
    alldata = np.zeros((nspec,nlam))

    linenum = stopline
    for sn in range(nspec):
        while ( (re.search("^ "+str(int(nlam)), txt[linenum]) is None) & (linenum<len(txt)-1) ): linenum+=1
        ndat, data, stopline = readsec(txt, linenum)
        linenum=stopline
        #print sn, linenum
        alldata[sn,:] = np.array(data)

    print ""

    # convert: L_sun/AA => * L_sun[erg/s] /(4.*pi*D[cm]**2) => erg/s/cm**2/AA (@10pc)
    # note that D=10pc!
    factor= (L_sun/(10.0*t.pc*100.0)**2.0) / (4.0*np.pi)
    alldata = alldata * factor

    spectra = t.spectrum(lamspec=alldata, lam=lams, age=ages, mass=mass, Z=Z, model="BC03", IMF=IMF, resolution=resolution)

    return spectra #nspec, ages, nlam, lams, alldata

        
def loadBC03ssps(sedpath="~/z/data/stellar_pops/BC03/models", \
                 tracks="Padova1994", imf="salpeter", glob="*hr*{52,62,72}*ssp.ised_ASCII", \
                 Zdict = {"m22":0.0001,"m32":0.0004,"m42":0.004,"m52":0.008, \
                 "m62":0.02, "m72":0.05}, Zind=3, verbose=False):

    """
    Author: Ryan Houghton (18/4/11)

    Purpose: Read in the SEDs from the standard BC03 SSP SEDs and convert to erg/s/cm**2/AA at D=10pc

    Inputs: 
       sedpath - the path to the directory containing the seds
       tracks  - the directory in the above that describes the tracks used
       imf     - the directory in the above that describes the IMF used
       glob    - the REGEXP to select the SEDs. Note 'hr', 'ised' may be useful
       Zdict   - a dictionary describing the link between 'm??' and the metalicity
       Zind    - if the sed filename is split by '_', then this is the index of the metallicity
                 code (e.g. 'm42') used in the above
       verbose - set to False for verbose output
    
    Returns:
       alldata - array of all spectra (nlam x nages x nZ) in units erg/s/cm**2/AA (at D=10pc)
       lams    - array of wavelength values
       ages    - array of ages of each spectrum
       Z       - array of metallicities

    Examples:
       > from BC03tools import loadBC03ssps
       > specs = loadBC03ssps()
       
    """

    files, nfiles = t.findfiles(expanduser(sedpath)+"/"+tracks+"/"+imf+"/"+glob)


    if verbose:
        print "Reading data from..."
        for fname in files: print fname

    # read in 1st file to get sizes
    spec = loadBC03spec(files[0], IMF=imf)
    Zcode = files[0].split("_")[Zind]
    if verbose: print "Read "+files[0]+" (Z="+str(Zdict[Zcode])+")"

    # init: note the apparent strange order: nlam last. This IS the most efficient way
    specs = []
    specs.append(spec)

    # now cycle through all files, loading one at a time,
    # using dict to get Z
    count=1
    for fname in files[1:]:
        
        spec = loadBC03spec(fname)
        Zcode = fname.split("_")[Zind]
        if verbose: print "Read "+fname+" (Z="+str(Zdict[Zcode])+")"
        specs.append(spec)
        count+=1

    return specs


def loadFilters(dir="~/z/data/stellar_pops/BC03/src/", file="filterfrm.res", \
                verbose=False):

    """
    Author: Ryan houghton (16/4/11)

    Purpose: Read in the BC03 filters into a dictionary

    Inputs:
       dir  - where the filter file is kept
       file - the name of the file
       verbose - if True, print more info about process

    Returns:
       filters - a dictionary of the filters in the file. Each key is linked to
                 a 2D array of lambda (AA) and transmission


    Examples:
       > from BC03tools import loadFilters
       > filters = loadFilters()
       > filters[filters.keys()[0]]  # print the first filter
    """

    # open file for reading, read, and then close file
    f = open(expanduser(dir)+file,'r')
    txt = f.readlines()
    f.close()

    # define num of lines of txt
    nlines=len(txt)

    # init
    count=1
    lams=[]
    trans=[]
    filters={}

    # get 1st filter name
    filt = s.join((txt[0].split())[1:])
    # read in 1st filter
    while (re.search("^#",txt[count])==None) & (count < nlines):
        lam, tran = txt[count].split()
        lams.append(float(lam))
        trans.append(float(tran))
        count+=1

    # now read in other filters in a loop
    while count < nlines:
        if re.search("^#",txt[count])==None:
            lam, tran = txt[count].split()
            lams.append(float(lam))
            trans.append(float(tran))
        else:
            # convert to arrays
            lams=np.array(lams)
            trans=np.array(trans)
            # check regularity of lam spacing
            #first = lams[0]
            #second= lams[1]
            #last  = lams[-1]
            #length= float(len(lams)-1)
            #delta = (last-first)/ length
            #if (second-first) != delta:
            #    print "FILTER SAMPLING NOT REGULAR: CORRECTING..."
            #    newlams = np.arange(first,last+delta,delta)
            #    newtrans = t.alinterp(lams,trans,newlams)
            #    
            #    lams=newlams
            #    trans=newtrans
            #    #raise "OH NO! Filter sampling is not regular" 
            #trans /= np.sum(trans) # normalise
                        
            filters[filt]=t.spectrum(lamspec=trans, lam=lams, filter=True, model="BC03 Filters")
            if verbose: print "Done "+str(filt)
            # start the next filter
            filt = s.join((txt[count].split())[1:])
            lams=[]
            trans=[]
            
        count+=1

    # add the last filter to the set
    lams=np.array(lams)
    trans=np.array(trans)
    filters[filt]=t.spectrum(lamspec=trans, lam=lams, filter=True, model="BC03 Filters")
    if verbose: print "Done "+str(filt)

    return filters


def loadVega(dir="~/z/data/stellar_pops/BC03/src/", \
             fname="A0V_KURUCZ_92.SED"):

    """
    Author: Ryan Houghton (16/4/11)

    Purpose: read the VEGA spec given in BC03 and return.

    Inputs:
       dir   - the dir where the file is found
       fname - the filename of the VEGA spec

    Returns:
       lams  - wavelength of the pixel (AA)
       flux  - flux in each pixel (erg/s/cm**2/AA)
    """

    lams, flux = np.loadtxt(expanduser(dir)+fname, skiprows=2, unpack=True)

    # convert flux from L_sun/AA to erg/s/cm**2/AA at D=10pc
    factor = (L_sun/(10.0*t.pc*100)**2.0) / (4.0*np.pi)
    flux  *= factor
    
    vega = t.spectrum(lams=lams,lamspec=flux, model="VEGA")
    
    return vega
        

def loadBC03indexfile(file, keyline=28, startline=31, verbose=False):
    """
    Purpose: To read in the BC03 spectral indicies from a single file

    Inputs:
       file      - a string filename
       keyline   - the line in the file that gives the names of the indicies (prefixed with '#')
       startline - the line in the file where the data starts flowing
       verbose   - if True, print info about process

    Returns:
       a dictionary of the indicies and for specific ages

    Example:

    > indicies = loadBC03indexfile()
    > pl.plot(indicies['log-age'], indicides['H\\beta'])
    
    """
    
    # load the first file to get the shape of things to come
    f = open(file)
    txt = f.readlines()
    f.close()

    # key descriptors for dictionary
    keys = txt[keyline].split()[1:]

    # number of rows and cols in file
    nrow = len(txt)
    ncol = len(keys)

    # init array for indicies
    indarray = np.zeros((ncol,nrow-startline)) 

    for row in range(startline,nrow):
        line = txt[row].split()
        for col in range(ncol):
            indarray[col,row-startline] = s.atof(line[col]) 
        
    indicies = {}

    for col in range(ncol):
        indicies[keys[col]] = np.squeeze(indarray[col,:])

    if verbose: print "Loaded indicies from file "+file
    return indicies


def loadBC03ind(indpath="~/z/data/stellar_pops/BC03/models", \
                tracks="Padova1994", imf="salpeter", glob="*hr*ssp.[6,7]lsindx_sed", \
                Zdict = {"m22":0.0001,"m32":0.0004,"m42":0.004,"m52":0.008, \
                "m62":0.02, "m72":0.05}, Zind=3, keyline=28, startline=31):

    """
    Purpose: Read in the spectral indicies from the BC03 models

    Inputs:

    Returns:
       a dictionary of the indicies (ordered at first in terms of metallicity and
       then interms of the different indicies - including log-age)

    Example:

    > indicies = loadBC03ind()
    > pl.plot(indicies['0.02']['log-age'],indicies['0.02']['Mg-b'])
    
    """

    files, nfiles = t.findfiles(expanduser(indpath)+"/"+tracks+"/"+imf+"/"+glob)

    indices = {}
    for n in range(nfiles):
        # get the metallicity
        Zcode = files[n].split("_")[Zind]
        key = str(Zdict[Zcode])
        if key in indices:
            indices[key].update(loadBC03indexfile(files[n], keyline=keyline, startline=startline))
        else:
            indices[key] = loadBC03indexfile(files[n], keyline=keyline, startline=startline)

    return indices

def testChabSalpIMF(f='f625w', z=0.545):
    '''
    Test how M/L varies for Chabrier and Salpeter IMFs

    Ratio is M2L_Salp/M2L_chab ~ 1.8, or M2L_chab/M2L_Salp ~ 0.55, but it changes with age.
    
    '''

    
    bcs_salp = bc.loadBC03ssps(imf='salpeter')
    bcs_chab = bc.loadBC03ssps(imf='chabrier')

    af=st.loadACSFilters()

    m2l_chab = bcs_chab[1].calcM2L(af[f], z=z)
    m2l_salp = bcs_salp[1].calcM2L(af[f], z=z)

    pl.plot(bcs_salp[0].age, m2l_salp/m2l_chab)
    pl.xlabel("Age")
    pl.ylabel(r'(M/L)$_{S}$ / (M/L)$_{C}$')
    
