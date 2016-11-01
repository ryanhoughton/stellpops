import numpy as np
import pylab as pl
import pdb
import specTools as st
from scipy import interpolate
from os.path import expanduser
import atpy as at

basedir="/home/houghton/z/data/stellar_pops/"

L_sun = 3.826E33 # the L_sun defined by BC03 in erg/s

def loadCD12ssps(sedpath=basedir+"CvD12_v1.2", ageglob="t??.?_solar.ssp", massfile="mass_ssp.dat", model="CD12", Z=0.02, verbose=True):
    """
    Author:  Ryan Houghton (18/09/12)

    Purpose: Load one of the Conroy / van Dokkum Z=SOLAR, varying age, SSP models, keeping:
       - IMF (string)
       - Age (Gyr)
       - Metallicity Z
       - Wavelength (AA)
       - Flux density (erg/s/AA/cm^2)
    
    """
    # get the AGE files to read in
    afiles, nafiles = st.findfiles(expanduser(sedpath)+"/"+ageglob)

    # get mass info
    minfo = at.Table(expanduser(sedpath)+"/"+massfile, type="ascii")
    imftags = minfo.col1
    masses = np.array([minfo.col2, minfo.col3, minfo.col4, minfo.col5, minfo.col6, minfo.col7])

    # convert: L_sun/micron => * L_sun[erg/s] / 1e4 / (4.*pi*D[cm]**2) => erg/s/cm**2/AA (@10pc)
    factor= (L_sun/1e4/(10.0*st.pc*100.0)**2.0) / (4.0*np.pi)

    # read SSP files one by one
    ages=[]
    Zs=[]
    imfs=[]
    flux1=[]
    flux2=[]
    flux3=[]
    flux4=[]
    flux5=[]
    wave=None
    for f in afiles:
        # load table
        ssp  = at.Table(f, type='ascii')
        # get age of file
        ages.append(float(f[-14:-10]))
        # get metallicity
        Zs.append(Z)
        # get imfs of fluxes
        imfs.append(imftags) 
        # get wave (once)
        if wave==None: wave = ssp.col1
        # get fluxes for each IMF
        flux1.append(ssp.col2*factor)
        flux2.append(ssp.col3*factor)
        flux3.append(ssp.col4*factor)
        flux4.append(ssp.col5*factor)
        flux5.append(ssp.col6*factor)
        if verbose: print("Loaded "+f)


    flux=[]
    flux.append(np.array(flux1))
    flux.append(np.array(flux2))
    flux.append(np.array(flux3))
    flux.append(np.array(flux4))
    flux.append(np.array(flux5))
    # now make spectra for age variations, one for each IMF
    specs=[]
    for q in range(5):
        spec = st.spectrum(lamspec=flux[q], lam=wave, age=ages, mass=masses[:,q], \
                          Z=Zs[q], IMF=imftags[q], model=model, wavesyst="vac")
        specs.append(spec)

    return specs

def loadCD12varelem():
    """
    RH 28/10/2016
    Load the CvD spectra with varyine element abundances.
    """
    return loadCD12spec(basedir+"CvD12_v1.2/"+"t13.5_varelem.ssp")

def loadCD12afe():
    """
    RH 1/11/2016
    Load the CvD spectra with varying [alpha/Fe]
    """
    s02 = loadCD12spec(basedir+"CvD12_v1.2/"+"t13.5_afe+0.2.ssp")
    s03 = loadCD12spec(basedir+"CvD12_v1.2/"+"t13.5_afe+0.3.ssp")
    s04 = loadCD12spec(basedir+"CvD12_v1.2/"+"t13.5_afe+0.4.ssp")
    return [s02,s03,s04]

def loadCD12spec(filepath):
    """
    Originally written by Simon Zieleniewski.
    Adapted by Ryan Houghton. 

    Function to read in CvD12 SSP files and return spectra as a
        spectrum class (created by RCWH).

    Inputs:
    =======
    - filepath: Path and filename string of file for CvD spectra

    Outputs:
    ========
    - spectrum: A spectrum class for the given SSP SED.
                Initialised with units (lambda=A, flux=erg/s/cm2/A) @ D=10pc


    e.g
    
    csalpha= loadCD12spec(basedir+"CvD12_v1.2/t13.5_varelem.ssp")
    csabun = loadCD12spec(basedir+"CvD12_v1.2/t13.5_varelem.ssp")

    """  

    dat = np.genfromtxt(filepath)

    #Get filename
    fname = filepath.split('/')[-1]

    #Wavelenghts in A
    lambs = dat[:,0].copy()

    #Set flux units to erg/s/cm**2/A at D = 10 pc. CvD flux in units of L_sun/um
    flux = np.transpose(dat[:,1:].copy())
    factor = (L_sun/1e4/(10.0*st.pc*100.0)**2.0) / (4.0*np.pi)
    flux *= factor

    #Age of spectra in Gyrs
    Age = [float(fname.split('_')[0].split('t')[1])]*flux.shape[0]
    

    #Interpolate to get linear dispersion
    newlambs = np.linspace(lambs[0], lambs[-1], len(lambs))
    finterp = interpolate.interp1d(lambs, flux, kind='linear', axis=-1)
    newflux = finterp(newlambs)

    #Get mass file
    masspath = filepath.split(fname)[0]
    masses_orig = np.loadtxt(masspath+'mass_ssp.dat', dtype=np.str)
    masses = np.copy(masses_orig)
    masspos = {13.5:6, 11.0:5, 9.0:4, 7.0:3, 5.0:2, 3.0:1}
    mass = np.array(masses[:,masspos[Age[0]]], dtype=np.float)

    #Depending on filename, spectra correspond to different IMFs, ages etc
    if 'solar' in fname:
        #IMFs = x=3.5, 3.0, 2.35, Chabrier, bottom-light
        IMFs = ['x = 3.5', 'x = 3.0', 'x = 2.35', 'Chabrier', 'bottom-light']
        return st.spectrum(lamspec=newflux, lam=newlambs, age=Age,
                          Z=0.2, IMF=IMFs, model='CvD12', mass=mass, wavesyst="vac")

    if 'afe' in fname:
        met = 0.0
        IMFs = ['x = 3.5', 'x = 3.0', 'x = 2.35', 'Chabrier', 'bottom-light']
        afes = {'afe': float(fname.split('+')[1][0:3])}
        return st.spectrum(lamspec=newflux, lam=newlambs, age=Age, alpha=afes['afe'], 
                          Z=met, IMF=IMFs, model='CvD12', 
                          mass=mass, wavesyst="vac") #userdict=afes,

    if 'varelem' in fname:
        IMF = 'Chabrier'
        uAge = list(set(Age))[0]
        met = 0.0
        massloc = np.where(masses[:,0]==IMF)[0]
        masses = masses[massloc[0],1:]
        mass = float(masses[masspos[uAge]-1])
        
        abunlist = {'abundances': ['[Na/Fe] = +0.3', '[Na/Fe] = -0.3','[Ca/Fe] = +0.15', '[Ca/Fe] = -0.15',
                '[Fe/H] = +0.3', '[Fe/H] = -0.3', '[C/Fe] = +0.15', '[C/Fe] = -0.15',
                '[a/Fe] = +0.2', '[as/Fe] = +0.2', '[N/Fe] = +0.3', '[N/Fe] = -0.3',
                '[Ti/Fe] = +0.3', '[Ti/Fe] = -0.3', '[Mg/Fe] = +0.3', '[Mg/Fe] = -0.3',
                '[Si/Fe] = +0.3', '[Si/Fe] = -0.3']}
        return st.spectrum(lamspec=newflux, lam=newlambs, age=Age,
                          Z=met, IMF=IMF, model='CvD12', userdict=abunlist,
                          mass=mass, wavesyst="vac")

    else:
        raise ValueError('Did not input correct CvD12 file [as of 03-04-14]')
