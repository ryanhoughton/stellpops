import numpy as np
import pylab as pl
import string as strg
import scipy.interpolate as ip
from numpy.polynomial import chebyshev as cheby
import BC03tools as bc
import M05tools as ma
import re
import os
import pdb
from os.path import expanduser
from nearest import nearest as nn


# globals

c = 299792458.0 # m/s
pc= 3.08568025E16
d_sun = 149597887.5e3 # m

class spectrum:
    """
    Author: Ryan Houghton (20/4/11)

    Purpose: A spectrum object to aid manipulation of wavelength/frequency and the associated
             flux values

    Inputs:
    This class MUST be initalised with EITHER
       lamspec - a NLAMxNSPEC numpy array with [:,0] being the wavelength array (AA) and [:,1:] being
                 the associated flux values (erg/s/cm**2/AA)
       muspec  - a NLAMxNSPEC numpy array with [:,0] being the frequency array (GHz) and [:,1:] being
                 the associated flux values (GJy = 1e-14 erg/s/cm**2/Hz)

    NOTE: 2xN seems a funny way to order the indicies but in numpy, this IS the efficient
    way to order the memory elements

    Definitions:
    Once you have initalised the class with either lamspec or muspec, you will find available
    the following variables (no matter if you initalised with lamspec or muspec):
       lam  - the wavelength array (AA)
       flam - the associated flux (erg/s/cm**2/AA)
       mu   - the frequency array (GHz)
       fmu  - the associated flux (GJy = 1e-14 erg/s/cm**2/Hz)

    Functions:
       calcABmag  - given a filter (spectrum class), return a magnitude on the AB system
       calcSTmag  - as above but for the ST system
       calcVEGAmag- as above but for the VEGA system; you must also supply a vega *spectrum*

    Notes:
       - You can define more than one flam/fmu for each lam/mu
       - Filters used in the functions should also be *spectrum* classes
       - Filters don't need to be sorted in mu or lam
       
    """

    global c, pc, d_sun
    
    def __init__(self, lam=None, lamspec=None, errlamspec=None, mu=None, muspec=None, errmuspec=None, \
                 age=None, mass=None, Z=None, IMF=None, filter=False, model=None, resolution=None, \
                 userdict=None):

        """
        Inputs:
           lam     - the 1D array of wavelengths (becomes self.lam)
           lamspec - the (1D, 2D or 3D) array of fluxes on the specified lambda grid (becomes self.flam)
           errlamspec - standard deviation for pixels in lamspec (becomes self.eflam)
           mu      - the 1D array of frequencies (becomes self.mu)
           muspec  - the (1D, 2D or 3D) array of fluxes on the specified mu grid (becomes self.fmu)
           errmuspec - standard deviation for pixels in muspec (becomes self.efmu)
           age     - the 1D array of ages (Gyr) for the spectra
           mass    - the 2D array of stellar masses (M_Sun) at each age for each spectrum
           Z       - the 1D array of metallicities for the spectra
           IMF     - string describing IMF type
           filter  - if the spectrum is a filter transmission (photon %), then set this to TRUE. This affects how we calc Fmu
           model   - model name (e.g. M05, BC03, CD12, ATRAN)
           resolution - 2 element list of
               [0]: spectral R=lam/dlam (FWHM),
               [1]: resolution FWHM (val/tuple/dict) in AA
                    - may be single value or tuple defining starting and final resolutions.
                    - may be a dict for 'visible' and 'nir' wavelengths if multiple libs used
               Either may be None, depending on how the models/libs are defined.
               
           userdict- add extra user-defined tags to the spectrum, such as airmass, water vapour etc.

        Notes:
        There are 3 ways to inialise a spectrum class:
            1. With a 1D lam/mu and a 1D spectrum (single spec mode).
                 - Magnitudes are returned as single values
            2. With a 1D lam/mu and a 2D spectrum array (NSPECxN{LAM/MU}) (multispec mode)
                 - In this case, you may specify AGE where len(age)=NSPEC.
                 - Magnitudes are returned as 1D arrays with NSPEC=len(age) elements
            #3. With a 1D lam/mu and a 3D spectrum array (NZxNAGExN{LAM/MU}) (multispec mode)
            #     - In this case, you can specify AGE and Z.
            #     - Magnitudes will be returned as 2D arrays with NZxNAGE elements
           
        """

        # check if initalised "correctly"
        if (lamspec is None) & (muspec is None):
            raise "Spectrum must be defined with either lamspec(AA) or muspec (Hz)"
        if (lamspec is not None) & (muspec is not None):
            raise "Please don't define BOTH lamspec and muspec: just one or t'other"

        # start defining spectrum parts
        if lamspec is not None:
            # check that lam has been given
            if lam is None: raise "If you give lamspec, you must also give lam"

            # make sure 2d
            flam = np.atleast_2d(lamspec)
            # get array size
            loc = flam.shape
            
            # check for bigger arrays
            if len(loc)> 2: raise "lamspec not understood"

            # get sizes
            nlam = loc[1]
            nspec= loc[0]

            self.lam  = lam
            self.flam = flam#.tolist()

            if errlamspec is not None:
                eflam = np.atleast_2d(errlamspec)
                eloc = eflam.shape
                self.eflam = eflam
                # sanity check
                assert np.all(loc==eloc), "Flux and error arrays appear different sizes..."

            self.calcmuspec(filter=filter)

        if muspec is not None:
            if mu is None: raise "If you give muspec, you must also give mu"

            # if 1D array, blk up to 2D
            fmu = np.atleast_2d(muspec)
            loc = fmu.shape
            if len(loc)!= 2: raise "muspec not understood"

            # get sizes
            nmu   = loc[1]
            nspec = loc[0]
            
            self.mu   = mu
            self.fmu = fmu#.tolist()

            if errmuspec is not None:
                # if 1D array, blk up to 2D
                efmu = np.atleast_2d(efmu)
                eloc = efmu.shape
                self.efmu = efmu
                # sanity check
                assert np.all(loc==eloc), "muspec and errmuspec seem to be different shapes."

            self.calclamspec(filter=filter)

        # add age info
        if age is not None:
            if(len(age)!=nspec): raise ValueError("NAGE != NSPEC?!")
            self.age=age
            self.logage = np.log10(age)

        # add stellar mass info for each age
        if mass is not None:
            #if(len(mass)!=nspec): raise "NMASS != NSPEC?!"
            self.mass=mass
            self.logmass = np.log10(mass)

        # add metallicitiy
        if Z is not None:
            if len(np.array([Z]).shape)!=1: raise ValueError("Metallicity Z must be a scalar")
            self.Z=Z

        # add IMF
        if IMF is not None:
            self.IMF=IMF

        # name of the model, e.g. BC03
        if model is not None:
            self.model=model

        # add the resolution in AA
        if resolution is not None:
            self.resolution=resolution 

        # add user dictionary for extra info
        if userdict is not None:
            keys = userdict.keys()
            for key in keys:
                setattr(self, key, userdict[key])

    def calcmuspec(self, filter=False):
        """
        Calculate a muspec (GJy @ GHz) from a lamspec (erg/s/cm**2/AA @ AA)
        """
        # (convert AA to GHz with c in m/s): (c[m/s] / (lam[AA] * 1e-10) / 1e9 = c/AA * 1e1)
        self.mu = c / (self.lam) * 1e1
        # erg/s/cm**2/AA => *LAM[AA]*LAM[M]  / C[M/S] * 1e23 = LAM[AA]*LAM[AA]/ c[m/s] * 1e13 => GJy=1e-14 erg/s/cm**2/Hz (Jy=1e-23erg/s/cm2/Hz=1e-26W/m2/Hz)

        self.fmu = []
        if filter:
            for flam in self.flam:
                self.fmu.append( flam )
            if hasattr(self, 'eflam'):
                self.efmu = []
                for eflam in self.eflam:
                    if eflam is not None:
                        self.efmu.append( eflam )
                    #else:
                    #    self.efmu.append(None)
        else:
            for flam in self.flam:
                self.fmu.append( flam * self.lam * (self.lam) / c * 1e4 )
            if hasattr(self, 'eflam'):
                self.efmu = []
                for eflam in self.eflam:
                    if eflam is not None:
                        self.efmu.append( eflam * self.lam * (self.lam) / c * 1e4 )
                    #else:
                    #    self.efmu.append(None)
            


    def calclamspec(self, filter=False):
        """
        Calculate a lamspec (erg/s/cm**2/AA @ AA) from a muspec (GJy @ GHz)
        """
        # (convert GHz to AA with c in m/s): c[m/s] / (mu[GHz]*1e9) * 1e10 = c[m/s] / mu[GHz] * 1e1
        self.lam = c / (self.mu) * 1e1
        # GJy => / LAM(AA), /LAM(M)  * C(M/S) => erg/s/cm**2/Hz

        self.flam = []

        if filter:
            for fmu in self.fmu:
                self.flam.append( fmu )
            if hasattr(self, 'efmu'):
                self.eflam = []
                for efmu in self.efmu:
                    if efmu is not None:
                        self.eflam.append( efmu )
                    else:
                        self.eflam.append(None)
            
        else:
            for fmu in self.fmu:
                self.flam.append( fmu / self.lam / (self.lam) * c / 1e4 )
            if hasattr(self, 'efmu'):
                self.eflam=[]
                for efmu in self.efmu:
                    if efmu is not None:
                        self.eflam.append( efmu / self.lam / (self.lam) * c / 1e4 )
                    else:
                        self.eflam.append(None)

    def rebinLam(self, dLam=1e-4, flux=False):
        """
        Rebin the spectrum to a fixed/minimum dLam
        """

        newflams=[]
        newlam = np.arange(nn(self.lam.min(),dLam,ceil=True), nn(self.lam.max(),dLam,floor=True)+dLam, dLam)
        for flam in self.flam:
            newflam = rebin(self.lam,flam,newlam,flux=False)
            newflams.append(newflam)

        self.lam=newlam
        self.flam=newflams
        self.calcmuspec()


    def calcABmag(self, filter, z=0.0, plot=False, bandcor=False):
        """
        Purpose: use a filter and a spectrum to calculate the AB magnitude. 

        Inputs:
           filter  - a Spectrum class of the filter profile
           z       - the redshift of the observation (applies cosmol stretch)
           plot    - if True, show the spec and filter profile
           bandcor - the K-correction usually consists of a term for the different spectral region
                     and a term for the squashing of the spectrum by the redshift. If you set *bandcor*
                     to False, we do not apply the latter term (brighten mags by 1+z) because this correction
                     can be easily applied to the measured magnitudes (i.e. you can make them fainter by 1+z)

        Returns:
           AB magnitude for each spectrum
           
        """

        # check object type of filter: make sure Spectrum
        if isinstance(filter,object)==False: raise "Filter is not a class"
        if len(filter.fmu) > 1: raise "Filter has more than one throughput curve (multispec)"
        
        # interpolate the spec onto the filter curve
        # find the first spec point before filter point
        startmu = np.max(self.mu[np.where(self.mu < np.min(filter.mu*(1.0+z)))[0]])
        stopmu  = np.min(self.mu[np.where(self.mu > np.max(filter.mu*(1.0+z)))[0]])

        # get the max number of points - either from number of spec
        # points or number of filter points (whichever greater)
        np_spec = float(len(np.where((self.mu > startmu) & (self.mu < stopmu))[0]))
        np_filt = float(len(filter.mu))
        npoints = max([np_spec*3,np_filt*3])

        # interpolate
        mu = np.linspace(startmu,stopmu,npoints)
        ithro = interpolate(filter.mu*(1.0+z), filter.fmu[0], mu)

        # cycle through spectra, calculating mag for each
        mag = []
        for fmu in self.fmu:
            ifmu = interpolate(self.mu, fmu, mu)
            
            # integrate using simple trapezium method
            flux = np.sum(ifmu*ithro/mu) / np.sum(ithro/mu) # this normalisation by the filter curve means mags are bandpass corrected (bandpass term of K-cor)
            
            # calc mag: but remember that flux in GJy (mu in GHz has no effect with trapezium integration)
            if bandcor:
                mag.append(-2.5*np.log10(flux) - 48.6 - 2.5*np.log10(1e-14) - 2.5*np.log10(1.0+z)) # this removes bandpass correction (makes mags brighter)
            else:
                mag.append(-2.5*np.log10(flux) - 48.6 - 2.5*np.log10(1e-14) )

            if plot:
                pl.plot(mu,ifmu/np.median(ifmu))
                pl.plot(mu,ithro/np.max(ithro))

        return np.array(mag)

    def quickABmag(self, filter, z=0.0, plot=False, bandcor=False):
        """
        Purpose: try to calculate magnitudes faster (for multispec)

        Inputs:
           filter -
           z      - the redshift of the observation (applies cosmol stretch)
           plot   -
           bandcor- the K-correction usually consists of a term for the different spectral region
                    and a term for the squashing of the spectrum by the redshift. If you set *bandcor*
                    to False, we do not apply the latter term (brighten mags by 1+z) because this
                    can be easily applied to the measured magnitudes

        Returns
           AB magnitude for each spectrum
        """

        # check object type of filter: make sure Spectrum
        if isinstance(filter,object)==False: raise "Filter is not a class"
        if len(filter.fmu) > 1: raise "Filter has more than one throughput curve (multispec)"
        
        # interpolate the filter curve onto the spec
        ithro = np.array(interpolate(filter.mu*(1.0+z), filter.fmu[0], self.mu, method=2))
        loc = np.where(ithro > 0.0)[0]
        lithro = ithro[loc]
        dmu = (self.mu[loc]-self.mu[loc+1]) # be sure to correct for uneven mu spacing
        lithro *= dmu/self.mu[loc] # divide by mu here as it occurs in integral - Hogg k-cor

        
        # calc mags
        specs = np.array(self.fmu)
        fluxes = np.dot(specs[:,loc],lithro) / np.sum(lithro) # this normalisation is essentially bandpass term of K-correction

        # calc mag: but remember that flux in GJy (mu in GHz has no effect with trapezium integration)
        if bandcor:
            mags = -2.5*np.log10(fluxes) - 48.6 - 2.5*np.log10(1e-14) - 2.5*np.log10(1.0+z) # this puts bandpass term back in (i.e. equiv to raw obs mag)
        else:
            mags = -2.5*np.log10(fluxes) - 48.6 - 2.5*np.log10(1e-14) 

        if plot:
            loc = np.where( (self.mu > np.amin(filter.mu)*(1.0+z)) & (self.mu < np.amax(filter.mu)*(1.0+z)) )
            s = np.sum(specs, axis=0)
            pl.plot(c/self.mu[loc]/1e9/1e-6, s[loc]/np.median(s[loc]))
            pl.plot(c/self.mu[loc]/1e9/1e-6, ithro[loc]/np.median(ithro[loc]))
            pl.axis([nn(c/(filter.mu[0]*(1.0+z))/1e9/1e-6,0.1,floor=True), nn(c/(filter.mu[-1]*(1.0+z))/1e9/1e-6,0.1,ceil=True), 0, 2])
            pl.axis('tight')
            pl.xlabel('Rest wavelength (microns)')
            
        return mags


    def calcSTmag(self, filter, plot=False):
        """
        Purpose: use a filter and a spectrum to calculate the ST magnitude. 

        Inputs:
           filter  - a Spectrum class of the filter profile
           plot    - if True, plot the (interpolated) spec and (interpolated) filter for each spectrum

        Returns:
           ST magnitude for each spectrum
           
        """
        # check object type of filter: make sure Spectrum
        if isinstance(filter,object)==False: raise "Filter is not a class"

        # interpolate the spec onto the filter curve
        # find the first spec point before filter point
        startlam = np.max(self.lam[np.where(self.lam < np.min(filter.lam))[0]])
        stoplam  = np.min(self.lam[np.where(self.lam > np.max(filter.lam))[0]])

        # get the max number of points - either from number of spec
        # points or number of filter points (whichever greater)
        np_spec = float(len(np.where((self.lam > startlam) & (self.lam < stoplam))[0]))
        np_filt = float(len(filter.lam))
        npoints = max([np_spec*3,np_filt*3])

        # interpolate
        lam = np.linspace(startlam,stoplam,npoints) 
        ithro = interpolate(filter.lam, filter.flam[0], lam)

        # cycle through all spectra, calculating mag for each
        mag=[]
        for flam in self.flam:
            iflam = interpolate(self.lam, flam, lam)

            # integrate using simple trapezium method
            flux = np.sum(iflam*ithro*lam) / np.sum(ithro*lam)

            # calc mag: remember that flux in erg/s/cm**2/AA which is official unit of ST mags (no correction)
            mag.append(-2.5*np.log10(flux) - 21.10)

            if plot:
                pl.plot(lam,iflam/np.median(iflam))
                pl.plot(lam,ithro/np.median(ithro))


        return np.array(mag)
    
  
    def calcVEGAmag(self, filter, vega):
        """
        Purpose: use a filter and a spectrum to calculate the magnitude. If
                 the vega spectrum is given, return a Vega magnitude, else
                 return an AB magnitude

        Inputs:
           filter  - a Spectrum class of the filter profile
           vega    - a Spectrum class of the Vega spectrum

        Returns:
           a single magnitude wrt Vega
    
        """
        
        selfmag = self.calcABmag(filter)
        vegamag = vega.calcABmag(filter)

        mag = selfmag-vegamag
        
        return mag

    def calcM2L(self, filter, z=0.0, bandcor=False):
        """
        Purpose: use a filter, SED and details of stellar masses of sed to calc
                 the mass-to-light ratio, relative to the Sun.

                 NB: mass INCLUDES remnants as this is BC03 standard I think...

        Inputs:
           filter  - a Spectrum class of the filter profile

        Returns:
           an array of M/L for each age, wrt the Sun
           
        """

        sol = loadHSTSolarSpec()
        lsun = sol.quickABmag(filter, z=z, bandcor=bandcor)
        msun = 1.0
        
        
        l = self.quickABmag(filter, z=z, bandcor=bandcor)
        m = self.mass
        m2l = (m/10.0**(-0.4*l)) / (msun/10.0**(-0.4*lsun))

        return m2l

    def gaussLamConvolve(self, sigma_lam, nsig=5.0, verbose=True):
        """
        Purpose: to convolve the spectrum with a Gaussian IN LAMBDA SPACE of known dispersion dispersion (AA)
                 NOTE: this replicates the effect of a poorer resolution spectrograph (dlam being fixed)

        Input:
           - sigma_lam : dispersion of gaussian in Angstroms
           - nsig  : width of kernel in sigmas

        """
        # put everything on a regular grid using finest resolution in spectrum
        min_dlam = np.min(self.lam[1:]-self.lam[0:-1])

        if sigma_lam < 2.0*min_dlam:
            raise ValueError("Dispersion too small for wavelength sampling")

        # interpolate onto a regular grid
        reg_lam = np.arange(self.lam.min(), self.lam.max(),min_dlam)
        reg_spec=[]
        count=0
        for flam in self.flam:
            reg_spec.append(interpolate(self.lam,flam,reg_lam, fill_value=np.nan, bounds_error=False, method=1, kind='linear'))
            if verbose: print "Interpolated spec "+str(count)+" of "+str(len(self.flam))+" onto regular wavelength grid"
            
        # make kernel
        krange = np.ceil((sigma_lam*nsig)/min_dlam) # in units of min_dlam
        kx = np.arange(-krange,krange)
        w = kx/(sigma_lam/min_dlam)
        ky = np.exp(-0.5*w**2.0)
        # normalise
        ky/=np.sum(ky)
        
        # do convolution
        count=0
        creg_spec=[]
        for rs in reg_spec:
            count=count+1
            creg_spec.append(np.convolve(rs,ky,'same'))
            if verbose: print "Convolved spec "+str(count)+" of "+str(len(reg_spec))

        # interpolate back to original sampling
        cspec=[]
        for crs in creg_spec:
            cspec.append(interpolate(reg_lam,crs,self.lam, fill_value=np.nan, bounds_error=False, method=1, kind='linear'))
            if verbose: print "Interpolated spec "+str(count)+" of "+str(len(creg_spec))+" onto regular wavelength grid"

        self.cflam=cspec

        return cspec

    def gaussVelConvolve(self, vel, sigma, h3h4=None, instsig=None, nsig=5.0, losvd=None, overwrite=True, verbose=True):
        """
        Purpose: to convolve the spectrum with a Gaussian of known velocity (V)
                 and width (SIGMA)

        Input:
           - vel  : velocity (use 0.0 for now)
           - sigma      

        """

        # get minimum vel resolution in spec and use this for dloglam
        velscale = np.min((self.lam[1:]-self.lam[:-1])/self.lam[:-1]) * c / 1e3 # km/s
        dloglam = np.log10(1.0 + velscale/c*1e3)
        nloglam = np.round((np.log10(self.lam.max())-np.log10(self.lam.min())) / dloglam)
        self.velscale=velscale
        # calc dlam at longest wavelengths - coarsest vel resolution
        self.loglam = 10.0**np.linspace(np.log10(self.lam[0]), np.log10(self.lam[-1]), nloglam )

        count=0
        self.floglam=[]
        for flam in self.flam:
            count=count+1
            self.floglam.append(interpolate(self.lam,flam,self.loglam, fill_value=np.nan, bounds_error=False, method=1, kind='linear'))
            if verbose: print "Interpolated spec "+str(count)+" of "+str(len(self.flam))

        if losvd == None: # speed up if losvd passed
            dv = np.ceil(nsig*sigma/velscale) 
            nv = 2*dv + 1
            v = np.linspace(dv,-dv,nv) 
            w = (v - vel/velscale) / (sigma/velscale)
            w2= w*w
            if h3h4 != None:
                h3=h3h4[0]
                h4=h3h4[1]
                poly = 1.0 + h3/np.sqrt(3.0)*(w*(2.0*w2-3.0)) + \
                       h4/np.sqrt(24.0)*(w2*(4.0*w2-12.0)+3.0)
            else:
                poly = np.ones(nv)

            losvd = np.exp(-0.5*w**2.0)/(np.sqrt(2.0*np.pi)*sigma/velscale) * poly 
            losvd=losvd/np.sum(losvd)

        count=0
        self.confloglam=[]
        self.conflam = []
        for floglam in self.floglam:
            count=count+1
            confloglam = np.convolve(floglam,losvd,'same')
            self.confloglam=confloglam
            # interpolate back onto origial grid
            self.conflam.append(interpolate(self.loglam, self.confloglam, self.lam, fill_value=np.nan, bounds_error=False, method=1, kind='linear'))
            if verbose: print "Convolved spec "+str(count)+" of "+str(len(self.floglam))

        if overwrite:
            self.flam = self.conflam

    def clipSpectralRange(self, minlam, maxlam):
        """
        Clip the spectral range to be between minlam and maxlam
        """

        newflams = []
        # set range
        loc = np.where((self.lam>minlam) & (self.lam<maxlam))[0]
        # clip specs
        for flam in self.flam:
            newflams.append(flam[loc])
        if hasattr(self, 'eflam'):
            neweflams=[]
            for eflam in self.eflam:
                neweflams.append(eflam[loc])
        
        # overwite old flams
        self.lam  = self.lam[loc]
        self.flam = newflams
        if hasattr(self, 'eflam'):
            self.eflam = neweflams
        self.calcmuspec()
        
    def calcIndex(self, index, disp=None, round_prec=8, verbose=False):
        '''Function to calculate absorption line indices.
        Function takes as argument an SED of form (wavlength col, data col)
        and a choice of line index (NaI, MgI, FeH, CaT, PaT).

        Inputs:

            
            index: an indlib class from indexTools.py
            disp: float - spectral dispersion [A/pixel]
            var_SED: variance spectrum corresponding to SED of the same form as SED.

        Outputs:

            ind: value of chosen index
            ind_var_tot: variance on index measurement.
        '''


        # check uniform spacing of data in wavelegth
        data_loc_start = np.where(self.lam > np.min(index['cont_start']))[0]
        data_start = data_loc_start[data_loc_start.argmin()]#self.lam[data_loc_start].argmin()
        data_loc_stop = np.where(self.lam < np.max(index['cont_stop']))[0]
        data_stop = data_loc_stop[data_loc_stop.argmax()]
        deltas = np.round(self.lam[data_start+1:data_stop]-self.lam[data_start:data_stop-1],round_prec)
        delta = np.median(deltas)
        if type(disp)==type(None): disp=delta
        assert np.all(deltas==delta), "Wavelength spacing not uniform"
        
        # init
        indices = np.zeros(len(self.flam),dtype=float)
        index_vars = np.zeros_like(indices)
        # loop over spectra and calc indices
        for spec in xrange(len(self.flam)):
            SED = np.column_stack((self.lam, self.flam[spec]))
            if hasattr(self, 'eflam'):
                var_SED = np.column_stack((self.lam, self.eflam[spec]**2.0)) # make variance array
            else:
                var_SED = None
                
            #Method for calculating pseudo-continuum using equations from
            #Cenarro et al (2001) Appendix A.
            eps1, eps2, eps3, eps4, eps5 = 0., 0., 0., 0., 0.
            if var_SED is None:
                for j in xrange(index['ncont']):
                    #Find first and last data indices within bandpass
                    a = np.where(SED[:,0] > index['cont_start'][j])[0][0]
                    b = np.where(SED[:,0] < index['cont_stop'][j])[0][-1]
                    #Determine which pixel is closest to start and end of bandpass
                    if (index['cont_start'][j] - SED[a-1,0]) < (SED[a,0] - index['cont_start'][j]):
                        a -= 1
                    if (index['cont_stop'][j] - SED[b,0]) > (SED[b+1,0] - index['cont_stop'][j]):
                        b += 1      

                    vals = np.copy(SED[a:b+1,:])

                    e1 = np.ones(len(vals[:,0]))
                    eps1 += np.sum(e1)            
                    e2 = vals[:,0]
                    eps2 += np.sum(e2)           
                    e3 = vals[:,0]**2
                    eps3 += np.sum(e3)            
                    e4 = vals[:,1]
                    eps4 += np.sum(e4)            
                    e5 = vals[:,0]*vals[:,1]
                    eps5 += np.sum(e5)

                delta = eps1*eps3 - eps2**2
                alpha1 = (1/float(delta))*(eps3*eps4 - eps2*eps5)
                alpha2 = (1/float(delta))*(eps1*eps5 - eps2*eps4)

            #If real spectrum with variance spectrum:
            elif var_SED is not None:
                for j in xrange(index['ncont']):
                    #Find first and last data indices within bandpass
                    a = np.where(SED[:,0] > index['cont_start'][j])[0][0]
                    b = np.where(SED[:,0] < index['cont_stop'][j])[0][-1]
                    #Determine which pixel is closest to start and end of bandpass
                    if (index['cont_start'][j] - SED[a-1,0]) < (SED[a,0] - index['cont_start'][j]):
                        a -= 1
                    if (index['cont_stop'][j] - SED[b,0]) > (SED[b+1,0] - index['cont_stop'][j]):
                        b += 1      

                    vals = np.copy(SED[a:b+1,:])
                    err_vals = np.copy(var_SED[a:b+1,:])

                    e1 = np.ones(len(vals[:,0]))/err_vals[:,1]
                    eps1 += np.sum(e1)            
                    e2 = vals[:,0]/err_vals[:,1]
                    eps2 += np.sum(e2)            
                    e3 = vals[:,0]**2/err_vals[:,1]
                    eps3 += np.sum(e3)            
                    e4 = vals[:,1]/err_vals[:,1]
                    eps4 += np.sum(e4)            
                    e5 = vals[:,0]*vals[:,1]/err_vals[:,1]
                    eps5 += np.sum(e5)

                delta = eps1*eps3 - eps2**2
                alpha1 = (1/float(delta))*(eps3*eps4 - eps2*eps5)
                alpha2 = (1/float(delta))*(eps1*eps5 - eps2*eps4)
            
##            #Test to check that continua works
##            xs = np.arange(cont[0], cont[-1], 1)
##            ys = alpha2*xs + alpha1
##            pl.plot(SED[:,0], SED[:,1], 'k-', xs, ys, 'r--')
##            pl.show()
                
            #Next stage, calculate index
            ind = 0.; ind_var_tot = 0.
            ind_var_1 = 0.; ind_var_2 = 0.
            for j in xrange(index['nfeat']):
                #Find first and last data indices within bandpass
                a = np.where(SED[:,0] > index['ind_start'][j])[0][0]
                b = np.where(SED[:,0] < index['ind_stop'][j])[0][-1]
                #Determine which pixel is closest to start and end of bandpass
                if (index['ind_start'][j] - SED[a-1,0]) < (SED[a,0] - index['ind_start'][j]):
                    a -= 1
                if (index['ind_stop'][j] - SED[b,0]) > (SED[b+1,0] - index['ind_stop'][j]):
                    b += 1
                #Multiplicative factors for start and end pixels
                Cstart_c = (SED[a,0] - index['ind_start'][j] + 0.5*disp)/disp
                Cend_c = (index['ind_stop'][j] - SED[b,0] + 0.5*disp)/disp

                Cvals = SED[a:b+1,:]
                
                contys = alpha2*Cvals[:,0] + alpha1
                Ccontvals = np.column_stack((Cvals[:,0],contys))

                array = (1-Cvals[:,1]/Ccontvals[:,1])
                array[0] *= Cstart_c
                array[-1] *= Cend_c
                value = disp*np.sum(array)
                ind += value

                ###ERRORS:
                if var_SED is not None:
                    #Index error:
                    Cerrvals = var_SED[a:b+1,:]
                    
                    ind_cont_var = np.zeros(len(Ccontvals),dtype=float)
                    #Loop over continuum feature pixels
                    for x in xrange(len(ind_cont_var)): 
                        #Calculate continuum variance:
                        Cvar_cont = 0.
                        for i in xrange(index['ncont']):
                            #Find first and last data indices within bandpass
                            ia = np.where(SED[:,0] > index['cont_start'][i])[0][0]
                            ib = np.where(SED[:,0] < index['cont_stop'][i])[0][-1]
                            #Determine which pixel is closest to start and end of bandpass
                            if (index['cont_start'][i] - SED[ia-1,0]) < (SED[ia,0] - index['cont_start'][i]):
                                ia -= 1
                            if (index['cont_stop'][i] - SED[ib,0]) > (SED[ib+1,0] - index['cont_stop'][i]):
                                ib += 1      
                            #Extract relevent region of SED and variance SED
                            contvals = np.copy(SED[ia:ib+1,:])
                            var_contvals = np.copy(var_SED[ia:ib+1,:])
                            #Apply equation
                            contsum = ((1/delta)*((1/var_contvals[:,1])*eps3 - (contvals[:,0]/var_contvals[:,1])*eps2)\
                                      + (Ccontvals[x,0]/delta)*((contvals[:,0]/var_contvals[:,1])*eps1\
                                                                - (1/var_contvals[:,1])*eps2))**2.*var_contvals[:,1]
                            Cvar_cont += np.sum(contsum)
                        ind_cont_var[x] = Cvar_cont

                    Carray_var_1 = (Ccontvals[:,1]**2*Cerrvals[:,1] + Cvals[:,1]**2*ind_cont_var)/\
                                 Ccontvals[:,1]**4.
                    Carray_var_1[0] *= Cstart_c
                    Carray_var_1[-1] *= Cend_c
                    ind_var_1 += np.sum(Carray_var_1)

                    #For part 2 of equation which includes covariance matrix:
                    #Loop over pixels in spectral feature
                    for y in xrange(len(Cvals[:,0])):
                        #Loop over first spectral feature
                        for z in xrange(index['nfeat']):
                            #Find first and last data indices within bandpass
                            iia = np.where(SED[:,0] > index['ind_start'][z])[0][0]
                            iib = np.where(SED[:,0] < index['ind_stop'][z])[0][-1]
                            #Determine which pixel is closest to start and end of bandpass
                            if (index['ind_start'][z] - SED[iia-1,0]) < (SED[iia,0] - index['ind_start'][z]):
                                iia -= 1
                            if (index['ind_stop'][z] - SED[iib,0]) > (SED[iib+1,0] - index['ind_stop'][z]):
                                iib += 1
                            #Multiplicative factors for start and end pixels
                            Fstart_c = (SED[iia,0] - index['ind_start'][z] + 0.5*disp)/float(disp)
                            Fend_c = (index['ind_stop'][z] - SED[iib,0] + 0.5*disp)/float(disp)
                            
                            Fvals = SED[iia:iib+1,:]
                            Fcontys = alpha2*Fvals[:,0] + alpha1
                            Fcontvals = np.column_stack((Fvals[:,0],Fcontys))
                            #Loop over pixels in second spectral feature
                            for zz in xrange(len(Fvals[:,0])):
                                cov_part = ((1/delta**2)*(eps1*eps3*eps3-eps2*eps2*eps3))+\
                                           ((1/delta**2)*(eps2*eps2*eps2-eps1*eps2*eps3))*(Cvals[y,0] + Fvals[zz,0])+\
                                           ((1/delta**2)*(eps1*eps1*eps3-eps1*eps2*eps2))*(Cvals[y,0]*Fvals[zz,0])
                                
                                part_2_num = (Cvals[y,1]*Fvals[zz,1]*cov_part)/(Ccontvals[y,1]**2*Fcontvals[zz,1]**2)
                                #edge pixel factors
                                if y == 0.:
                                    part_2_num *= Cstart_c
                                if zz == 0.:
                                    part_2_num *= Fstart_c
                                if y == (len(Cvals[:,0])-1):
                                    part_2_num *= Cend_c
                                if zz == (len(Fvals[:,0])-1):
                                    part_2_num *= Fend_c

                                ind_var_2 += part_2_num

                    ind_var_tot = disp**2*(ind_var_1 + ind_var_2)

            if verbose: print index['name']+' = %.3f pm %.3f Angstroms' % (ind, np.sqrt(ind_var_tot))
            indices[spec] = ind
            index_vars[spec] = ind_var_tot

        rlist = indices
        if var_SED is not None:
            rlist = [indices, np.sqrt(index_vars)]
        
        return rlist #indices#, index_vars


    def CaT_star(self, disp, var_SED=None):
        '''Function to calculate CaT* index

        Inputs:

            SED - spectrum
            disp: spectral dispersion [A/pixel]
            var_SED: variance spectrum

        Outputs:

            CaT_star: CaT* index value
            CaT_star_var: index variance
        '''

        cat = self.irindex(disp, index='CaT', var_SED=None)
        pat = self.irindex(disp, index='PaT', var_SED=None)

        cat_star = cat - 0.93*pat
 #       cat_star_var = cat[1] + 0.93**2*pat[1]

        return cat_star


    def normaliseSpec(self, polyOrder=3, indLib=None, maxErrVal=999.0, overwrite=True):
        """
        Normalise spectrum continuum to one, avoiding abs lines

        """

        # scale wave to be between -1 and 1
        x = 2.0*(self.lam-self.lam.min())/(self.lam.max()-self.lam.min()) - 1.0

        nflam  = []
        pfits  = []
        
        # sort out err 
        if hasattr(self, 'eflam'):
            neflam = []
            eflams = self.eflam
        else:
            eflams = [None for flam in self.flam]
        # loop over spectra, normalising
        for flam, eflam in zip(self.flam, eflams):
            # init
            xfit = np.copy(x)
            wfit = np.copy(self.lam)
            yfit = np.copy(flam)
            if eflam is not None: efit = np.copy(eflam)

            # remove nans and err<0.0
            if eflam is not None:
                loc=np.where((~np.isfinite(flam)) | (~np.isfinite(eflam)) | (eflam<=0.0))
            else:
                loc=np.where(~np.isfinite(flam))  
            good = np.ones_like(self.lam, dtype=np.bool)
            good[loc]=False
            keep = np.where(good)
            xfit = xfit[keep]; wfit = wfit[keep]; yfit = yfit[keep];
            if eflam is not None: efit=efit[keep]

            # remove abs lines from poly fit
            if type(indLib)!=type(None):
                locs=[]
                for ind in indLib.names:
                    loc = np.where((wfit>getattr(indLib,ind)['ind_start']) & (wfit<getattr(indLib,ind)['ind_stop']))[0]
                    if len(loc)>0:
                        locs.extend(loc)
                good = np.ones_like(xfit, dtype=np.bool)
                good[locs] = False
                keep = np.where(good)
                xfit = xfit[keep]
                wfit = wfit[keep]
                yfit = yfit[keep]
                if eflam is not None: efit=efit[keep]

            # fit poly
            if eflam is not None:
                coef = cheby.chebfit(xfit, yfit, deg=polyOrder, w=1.0/efit**2.0) # variance weight poly fit
            else:
                coef = cheby.chebfit(xfit, yfit, deg=polyOrder)

            polyFit = cheby.chebval(x,coef)

            # normalise 
            nflam.append(flam/polyFit)
            pfits.append(polyFit)

            # return extra arrays
            if eflam is not None:
                neflam.append(eflam/polyFit)
                
        self.nflam = np.array(nflam)
        self.pfits = np.array(pfits)
        if hasattr(self, 'eflam'): self.enflam = np.array(neflam)
        
        if overwrite:
            self.flam = self.nflam
            if hasattr(self, 'eflam'): self.eflam = self.neflam
    
        

###################################### END OF SPECTRUM CLASS #####################################


def linterp(x,y,x0):
    """
    Author: Ryan Houghton (18/4/11)

    Purpose: Linearly interpolate value at x0 from values y at x

    Inputs:
       x - array of x coords
       y - array of y coords
       x0- x coord at which to calculate value (NOT an array)

    Returns:
       y value corresponding to x0
    """

    if len([x0]) > 1: raise "Arrays not allowed for x0"
    
    a, b = 0, len(x)-1
    sgn=1
    if(x[0]>x[-1]): sgn=-1
    while (b-a)>1:
        m=(a+b)/2
        if sgn*x[m]>=sgn*x0: b=m
        else: a=m
    return y[a]+(y[a+1]-y[a])/(x[a+1]-x[a])*(x0-x[a])


def alinterp(x,y,newx):
    """
    Author: Ryan Houghton (18/4/11)

    Purpose: Linearly interpolate values at xnew from values y at x

    Inputs:
       x - array of x coords
       y - array of y coords
       xnew - x coords for which values to be calculated (array allowed)

    Returns
       y values cooresponding to xnew
    """
    return np.array( [linterp(x,y,xval) for xval in newx] )


def interpolate(x,y,newx, fill_value=0.0, bounds_error=False, method=1, kind='linear'):
    """
    Just a wrapper for my favourite choice of interpolation
    """

    if method==1:
        newy = alinterp(x,y,newx) # this does a poor job when 
    elif method==2:
        loc = x.argsort()
        x.sort()
        y=y[loc]
        f = ip.interp1d(x, y, fill_value=fill_value, bounds_error=bounds_error, kind=kind)
        newy = f(newx)
    else:
        raise "Method not understood"
    
    return newy

def findfiles(path, glob='', vglob=None):
    """
    Purpose: returns all the files that match the *path* (REGEXP)

    - glob  : regexp for what you want
    - vglob : string for what you don't want

    e.g to find all files beginning with "IFU" but rejecting files with "_s" in:

    > from stellpops import specTools as st
    > files, nfiles = st.findfiles("./", glob="IFU*", vglob="_s")
    """
    
    # get all files
    if vglob is None:
        (si,so,se) = os.popen3("ls "+path+glob)
    else:
        (si,so,se) = os.popen3("ls "+path+glob+" | grep -v "+vglob)
    si.close()
    se.close()
    files = so.readlines()
    so.close()

    # strip the \n from the end
    fnames=[]
    for name in files: fnames.append(name.split()[0])

    # total number found
    nfiles = len(fnames)

    return fnames, nfiles

def loadMyFilters(dir="~/z/data/stellar_pops/myfilters/", file="filters.dat", \
                  verbose=False):

    """
    Author: Ryan houghton (16/4/11)

    Purpose: Read  my filters into a dictionary

    Inputs:
       dir  - where the filter file is kept
       file - the name of the file
       verbose - if True, print info about process

    Returns:
       filters - a dictionary of the filters in the file. Each key is linked to
                 a 2D array of lambda (AA) and transmission


    Examples:

       > from specTools import loadFilters
       > filters = loadMyFilters()
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
    nfilt = strg.atof(txt[0].split()[0])

    # read in the filters
    count=1
    while count < nlines:

        # get filter name and number of lines in filter definition
        filt = txt[count].split("'")[1]
        nlam = strg.atoi(txt[count].split()[0])

        # read in the wavelength and throughput info into Spectrum()
        lamspec=np.zeros((nlam))
        lam    =np.zeros((nlam)) 
        for n in range(count+1,count+nlam+1):

            line = txt[n].split()
            l = strg.atof(line[0])
            f = strg.atof(line[1])
            lam[n-count-1]=l
            lamspec[n-count-1]=f
            

        s=spectrum(lam=lam, lamspec=lamspec, filter=True)
        filters[filt]=s
        count+=nlam+1
        if verbose: print "Loaded "+filt

    return filters

    
def loadGMOSFilters(dir="~/z/data/stellar_pops/myfilters/gmos", \
                    glob="gmos*.txt"):
    """
    Purpose: load up the GEMINI/GMOS filter set (based on SDSS but
             slightly different, see Jorgensen 2009)

    Inputs:
       dir  - the full path to the directory containing the filters
       glob - the REGEXP to select the desired files. 

    Example:
       > from stellpops import specTools as st
       > filts = st.loadGMOSFilters()

       
    """


    # get files

    files = findfiles(expanduser(dir)+"/"+glob)[0]

    # read in filters/files
    filts={}
    for fname in files:

        lam, flux = np.loadtxt(fname, unpack=True, skiprows=16)
        nlam = len(lam)
        # get filter name from filename
        elems = (fname.split("/")[-1]).split("_")
        filt = elems[2]+"_"+elems[3].split(".")[0]

        # make dict of spectrum classes: x10.0 to convert nm to AA
        filts[filt]=spectrum(lam=lam*10.0, lamspec=flux, filter=True)

    return filts


def loadACSFilters(dir="~/z/data/stellar_pops/myfilters/hst/acs/", \
                    glob="*.dat"):
    """
    Purpose: load up the HST/ACS/WFPC2 filter set 

    Inputs:
       dir  - the full path to the directory containing the filters
       glob - the REGEXP to select the desired files. 

    Example:
       > from stellpops import specTools as st
       > filts = st.loadACSFilters()

       
    """


    # get files

    files = findfiles(expanduser(dir)+"/"+glob)[0]

    # read in filters/files
    filts={}
    for fname in files:

        lam, flux = np.loadtxt(fname, unpack=True, skiprows=16)
        nlam = len(lam)
        # get filter name from filename
        elems = (fname.split("/")[-1]).split("_")
        filt = elems[1]

        # make dict of spectrum classes
        filts[filt]=spectrum(lam=lam, lamspec=flux, filter=True)

    return filts

def load2MASSFilters(dir="~/z/data/stellar_pops/myfilters/2mass/"):
    """
    RH 24/7/14

    Load up the official 2mass filter functions
    """
    Jfile="J_rsr.dat"
    Hfile="H_rsr.dat"
    Kfile="Ks_rsr.dat"
    specs={}
    
    # open files for reading, read, and then close file
    f = open(expanduser(dir)+Jfile,'r')
    txt = f.readlines()
    f.close()
    lam,flux=[],[]
    for t in txt:
        st=t.split()
        lam.append(st[0])
        flux.append(st[1])
    lam = np.array(lam,dtype=np.float)*1e4
    flux=np.array(flux,dtype=np.float)
    flux/=np.sum(flux)
    specs['J_2mass'] = spectrum(lam=lam, lamspec=flux, filter=True)
    
    f = open(expanduser(dir)+Hfile,'r')
    txt = f.readlines()
    f.close()
    lam,flux=[],[]
    for t in txt:
        st=t.split()
        lam.append(st[0])
        flux.append(st[1])
    lam = np.array(lam,dtype=np.float)*1e4
    flux=np.array(flux,dtype=np.float)
    flux/=np.sum(flux)
    specs['H_2mass'] = spectrum(lam=lam, lamspec=flux, filter=True)

    f = open(expanduser(dir)+Kfile,'r')
    txt = f.readlines()
    f.close()
    lam,flux=[],[]
    for t in txt:
        st=t.split()
        lam.append(st[0])
        flux.append(st[1])
    lam = np.array(lam,dtype=np.float)*1e4
    flux=np.array(flux,dtype=np.float)
    flux/=np.sum(flux)
    specs['Ks_2mass'] = spectrum(lam=lam, lamspec=flux, filter=True)

    return specs

def loadHSTFilters(dir="~/z/data/stellar_pops/myfilters/hst/", file="filter_HST.res"):
    """

    RH 3/10/13

    Purpose: read in the library of HST filter curves from the MILES website (for Simon)

    Inputs:
       dir  - where the txt file is located
       file - the name of the text file
    
    """
    # open file for reading, read, and then close file
    f = open(expanduser(dir)+file,'r')
    txt = f.readlines()
    f.close()

    # find the filter breaks and add filters to dict
    nlines=len(txt)
    filters={}
    count=0

    fname=strg.lower(txt[0].split()[1]) # use lowercase to be consistent with previous functions (loadACSFilters)
    filt=[]
    for t in txt[1:]:
        if t[0]=='#': # if new filter starting
            fdat=np.array(filt).T
            filters[fname]=spectrum(lam=fdat[0], lamspec=fdat[1], filter=True)
            filt=[]
            fname=strg.lower(t.split()[1])
            print "starting read of "+fname
        else:
            dat=t.split()
            filt.append([strg.atof(dat[0]),strg.atof(dat[1])])
            
        count=count+1
    
    fdat=np.array(filt).T
    filters[fname]=spectrum(lam=fdat[0], lamspec=fdat[1], filter=True)
    print "Loaded "+str(len(filters))+" from "+str(count)+" lines of text"
    return filters


def loadHSTFiltersSTSCI(dir="~/Data/stelpopmods/Filters/hst/", glob="*.txt"):
    """

    RH & SZ 31/10/13

    Purpose: read in the library of HST filter curves from the STSci website

    Inputs:
       dir  - where the txt file is located
       file - the name of the text file
    
    """
    # get files

    files = findfiles(expanduser(dir)+"/"+glob)[0]

    print files

    # read in filters/files
    filts={}
    for fname in files:

        lam, flux = np.loadtxt(fname, unpack=True)
        nlam = len(lam)
        # get filter name from filename
        #elems = (fname.split("/")[-1]).split("_")
        elems = (fname.split("/")[-1])
        #filt = elems[1]
        filt = elems[0:5]

        # make dict of spectrum classes
        filts[filt]=spectrum(lam=lam, lamspec=flux, filter=True)

    return filts
    

def calcKcorAB(spec, redshift, filter1, filter2=None, quick=False, bandcor=True):
    """
    Purpose: calculate the K-correction (AB mags) given a spectrum, a redshift and a filter, filter1.
             Optionally, calculate the K-correction relative to a different filter, filter2.

             Returns K such that mag_0 = mag_z + K
             
    """

    # some quick checks
    #if not(isinstance(spec, spectrum)): raise "Spec is not a Spectrum"
    #if not(isinstance(filter1, spectrum)): raise "Filter1 is not a Spectrum"
    #if (not(isinstance(spec, spectrum)) &  (not(filter2==None))): raise "Filter2 is not a Spectrum, nor is it None type"

    if filter2 is not None:
        filter3=filter2
    else:
        filter3=filter1
    
    if quick:
        mag_0 = spec.quickABmag(filter3,z=0.0)
        mag_z = spec.quickABmag(filter1,z=redshift, bandcor=bandcor)
    else:
        mag_0 = spec.calcABmag(filter3,z=0.0)
        mag_z = spec.calcABmag(filter1,z=redshift, bandcor=bandcor)

    return mag_z - mag_0 # return the standard colour (blue-red)


def loadHSTSolarSpec(dir="~/z/data/stellar_pops/HSTcals/", file="sun_reference_stis_001.ascii"):
    """
    Author: Ryan Houghton 24/5/2011
    Purpose: Read in the HST calibration solar spectrum (http://www.stsci.edu/hst/observatory/cdbs/calspec.html)

    Example
    -------

    > import specTools as st
    > solar = st.loadHSTSolarSpec()

    """

    

    lam, flux = np.loadtxt(expanduser(dir)+"/"+file, unpack=True)

    AbsFlux = flux * (d_sun/10.0/pc)**2.0

    solar = spectrum(lamspec=AbsFlux, lam=lam, model="HSTcal", Z=0.02)

    return solar


def loadHSTVegaSpec(dir="~/z/data/stellar_pops/HSTcals/", file="alpha_lyr_stis_005.ascii"):
    """
    Author: Ryan Houghton 24/5/2011
    Purpose: Read in the HST calibration Vega spectrum (http://www.stsci.edu/hst/observatory/cdbs/calspec.html)

    Example
    -------

    > import specTools as st
    > vega = st.loadHSTSolarSpec()

    """

    lam, flux = np.loadtxt(expanduser(dir)+"/"+file, unpack=True)

    AbsFlux = flux * (7.76/10.0)**2.0

    vega = spectrum(lamspec=flux, lam=lam, model="HSTcal", Z=0.02)

    return vega

def loadSpecs(verbose=False):
    """
    Quickly load the BC03 and M05 specs
    """

    bcs = bc.loadBC03ssps(verbose=verbose)
    mas = ma.loadM05ssps(verbose=verbose)
    specs = []
    for s in bcs: specs.append(s)
    for s in mas: specs.append(s)

    return specs

def rebin(x, y, newx, verbose=True, flux=True):
    """
    RH 23/4/15

    Taken from fitStar.py

    Rebin a spectrum conserving flux (unless Flux==False)
    
    """
    
    # get length of inputs
    nold = len(x)
    nnew = len(newx)
    
    # get element sizes
    doldx = np.fabs(np.median(x-np.roll(x,1)))
    dnewx = np.fabs(np.median(newx-np.roll(newx,1)))

    # get array limits and edges (in units of element sizes)
    oldlimits = np.array([x.min(), x.max()]) / doldx + np.array([-0.5,0.5])
    oldedges  = np.linspace(oldlimits[0], oldlimits[1],nold+1)

    newlimits = np.array([newx.min(), newx.max()]) / doldx + np.array([-0.5,0.5])*dnewx/doldx
    newedges  = np.linspace(newlimits[0], newlimits[1],nnew+1)

    # make an array in index new bins to old bins
    binindex = np.floor(newedges-oldlimits[0])
    binindex[np.where(binindex<=0)]=0
    binindex[np.where(binindex>=(nold-1))]=nold-1

    # calc the new spec
    newy = np.zeros(nnew)
    for ni in range(nnew):
        newy[ni] = np.sum(y[binindex[ni]:binindex[ni+1]+1]) \
                   - (newedges[ni] - oldedges[binindex[ni]]) * y[binindex[ni]] \
                   - (oldedges[binindex[ni+1]+1] - newedges[ni+1]) * y[binindex[ni+1]]


    # undo flux conservation (becomes NEAREST NEIGHBOUR interpolation)
    if flux==False:
        newy /= newedges[1:] - newedges[0:-1]

    return newy

def vac2air(wave_vac, useIDLcoef=True, verbose=False):
    """
    RH 14/3/15

    Code adapted from NASA IDL version, but with alt method/polynomial too.

    Input/Output in AA
    
    """
    wave_vac = np.atleast_1d(wave_vac)
    wave_air = np.copy(np.array(wave_vac))
    g = np.where(wave_vac> 2000)[0] # Only modify above 2000 A
    ng = len(g)
    if ng>0:
        sigma2 = (1e4/wave_vac[g] )**2.0   # Convert to wavenumber squared
        # Compute conversion factor
        if useIDLcoef:
            fact = 1.0 +  5.792105e-2/(238.0185e0 - sigma2) + 1.67917e-3/( 57.362e0 - sigma2)
        else:
            fact = 1.0 + 6.4328e-5 + 2.94981e-2/(146.0 - sigma2) + 2.5540e-4/( 41.0 - sigma2)
        wave_air[g] = wave_vac[g]/fact  # Convert wavelengths
    return wave_air
    
def air2vac(wave_air, useIDLcoef=True, verbose=False):
    """
    RH 14/3/15

    Code adapted from NASA IDL version, but with alt method/polynomial too.

    Input/Output in AA
    """
    wave_air = np.atleast_1d(wave_air)
    wave_vac = np.copy(wave_air)
    g = np.where(wave_air>2000)[0]          # Only modify above 2000 A
    ng= len(g)
    if ng>0:
        sigma2 = (1e4/wave_air[g])**2.0     # Convert to wavenumber squared
        # Compute conversion factor
        if useIDLcoef:
            fact = 1.0 + 5.792105e-2/(238.0185e0 - sigma2) + 1.67917e-3/( 57.362e0 - sigma2)
        else:
            fact = 1.0 + 6.4328e-5 + 2.94981e-2/(146.0 - sigma2) + 2.5540e-4/( 41.0 - sigma2)
        wave_vac[g] = wave_air[g]*fact      # Convert Wavelength
    return wave_vac

