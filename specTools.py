import numpy as np
import pylab as pl
import string as strg
import scipy.interpolate as ip
from numpy.polynomial import chebyshev as cheby
#from stellpops import BC03tools as bc
#from stellpops import M05tools as ma
import re
import os
import copy as cp
import warnings as warn
import pdb
from os.path import expanduser
from nearest import nearest as nn

# globals
c = 299792458.0 # m/s
pc= 3.08568025E16
d_sun = 149597887.5e3 # m
sigmaInFWHM = np.sqrt(8.*np.log(2.0))
L_sun = 3.826E33 # the L_sun defined by BC03 in erg/s

class spectrum:
    """
    Author: Ryan Houghton (20/4/11)

    Purpose: A spectrum object to aid manipulation of wavelength/frequency and the associated
             flux values

             Originally developed to calculate magnitudes of various SPS models at different ages/metallicities.

             Now adapted to also calculate spectral indices over specific absorption lines. 

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
                 age=None, mass=None, alpha=None, Z=None, IMF=None, filter=False, model=None, resolution=None, wavesyst=None, \
                 userdict=None, debug=False):

        """
        Inputs:
        -------
           lam     - the 1D array of wavelengths (becomes self.lam)
           lamspec - the (1D, 2D or 3D) array of fluxes on the specified lambda grid (becomes self.flam)
           errlamspec - standard deviation for pixels in lamspec (becomes self.eflam)
           mu      - the 1D array of frequencies (becomes self.mu)
           muspec  - the (1D, 2D or 3D) array of fluxes on the specified mu grid (becomes self.fmu)
           errmuspec - standard deviation for pixels in muspec (becomes self.efmu)
           age     - the 1D array of ages (Gyr) for the spectra
           mass    - the 1D array of stellar masses (M_Sun) at each age for each spectrum
           alpha   - the 1D array of [alpha/Fe] values corresponding to each spectrum
           Z       - the 1D array of metallicities for the spectra
           IMF     - string describing IMF type
           filter  - if the spectrum is a filter transmission (photon %), then set this to TRUE. This affects how we calc Fmu
           model   - model name (e.g. M05, BC03, CD12, ATRAN)
           resolution - 2 element list of
               [0]: spectral R=lam/dlam (FWHM),
               [1]: resolution FWHM (val/tuple/dict) in AA
                    - may be single value for the FWHM in AA
                    - may be list/tuple defining starting and final resolution FWHM in AA.
                    - may be a dict for 'visible' and 'nir' wavelengths if multiple libs used; the dict value may be single (FWHM),
                      or tuple as above
               Either may be None, depending on how the models/libs are defined.
               
               SEE calcResolution BELOW FOR MORE INFO!!!
           wavesyst - you MUST specify if the wavelengths are in the AIR or VAC system. 
           userdict- add extra user-defined tags to the spectrum, such as airmass, water vapour etc.

        Notes:
        ------
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
        
        # set which attributes may be arrays, excluding flam and fmu (which obviously can be)
        self.arrayAttrs=['age', 'Z', 'alpha', 'mass', 'IMF']
        self.unrastered=False # Spectra always loaded in an unrasterised, multi-D manner
        
        # sanity check: check if initalised "correctly"
        assert (lamspec is not None) or (muspec is not None), \
               "Spectrum must be defined with either lamspec(AA) or muspec (Hz)"
        assert not ((lamspec is not None) and (muspec is not None)), \
               "Please don't define BOTH lamspec and muspec: just one or the other"

        ##########################################################
        # Define spectrum part that was given
        ##########################################################
        if lamspec is not None:
            self._processLamspec(lam, lamspec, errlamspec=errlamspec)
        elif muspec is not None:
            self._processMuspec(mu, muspec, errmuspec=errmuspec)

        if debug: pdb.set_trace()
        #############################################################
        # Define spectrum attributes
        #############################################################
        
        # add age info
        if age is not None:
            # scalar or array
            self.age = singleOrList2Array(age)
            # check dims
            self.checkDims("age")
            # add useful logage
            self.logage = np.log10(self.age)
        else:
            self.age = None

        # add stellar mass info, in same way added age info
        if mass is not None:
            self.mass = singleOrList2Array(mass)
            self.checkDims("mass")
            self.logmass = np.log10(self.mass)
        else:
            self.mass = None

        # add alpha info in similar way as added age
        if alpha is not None:
            self.alpha = singleOrList2Array(alpha)
            self.checkDims("alpha")
        else:
            self.alpha = None

        # add metallicitiy in same way to alpha
        if Z is not None:
            self.Z = singleOrList2Array(Z)
            self.checkDims("Z")
        else:
            self.Z = None

        # add IMF in same way 
        if IMF is not None:
            self.IMF = singleOrList2Array(IMF)
            self.checkDims("IMF")
        else:
            self.IMF = None

        # name of the model, e.g. BC03: allow only scalar values here
        if model is not None:
            assert np.isscalar(model), "Model not scalar"
            self.model = model 
        else:
            self.model = None

        # add the resolution in AA, allow only 2-elem-list here
        if resolution is not None:
            assert len(resolution)==2, "Resolution not understood"
            self.resolution = resolution 
        else:
            self.resolution = None

        # add VAC or AIR - specific strings allowed, or give warning
        if wavesyst is not None:
            if (wavesyst=="vac" or wavesyst=="VAC" or wavesyst=="Vac"):
                self.wavesyst="vac"
            elif (wavesyst=="air" or wavesyst=="AIR" or wavesyst=="Air"):
                self.wavesyst="air"
            else:
                raise ValueError("wavesyst not understood. Should be air or vac.")
        else:
            warn.warn("You failed to specify if the wavelength is defined in AIR or VAC units.")
            self.wavesyst=None

        # add user dictionary for extra info
        if userdict is not None:
            self.__userdict__ = userdict
            keys = userdict.keys()
            for key in keys:
                # add info as scalar or array
                setattr(self, key, singleOrList2Array(userdict[key]))
                # make sure the dims agree
                self.checkDims(key)
        else:
            self.__userdict__ = None

        # apply the filter keyword
        if filter:
            self.__filter__=filter
        else:
            self.__filter__=False


        ###############################################################################################################
        # Define the spectrum part that wasn't given - this needs to be here due to the way unraster and reraster work
        ###############################################################################################################
        if debug: pdb.set_trace()
        if hasattr(self, "flam"):
            # calc the corresponding muspec
            self.calcmuspec(filter=filter, debug=debug)
        elif hasattr(self,"fmu"):
            # calc the corresponding lamspec
            self.calclamspec(filter=filter)
        if debug: pdb.set_trace()

    def _processLamspec(self, lam, lamspec, errlamspec=None):
        """
        To be usually called by __init__: process the lam and flam arrays 
        """
        # check that lam has been given
        if lam is None: raise "If you give lamspec, you must also give lam"
        
        # make sure 2d
        flam = np.atleast_2d(singleOrList2Array(lamspec))
        # get array size
        loc = flam.shape

        # assign dims of non-spectral component
        self.dims = list(loc[:-1])
        self.nspec = self.getSize()
        # assign dims of spectral component
        self.nlam = loc[-1]

        # assign data to class
        self.lam  = np.array(lam)
        self.flam = flam

        # deal with error spectra
        if errlamspec is not None:
            eflam = np.atleast_2d(errlamspec)
            eloc = eflam.shape
            self.eflam = eflam
            # sanity check
            assert np.all(loc==eloc), "Flux and error arrays appear different sizes..."
        else:
            self.eflam = None

    def _processMuspec(self, mu, muspec, errmuspec=None):
        """

        """
        # check that mu has been given
        if mu is None: raise "If you give muspec, you must also give mu"
        # make sure 2D
        fmu = np.atleast_2d(muspec)
        # get array size
        loc = fmu.shape

        # assign dims of non-spectral component
        self.dims = list(loc[:-1])
        self.nspec = self.getSize()
        # assign dims of spectral component
        self.nmu = loc[-1]

        # assign data to class
        self.mu  = np.array(mu)
        self.fmu = np.atleast_2d(singleOrList2Array(fmu))

        # deal with error spectra
        if errmuspec is not None:
            efmu = np.atleast_2d(errmuspec)
            eloc = efmu.shape
            self.efmu = efmu
            # sanity check
            assert np.all(loc==eloc), "Flux and error arrays appear different sizes..."
        else:
            self.efmu=None
            
        
    def checkDims(self, attr):
        """
        Check the dimensions of the given attribute match the dimensions of flam and fmu
        """
        checkDims(getattr(self, attr), attr, self.dims)
        
    # define a copy function
    def copy(self, loc=None):
        """
        Make an copy of the spectrum, so that if you change the copy, it DOESNT change the original.
        """
        if loc is None:
            # simple copy of everything
            dself = cp.deepcopy(self.__dict__)
            newspec = spectrum(lam=dself['lam'], lamspec=dself['flam'], errlamspec=dself['eflam'], \
                               age=dself['age'], mass=dself['mass'], alpha=dself['alpha'], Z=dself['Z'], \
                               IMF=dself['IMF'], filter=dself['__filter__'], model=dself['model'], \
                               resolution=dself['resolution'], wavesyst=dself['wavesyst'], \
                               userdict=dself['__userdict__'])
        else:

            # simple copy of everything
            dself = cp.deepcopy(self.__dict__)

            coreLocArgs = ['flam', 'eflam', 'mass', 'age', 'alpha', 'Z', 'IMF']
            # now check indiv args
            for cla in coreLocArgs:
                # for each core arg with loc indexing, try and index it. If you fail, stick with original
                if dself[cla] is not None:
                    try:
                        dself[cla] = dself[cla][loc]
                    except:
                        pass
            # check userdict
            if dself['__userdict__'] is not None:
                keys = dself['__userdict__'].keys()
                for k in keys:
                    # for each key, try and address it with loc. Don't fret if not possible.
                    try:
                        dself['__userdict__'][k] = singleOrList2Array(dself['__userdict__'][k])[loc]
                    except:
                        pass
                    
            newspec = spectrum(lam=dself['lam'], lamspec=np.atleast_2d(dself['flam']), errlamspec=dself['eflam'], \
                               age=dself['age'], mass=dself['mass'], alpha=dself['alpha'], Z=dself['Z'], \
                               IMF=dself['IMF'], filter=dself['__filter__'], model=dself['model'], \
                               resolution=dself['resolution'], wavesyst=dself['wavesyst'], \
                               userdict=dself['__userdict__'])
            
        return newspec


    def getOrderOfDims(self):
        """
        For a multiD format, the spectra are ordered in a multi dimensional array.
        This functions returns which parameter is associated with each dimension, e.g. age, alpha, Z, etc.

        """

        # calc uniq values 
        uvals = self.calcUniqAttrs(augment=False)
        # get dims
        loc = np.array(self.dims)
        # loop over dims, finding which one has same length as each array attrs
        dims = []
        for uv in uvals.keys():
            dim = np.where(loc==uvals[uv].size)[0]
            if len(dim)!=1:
                warn.warn("Failed to find dimension for "+uv)
            else:
                dims.append(dim[0])
        dims = np.array(dims)
        dims = dims.max()-dims # last value in dims is 1st dimension
        sdims = dims.argsort()
        dimOrder = np.array(uvals.keys())[sdims]

        return dimOrder

    def calcUniqAttrs(self, augment=True):
        """
        Calculate the unique values for all array attributes
        """

        uvals = {}
        # loop over array attributes
        for attr in self.arrayAttrs:
            # make sure values set
            if getattr(self,attr) is not None:
                # find uniq vals
                temp = self._uniqueAttr(attr, augment=False)
                uvals.update({attr:temp})
        # return or augment
        if augment:
            self.uniqueAttrs = uvals
        else:
            return uvals

    def _uniqueAttr(self, attr, augment=False):
        """
        Calculate the unique values of the given attribute
        """

        # get attr vals
        attrVals = getattr(self, attr)
        # sanity check
        assert attrVals is not None, "Cannot find any values for attr:"+attr
        # make uniq array
        uattr = np.array(list(set(list(attrVals.ravel()))))
        uattr.sort()
        if augment:
            setattr(self,"u"+attr, uattr)
        else:
            return uattr

    def calcResolution(self, wave):
        """
        RH 21/10/16

        Calculating the spectral resolution for different libraries at different wavelengths is complicated.
        This function attempts to use the 'self.resolution' variable with a given wavelength to return a resolution
        (R=lam/dlam) in AA.

        There are strict rules for the format of self.resolution:
        1) It must be a 2-element list
        2) The first element refers to a single R=lam/dlam value for ALL wavelengths.
           In this case, the second element is None!
        3) The second element can be a two element list or tuple.
           In this case, the first element of the list should be None!
           E.g. [None, [3.1, 3.4]] as in BC03
           
           The two numbers refer to the starting and finishing FWHM resolution in AA.
           The value at a given wavelength is then calculated from these values by linear interpolation.
        4) The second element can be a dictionary.
           Then the fist element should be None!
           E.g. [None, {'3200,9500':[3.1,3.5]}]

           The dictionary keyword gives starting and finishing wavelengths in AA
           The dictionary value can be a single value - FWHM in AA
           or it can be a 2-element list/tuple defining the starting and finishing resolution FWHMs in AA.
           As before, linear interpolation is used to estimate the resolution at the given wavelength.

        Note that there is currently no way to define R=lam/dlam for different wavelength ranges - only FWHM resolutions.
           

        Inputs:
        - wave : the wavelength of interest

        Outputs:
        - FWHM : the resolution R=lam/dlam FWHM at the given wavelength 
        
        """
        # make local copy for speed
        R=self.resolution[0]
        FWHM=self.resolution[1]

        # helper functions
        def interpFWHM(FWHM, wave, minlam=None, maxlam=None):
            # using starting and ending resolutions,
            # interpolate between
            if minlam is None: minlam = self.lam.min()
            if maxlam is None: maxlam = self.lam.max()
            lamrange=maxlam-minlam
            FWHMrange = FWHM[1]-FWHM[0]
            resAA = (wave-minlam)/lamrange * FWHMrange + FWHM[0]
            return resAA

        ##########################
        # start if, elif, etc loop
        if R is not None:
            # simple case - single R value
            resAA = wave/R
        elif (isinstance(FWHM,float) or isinstance(FWHM,int)):
            # another simple case - single FWHM given
            resAA = FWHM
        elif ((isinstance(FWHM, list) or isinstance(FWHM,tuple)) and (len(FWHM)==2)):
            # linearly interpolate between two values
            resAA = interpFWHM(FWHM,wave)
        elif (isinstance(FWHM,dict)):
            # dict
            keys = FWHM.keys()
            nkeys=len(keys)
            # cycle through keys / wavelength ranges
            klams=[]
            for k in keys:
                # get wavelength range start/stop
                slams = k.split(",")
                lams=[]
                for l in slams:
                    lams.append(float(l))
                klams.append(np.array(lams))
            # now figure out which key is relevant wavelength
            thisone=np.zeros(nkeys, dtype=np.int)
            for n, kl in enumerate(klams):
                if (wave >= kl.min()) & (wave <= kl.max()):
                    thisone[n]=1
            assert np.sum(thisone)==1, "Failed to find a suitable resolution for this wavelength"
            # get key location
            loc=np.where(thisone==1)[0]
            k = keys[loc]
            # start new if elif else loop
            if isinstance(FWHM[k],float) or isinstance(FWHM[k],int):
                # single value of resAA
                resAA = FWHM
            elif ((isinstance(FWHM[k], list) or (isinstance(FWHM[k],tuple))) and (len(FWHM[k])==2)):
                # interp
                resAA = interpFWHM(FWHM[k],wave, minlam=klams[loc].min(), maxlam=klams[loc].max())
            else:
                raise ValueError("Failed to understand dict type for this key")
        else:
            raise ValueError("Resolution type not understood.")
        
        # return resolution FWHM in AA
        return resAA

    def unraster(self):
        """
        Convert the flam and fmu arrays into simple 2D arrays with the first dim equal to wavelength
        and the second dim indexing the spectrum number.

        """
        # unrasterise spectra - flam and fmu
        self.originalDims = cp.copy(self.dims)
        self.originalFDims = cp.copy(self.flam.shape)
        self.unrasteredAttrs=[] # init a tally of which attrs were unrastered

        # check if spectrum already just 2D
        if len(self.originalFDims)==2:
            warn.warn("Spectrum already in 2D format, nothing to unraster")
        else:
            # make sure flam or fmu exists before reshaping - as we use this function in calclamspec/calcmuspec, where
            # the other type hasn't been created yet
            if hasattr(self,"flam"):
                self.flam = self.flam.reshape((-1,self.nlam))
                if self.eflam is not None:
                    self.eflam = self.eflam.reshape((-1, self.nlam))
            if hasattr(self,"fmu"):
                self.fmu = self.fmu.reshape((-1,self.nmu))
                if self.efmu is not None:
                    self.efmu = self.efmu.reshape((-1, self.nmu))

            self.dims = cp.copy(self.flam.shape[:-1])
            self.fdims = cp.copy(self.flam.shape)

            # can't do this - in case fmu not defined yet
            #assert np.all(np.equal(np.array(self.flam.shape[:-1]), np.array(self.fmu.shape[:-1]))), \
            #       "Unrasterised flam and fmu have different shapes?"

            # unrasterise known (possibly) array attributes
            for attr in cp.copy(self.arrayAttrs):
                self._unrasterAttr(attr)
            # attempt to unrasterise userdict attributes
            if self.__userdict__ is not None:
                for attr in __userdict__.keys():
                    try:
                        self._unrasterAttr(attr)
                    except:
                        warn.warn("Failed to unrasterise user attribute "+str(attr)+".")
            # define the dimensions of the unrastered types
            self.unrasterDims = cp.copy(self.dims)
            self.unrasterFDims = cp.copy(self.fdims)
            
            
            # finally, log the current format as unrasterised
            self.unrastered=True


    def _unrasterAttr(self,attr):
        """
        Helper function to unrasterise attributes such as age, Z, alpha, mass, etc
        """

        # only unraster if not scalar attribute
        if (not np.isscalar(getattr(self,attr))) and (getattr(self,attr) is not None):
            # sanity check
            ndims = len(self.originalDims)
            attrShape = getattr(self, attr).shape
            nAdims = len(attrShape)
            # check for array attribute
            if nAdims == ndims:
                firstCheck = np.all(np.equal(attrShape, self.originalDims)) # attr has shape like age, Z, etc
                assert firstCheck, \
                   attr+" attribute does not have same dimensions as array attribute."
                setattr(self, attr, getattr(self,attr).reshape(self.dims))
                # keep a tally of which attrs were unrastered
                self.unrasteredAttrs.append(attr)
            # check for flux attribute
            elif nAdims==ndims+1:
                secondCheck = np.all(np.equal(attrShape, self.originalFDims))
                assert secondCheck, \
                   attr+" attribute does not have same dimensions as spectra."
                setattr(self, attr, getattr(self,attr).reshape(self.fdims))
                # keep a tally of which attrs were unrastered
                self.unrasteredAttrs.append(attr)
            else:
                raise ValueError("Attr dimensions not understood.")

    def reraster(self, checkSanity=True, debug=False, verbose=False):
        """
        Convert the flam and fmu spectral arrays back into their original, mutli-D formats.
        """
        if debug: pdb.set_trace()

        # check if spectrum is inherently 2D
        if len(self.originalFDims)==2:
            warn.warn("spectum inherently 2D, nothing to reraster")
        else:
            # sanity check
            assert self.unrastered, "Spectra are not unrasterised, so cannot re-rasterise"
            # re-rasterise the flam and fmu - they should ALWAYS exist together at this point.
            newLamShape = list(cp.copy(self.originalDims))
            newLamShape.append(self.nlam)
            tnewLamShape = tuple(newLamShape)
            self.flam = self.flam.reshape(tnewLamShape)
            newMuShape = list(cp.copy(self.originalDims))
            newMuShape.append(self.nmu)
            tnewMuShape = tuple(newMuShape)
            self.fmu = self.fmu.reshape(tnewMuShape)
            # don't forget about the error arrays
            if self.eflam is not None:
                assert self.efmu is not None, "Why does eflam exist, but not efmu?"
                self.eflam = self.eflam.reshape(tnewLamShape)
                self.efmu = self.efmu.reshape(tnewMuShape)

            # sanity check
            if checkSanity:
                assert np.all(np.equal(np.array(self.flam.shape[:-1]), np.array(self.fmu.shape[:-1]))), \
                       "Re-rasterised flam and fmu have different shapes?"
            self.unrastered=False

            # re-rasterise all unrastered attrs
            for attr in cp.copy(self.unrasteredAttrs):
                self._rerasterAttr(attr, debug=debug)
                if verbose: print "Rerastered "+attr

            if debug: pdb.set_trace()
            # sanity check - no more unrastered attrs
            assert len(self.unrasteredAttrs)==0, "Still have unrastered attributes?"
            # reset dims - don't use originalDims here, as flux dims might have changed between rastering (interpOnNewLam)
            self.fdims = cp.copy(self.flam.shape)
            self.dims = cp.copy(self.fdims[:-1])
            # delete the redundant var
            delattr(self, "unrasterDims")
            delattr(self, "unrasterFDims")
            
            # set as re-rastered
            self.unrastered=False
        # delete more redundant vars, which were needed for 2D spec case
        delattr(self, "originalDims")
        delattr(self, "originalFDims")
        delattr(self, "unrasteredAttrs")

    def _rerasterAttr(self, attr, debug=False):
        """
        Helper function to re-rasterise the array attributes
        """
        # only unrasterise if not scalar attribute.
        # Although this should never be a problem, as _rerasterise is called with args from self.unrasteredAttrs
        
        if (not np.isscalar(getattr(self,attr))) and (getattr(self,attr) is not None):
            if debug: pdb.set_trace()
            ndims = len(self.unrasterDims)
            attrShape = getattr(self, attr).shape
            nAdims = len(attrShape)
            # check for array attrs
            if nAdims==ndims:
                # attr is like an age, Z or alpha: no wavelength info
                # sanity check
                firstCheck = np.all(np.equal(attrShape,self.unrasterDims)) # attr has shape like age, Z, etc
                assert firstCheck, \
                   attr+" attribute does not have same dimensions as array attribute."
                setattr(self, attr, getattr(self,attr).reshape(self.originalDims))
                # remove this attr from the unraster list
                self.unrasteredAttrs.remove(attr)
            elif nAdims == ndims+1:
                # attr has wavelength info, like flam, eflam etc
                secondCheck = np.equal(attrShape[0], self.unrasterFDims[0]) # check only product of dims - not nlam
                assert secondCheck, \
                   attr+" attribute does not have same dimensions as spectra."
                newFDims = list(cp.copy(self.originalDims))
                newFDims.append(attrShape[1])
                setattr(self, attr, getattr(self,attr).reshape(newFDims)) # Don't use self.fdims or self.originalFDims here (e.g. interpOnNewLam)
                # remove this attr from the unraster list
                self.unrasteredAttrs.remove(attr)
            else:
                raise ValueError("Attr dimensions not understood.")
            
    def getShape(self):
        """
        Return the shape of the flam/fmu arrays, ignoring the spectral dim
        """
        return self.dims

    def getSize(self):
        """
        Return the total number of spectra
        """
        return np.prod(self.dims)

    def calcmuspec(self, filter=False, debug=False):
        """
        Calculate a muspec (GJy @ GHz) from a lamspec (erg/s/cm**2/AA @ AA)
        """
        # (convert AA to GHz with c in m/s): (c[m/s] / (lam[AA] * 1e-10) / 1e9 = c/AA * 1e1)
        self.mu = c / (self.lam) * 1e1
        # erg/s/cm**2/AA => *LAM[AA]*LAM[M]  / C[M/S] * 1e23 = LAM[AA]*LAM[AA]/ c[m/s] * 1e13 => GJy=1e-14 erg/s/cm**2/Hz (Jy=1e-23erg/s/cm2/Hz=1e-26W/m2/Hz)

        if debug: pdb.set_trace()
        # unraster - convert flam into 2D
        self.unraster()
        
        self.fmu = [] # init
        if filter:
            for flam in self.flam:
                self.fmu.append(flam)
            self.fmu=np.atleast_2d(singleOrList2Array(self.fmu))
        else:
            for flam in np.atleast_2d(self.flam):
                self.fmu.append( flam * self.lam * (self.lam) / c * 1e4 )
            self.fmu = np.atleast_2d(singleOrList2Array(self.fmu)) # at least 2D in case single spec
            if self.eflam is not None: 
                self.efmu = []
                for eflam in self.eflam:
                    if eflam is not None:
                        self.efmu.append( eflam * self.lam * (self.lam) / c * 1e4 )
                self.efmu = np.atleast_2d(singleOrList2Array(self.efmu)) # at least 2D in case single errspec
            else:
                self.efmu = None

        # assign dims of spectral component
        loc = self.fmu.shape
        self.nmu = loc[-1]

        # put back into multi-D format
        self.reraster(debug=debug)

    def calclamspec(self, filter=False):
        """
        Calculate a lamspec (erg/s/cm**2/AA @ AA) from a muspec (GJy @ GHz)
        """
        # (convert GHz to AA with c in m/s): c[m/s] / (mu[GHz]*1e9) * 1e10 = c[m/s] / mu[GHz] * 1e1
        self.lam = c / (self.mu) * 1e1
        # GJy => / LAM(AA), /LAM(M)  * C(M/S) => erg/s/cm**2/Hz

        # unraster - convert fmu into 2D
        self.unraster()

        self.flam = [] # init
        if filter:
            for fmu in self.fmu:
                self.flam.append(fmu)
            self.flam = np.atleast_2d(singleOrList2Array(self.flam))
        else:
            for fmu in self.fmu:
                self.flam.append( fmu / self.lam / (self.lam) * c / 1e4 )
            self.flam = np.atleast_2d(singleOrList2Array(self.flam)) # at least 2D in case single spec
            if self.efmu is not None: 
                self.eflam=[]
                for efmu in self.efmu:
                    if efmu is not None:
                        self.eflam.append( efmu / self.lam / (self.lam) * c / 1e4 )
                self.eflam = np.atleast_2d(singleOrList2Array(self.eflam)) # at least 2D in case single errspec
            else:
                self.eflam = None

        # assign dims of spectral component
        loc = self.flam.shape
        self.nlam = loc[-1]

        # put back into multi-D format
        self.reraster()

    def rebinLam(self, dLam=1e-4, flux=False):
        """
        Rebin the spectrum to a fixed/minimum dLam

        Flux - do rebinning by conserving flux
        """

        newflams=[]
        newlam = np.arange(nn(self.lam.min(),dLam,ceil=True), nn(self.lam.max(),dLam,floor=True)+dLam, dLam)
        self.unraster()
        for flam in self.flam:
            newflam = rebin(self.lam,flam,newlam,flux=False)
            newflams.append(newflam)

        self.lam=singleOrList2Array(newlam)
        self.flam=np.atleast2d(singleOrList2Array(newflams))
        self.reraster()
        self.calcmuspec()

    def interpOnNewLam(self, newlam, imethod='linear'):
        """
        Interpolate the spectra onto a new lambda grid.

        """
        self.newlam = newlam
        # put into 2D mode
        self.unraster()
        newflam = [] # init
        # loop over spectra, interpolating onto new grid
        for fl in self.flam:
            newflam.append(interpolate(self.lam, fl, newlam, method=2, kind=imethod))
        # put as new attrib
        self.newflam = np.atleast_2d(singleOrList2Array(newflam))
        # add to unrastered list
        self.unrasteredAttrs.append("newflam")
        
        # process error array in same way if exists
        if self.eflam is not None:
            neweflam = []
            for ef in self.eflam:
                neweflam.append(interpolate(self.lam, ef, newlam, method=2, kind=imethod))
            self.neweflam = np.atleast_2d(singleOrList2Array(neweflam))
            self.unrasteredAttrs.append("neweflam")
            assert np.all(self.newflam.shape==self.neweflam.shape), 'new flux and error arrays different sizes?'
        else:
            self.neweflam=None
            neweflam=None

        # update unrastered Flux dims
        #newFDims = self.newflam.shape
        #assert self.unrasterFDims[0]==newFDims[0], 'Changed dimensions, somehow.'
        #self.unrasterFDims = newFDims
        #self.nlam = newFDims[-1] # be sure to update nlam here as this is used to construct dims in reraster
        # reraster, incl newlam (and neweflam), but with no checks against sizes of flam and fmu
        self.reraster(checkSanity=False)
        # copy back & process
        self._processLamspec(self.newlam, self.newflam, errlamspec=self.neweflam) # need self.[vars] here as been rerastered (except for newlam)
        # construct new muspec to be consistent
        self.calcmuspec(filter=filter)

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
        self.unraster()
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

        self.reraster()
        
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
        self.unraster()
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

        # take care of rastering
        self.mags = mags
        self.unrasteredAttrs.append("mags")
        self.reraster()
            
        return self.mags

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
        mags=[]
        self.unraster()
        for flam in self.flam:
            iflam = interpolate(self.lam, flam, lam)

            # integrate using simple trapezium method
            flux = np.sum(iflam*ithro*lam) / np.sum(ithro*lam)

            # calc mag: remember that flux in erg/s/cm**2/AA which is official unit of ST mags (no correction)
            mags.append(-2.5*np.log10(flux) - 21.10)

            if plot:
                pl.plot(lam,iflam/np.median(iflam))
                pl.plot(lam,ithro/np.median(ithro))

        mags = np.array(mags)
        self.unrasteredAttrs.append("mags")
        self.reraster()

        return self.mags
    
  
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

        mags = selfmag-vegamag
        
        return mags

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
        setattr(self,"Lsuns", 10.0**(-0.4*l) / 10.0**(-0.4*lsun))
        
        m2l = (m/10.0**(-0.4*l)) / (msun/10.0**(-0.4*lsun))

        return m2l

    def gaussLamConvolve(self, sigma_lam, nsig=5.0, overwrite=True, verbose=True):
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
        self.unraster()
        
        for flam in self.flam:
            count+=1
            reg_spec.append(interpolate(self.lam,flam,reg_lam, fill_value=np.nan, \
                                        bounds_error=False, method=1, kind='linear'))
            if verbose: print "Interpolated spec "+str(count)+" of "+str(len(self.flam))+" onto regular wavelength grid"
            
        # make kernel
        krange = np.ceil((sigma_lam*nsig)/min_dlam) # in units of min_dlam
        kx = np.arange(-krange,krange+1) # need +1 here otherwise kernel is not even and it adds shift
        w = kx/(sigma_lam/min_dlam)
        ky = np.exp(-0.5*w**2.0)
        # normalise
        ky/=np.sum(ky)

        # do convolution
        count=0
        creg_spec=[]
        for rs in reg_spec:
            count+=1
            creg_spec.append(np.convolve(rs,ky,'same'))
            if verbose: print "Convolved spec "+str(count)+" of "+str(len(reg_spec))

        # interpolate back to original sampling
        cspec=[]
        count=0
        for crs in creg_spec:
            count+=1
            cspec.append(interpolate(reg_lam,crs,self.lam, fill_value=np.nan, bounds_error=False, method=1, kind='linear'))
            if verbose: print "Interpolated spec "+str(count)+" of "+str(len(creg_spec))+" onto regular wavelength grid"

        cspec = np.atleast_2d(singleOrList2Array(cspec))
        self.cflam = cspec

        if overwrite:
            self.flam = np.copy(self.cflam)
        
        self.unrasteredAttrs.append("cflam")
        self.reraster()

        return cspec

    def calcLogLamGrid(self, velscale=None, reraster=True, verbose=False):
        """
        Put the spectra onto a loglam grid
        """
        
        # get minimum vel resolution in spec and use this for dloglam
        #if velscale is None:
        #    velscale = np.mean((self.lam[1:]-self.lam[:-1])/self.lam[:-1]) * c / 1e3 #
        
        
        # calc regular loglam grid
        if velscale is None:
            self.nloglam = self.nlam 
            self.loglam = np.exp(np.linspace(np.log(self.lam[0]), np.log(self.lam[-1]), nloglam ))
            dloglam = (np.log(self.lam[-1])-np.log(self.lam[0]))/(nloglam-1.)
            self.velscale = (np.exp(dloglam)-1.) * c
        else:
            dloglam = np.log(1.0 + velscale/c*1e3)
            self.loglam = np.exp(np.arange(np.log(self.lam[0]), np.log(self.lam[-1]), dloglam))
            self.nloglam = len(self.loglam)

        count=0
        self.floglam=[]
        self.unraster()
        
        for flam in self.flam:
            count=count+1
            # interpolate onto loglam grid
            self.floglam.append(interpolate(self.lam,flam,self.loglam, fill_value=np.nan, \
                                            bounds_error=False, method=1, kind='linear'))
            if verbose: print "Interpolated spec "+str(count)+" of "+str(len(self.flam))
        self.floglam = np.atleast_2d(singleOrList2Array(self.floglam))

        if reraster:
            # manually re-raster floglam because it might have different nlam to original spec
            self.floglam = np.array(self.floglam)
            loglamDims = list(cp.deepcopy(self.originalDims))
            loglamDims.extend([self.nloglam])
            self.floglam = self.floglam.reshape(loglamDims)
            self.reraster()
        
    def gaussVelConvolve(self, vel, sigma, h3h4=None, correct4InstRes=True, nsig=5.0, losvd=None, overwrite=True, verbose=True):
        """
        Purpose: to convolve the spectrum with a Gaussian of known velocity (V)
                 and width (SIGMA)

        Input:
           - vel  : velocity (use 0.0 for now)
           - sigma: dispersion in km/s      

        """
        global c

        # get minimum vel resolution in spec and use this for dloglam
        velscale = np.mean((self.lam[1:]-self.lam[:-1])/self.lam[:-1]) * c / 1e3 #np.min((self.lam[1:]-self.lam[:-1])/self.lam[:-1]) * c / 1e3 # km/s
        dloglam = np.log10(1.0 + velscale/c*1e3)
        nloglam = self.nlam #np.round((np.log10(self.lam.max())-np.log10(self.lam.min())) / dloglam)
        self.velscale=velscale
        
        # calc regular loglam grid
        self.loglam = 10.0**np.linspace(np.log10(self.lam[0]), np.log10(self.lam[-1]), nloglam )

        count=0
        self.floglam=[]
        self.unraster()
        
        for flam in self.flam:
            count=count+1
            # interpolate onto loglam grid
            self.floglam.append(interpolate(self.lam,flam,self.loglam, fill_value=np.nan, \
                                            bounds_error=False, method=1, kind='linear'))
            if verbose: print "Interpolated spec "+str(count)+" of "+str(len(self.flam))
        self.floglam = np.atleast_2d(singleOrList2Array(self.floglam))

        # if kernel not passed, create it.
        if losvd == None: # speed up if losvd passed
            
            # calc required kernel dispersion
            if correct4InstRes:
                # correct for instrumental dispersion
                meanLam = np.mean(np.array([self.lam[0],self.lam[-1]])) # mean wavelength
                resAA = self.calcResolution(meanLam) # FWHM res at this wavelength
                sigmaSpec = resAA/meanLam * c / 1e3 / np.sqrt(8.*np.log(2.)) # equivalent vel disp
                # sanity check
                assert sigma > sigmaSpec, "Cannot convolve to this dispersion, spectral resolution too coarse."
                # do calc
                sigmaKernel = np.sqrt(sigma**2.0 - sigmaSpec**2.0)
            else:
                sigmaKernel = sigma
            
            dv = np.ceil(nsig*sigmaKernel/velscale) 
            nv = 2*dv + 1
            v = np.linspace(dv,-dv,nv) 
            w = (v - vel/velscale) / (sigmaKernel/velscale)
            w2= w*w
            if h3h4 != None:
                h3=h3h4[0]
                h4=h3h4[1]
                poly = 1.0 + h3/np.sqrt(3.0)*(w*(2.0*w2-3.0)) + \
                       h4/np.sqrt(24.0)*(w2*(4.0*w2-12.0)+3.0)
            else:
                poly = np.ones(nv)

            losvd = np.exp(-0.5*w**2.0)/(np.sqrt(2.0*np.pi)*sigmaKernel/velscale) * poly 
            losvd=losvd/np.sum(losvd)

        count=0
        self.confloglam=[]
        self.conflam = []
        for floglam in self.floglam:
            count=count+1
            confloglam = np.convolve(floglam,losvd,'same')
            self.confloglam.append(confloglam)
            # interpolate back onto origial grid
            self.conflam.append(interpolate(self.loglam, confloglam, self.lam, fill_value=np.nan, bounds_error=False, method=1, kind='linear'))
            if verbose: print "Convolved spec "+str(count)+" of "+str(len(self.floglam))

        self.confloglam = np.atleast_2d(singleOrList2Array(self.confloglam))
        self.conflam = np.atleast_2d(singleOrList2Array(self.conflam))
        
        
        if overwrite:
            self.flam = self.conflam

        self.unrasteredAttrs.extend(["confloglam", "conflam"])
        self.reraster()

        return self.conflam

    def clipSpectralRange(self, minlam, maxlam):
        """
        Clip the spectral range to be between minlam and maxlam
        """

        newflams = []
        # set range
        loc = np.where((self.lam>minlam) & (self.lam<maxlam))[0]
        # clip specs
        self.unraster()
        
        for flam in self.flam:
            newflams.append(flam[loc])
        if self.eflam is not None: 
            neweflams=[]
            for eflam in self.eflam:
                neweflams.append(eflam[loc])
        
        # overwite old flams
        self.lam  = self.lam[loc]
        self.flam = np.atleast_2d(singleOrList2Array(newflams))
        if self.eflam is not None: 
            self.eflam = np.atleast_2d(singleOrList2Array(neweflams))

        self.nlam = loc.shape[0]

        self.reraster()

        # update muspec
        self.calcmuspec()
        
    def calcIndex(self, index, disp=None, round_prec=10, method=0, verbose=False):
        """
        Function to calculate absorption line index strength.

        Inputs:
            index: an ind class from indexTools.py
            disp: float - spectral dispersion [A/pixel] (i.e. pixel sampling),
                          NOT anything to do with spectral (FWHM) resolution!!!
                          NOT anything to do with velocity dispersion !!!
            round_prec: precision for rounding wavelengths (Method=4)
            Method: Specify the index calculation method
                    0 = Traditional/Simple method using means (no variance weighting)
                    1 = Simple method using median (no variance weighting)
                    3 = Linear fit, no variance weighting
                    4 = Full inverse varaince weighting and propagation of errors (Cenarro 2001)

        Outputs:
            ind: value of chosen index
            ind_var (optional): variance on index measurement if spectrum class has variance defined
            
        """
        # check that both the spectrum and the index are defined in the same wavelength system
        if (self.wavesyst is not None) & (index['wavesyst'] is not None):
            if self.wavesyst!=index['wavesyst']:
                raise ValueError("Spectrum and index defined on different wavelength systems")
        else:
            warn.warn("Unsure of wavelength systems for index and/or spectrum (air or vac?).")

        # check the resolution
        if index['resol'] is not None:
            outputFWHM = index['resol']
            newSpec = cutAndGaussLamConvolve(self, index, outputFWHM, verbose=False)
        else:
            newSpec = self
        
        # loop through various methods
        if method==0:
            vals = calcSimpleIndex(newSpec, index, contMethod='mean', verbose=verbose)
        elif method==1:
            vals = calcSimpleIndex(newSpec, index, contMethod='median', verbose=verbose)
        elif method==2:
            raise ValueError("Code not written")
        elif method==3:
            vals = calcCenarroIndex(newSpec, index, disp=disp, round_prec=round_prec, verbose=verbose)
        else:
            raise ValueError("Method not understood")

        return np.array(vals)

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


    def normaliseSpec(self, polyOrder=3, indLib=None, maxErrVal=999.0, overwrite=True, \
                      keepPoly=False):
        """
        Normalise spectrum continuum to one, avoiding abs lines

        """

        nfs=[]
        nefs=[]
        pfs=[]

        self.unraster()
        nspec = self.nspec

        # loop over specs, normalising
        for ni in xrange(nspec):
            # take local copies of flam and eflam (if defined)
            f=self.flam[ni]
            if self.eflam is not None: #hasattr(self,"eflam"):
                ef = self.flam[ni]
            else:
                ef = None
            # apply normalisation
            rl = normaliseSpec(self.lam, f, eflam=ef, polyOrder=polyOrder, indLib=indLib, \
                              maxErrVal=maxErrVal, returnPoly=keepPoly)
            nfs.append(rl[0])
            # sort additional outputs
            listindex=0
            if ef is not None:
                listindex+=1
                nefs.append(rl[listindex]) # error
            if keepPoly:
                listindex+=1
                pfs.append(rl[listindex])  # polynomial fit

        # save copied to object
        self.nflam = np.atleast_2d(singleOrList2Array(nfs))
        
        if len(nefs)!=0:
            self.neflam = np.atleast_2d(singleOrList2Array(nefs))
            self.unrasteredAttrs.append('neflam')
        if keepPoly:
            self.pfits = np.atleast_2d(singleOrList2Array(pfs))
            self.unrasteredAttrs.append('pfits')
        
        if overwrite:
            self.flam = self.nflam
            if len(nefs)!=0: self.eflam = self.neflam

        self.reraster() # put into multiD format


###################################### END OF SPECTRUM CLASS #####################################

def checkDims(var, varname, parentShape):
    """
    Check if the dimensions of the attribute match the dimension of the parent
    """
    assert (np.isscalar(var)) or (np.all(np.equal(np.array(var).shape,parentShape))), varname+" dimensions not understood."

def singleOrList2Array(invar):
    """
    If a single value, leave. If a list, convert to array. If array, leave as array.
    But for size-1 arrays/lits, convert back to scalar.
    """

    if isinstance(invar, list):
        # convert to array unless size-1
        rval = np.squeeze(np.array(invar))
        if rval.size==1:
            rval = np.asscalar(rval)
    elif isinstance(invar, np.ndarray):
        # leave except if size-1
        rval = np.squeeze(invar)
        if rval.size==1:
            rval = np.asscalar(rval)
    else:
        # leave
        rval = invar
    # return
    return rval
        

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

def cutAndGaussLamConvolve(longSpec, index, outputFWHM, currentFWHM=None, nSig=5.0, verbose=True):
    """
    RH 18/10/2016

    In order to measure Lick/IDS indices, we need to convolve the spectra up to a particular resolution
    and measure the index strength.

    This code cuts out a part of the spectrum (red and blue continuum limits +- 5 sigma) and convolves
    the spectrum with the given dispersion. It then returns this new, shorter, convolved spectrum.

    Puzia et al 2004 Eq in Sec4.1 is useful

    All convolutions done in lambda space, not velocity space.

    Inputs:
     - longSpec   : the spectrum which is to be cut and convolved
     - index      : the index definition, used to cut spec to right size
     - currentFWHM: the current spectral resolution (FWHM) in AA
     - outputFWHM : the desired output resolution (FWHM) in AA
     - verbose    : print stuff, can get anoying

    """
    global sigmaInFWHM
    
    # sort inputs
    if currentFWHM is None:
        meanWave = np.mean(np.array([index['ind_start'],index['ind_stop']]))
        currentFWHM = longSpec.calcResolution(meanWave)

    # calc convolution kernel width, using Puzia+04 Eq
    kernelFWHM = np.sqrt(outputFWHM**2.0 - currentFWHM**2.0)
    kernelSigma=kernelFWHM/sigmaInFWHM
    assert outputFWHM > currentFWHM, "Cannot convole to a lower resolution than already have."
    if verbose: print "In FWHM="+str(currentFWHM)+", out FWHM="+str(outputFWHM)+", kernel FWHM="+str(kernelFWHM)

    # define cut spectrum edges
    cutlow = (index['blue_start'] - nSig*kernelSigma)
    cuthigh= (index['red_stop']   + nSig*kernelSigma)

    # cut spec
    cutSpec = longSpec.copy()
    cutSpec.clipSpectralRange(cutlow, cuthigh)

    # convolve
    cutSpec.gaussLamConvolve(kernelSigma, overwrite=True, verbose=verbose)
    # return
    return cutSpec


def cutAndGaussVelConvolve(longSpec, index, outputSigma, currentFWHM=None, doPlot=False, verbose=True):
    """
    RH 18/10/2016

    In order to measure indices at different velocity dispersions, we need to cut and convolve a
    spectrum with a Gaussian kernel in velocity space.

    This code cuts out a part of the spectrum (red and blue continuum limits +- 5 sigma) and convolves
    the spectrum with the given dispersion, in velocity space. It then returns this new, shorter, convolved spectrum.

    Puzia et al 2004 Eq in Sec4.1 is useful

    All convolutions done in velocity space, not lambda space

    Inputs:
     - longSpec   : the spectrum which is to be cut and convolved
     - index      : the index definition, used to cut spec to right size
     - outputSigma: the desired output velocity dispersion in km/s
     - currentFWHM: the current spectral resolution (FWHM) in AA
     - verbose    : print stuff, can get anoying

    """
    global sigmaInFWHM, c
    meanWave = np.mean(np.array([index['ind_start'],index['ind_stop']]))
    # sort inputs
    if currentFWHM is None:
        currentFWHM = longSpec.calcResolution(meanWave)
        
    # calc convolution kernel width, using Puzia+04 Eq
    outputFWHM = outputSigma*1e3/c * meanWave * sigmaInFWHM
    kernelFWHM = np.sqrt(outputFWHM**2.0 - currentFWHM**2.0)
    kernelSigma= kernelFWHM/sigmaInFWHM
    assert outputFWHM > currentFWHM, "Cannot convole to a lower resolution than already have."
    if verbose: print "In FWHM="+str(currentFWHM)+", out FWHM="+str(outputFWHM)+", kernel FWHM="+str(kernelFWHM)

    # define cut spectrum edges
    cutlow = (index['blue_start'] - 5.0*kernelSigma)
    cuthigh= (index['red_stop']   + 5.0*kernelSigma)

    ## cut spec
    cutSpec = longSpec.copy()
    cutSpec.clipSpectralRange(cutlow, cuthigh)

    # convolve
    velSigma = kernelSigma/meanWave*c/1e3
    cutSpec.gaussVelConvolve(0.0, velSigma, correct4InstRes=False, overwrite=True, verbose=verbose)

    # do plot if asked
    if doPlot:
        pl.figure(figsize=(12,6))
        pl.subplot(121)
        lnorm = np.median(longSpec.flam[0][cloc])
        pl.plot(longSpec.lam, longSpec.flam[0]/lnorm, "k-")
        pl.plot(cutSpec.lam, cutSpec.flam[0]/lnorm, "r-")
        pl.fill_between([index['ind_start'][0], index['ind_stop'][0]], 0., [2.0, 2.0], color="blue", alpha=0.2)
        pl.fill_between([index['blue_start'], index['blue_stop']], 0., [2.0, 2.0], color="red", alpha=0.2)
        pl.fill_between([index['red_start'], index['red_stop']], 0., [2.0, 2.0], color="red", alpha=0.2)
        pl.axis([index['blue_start']*0.995, index['red_stop']*1.005, 0.5, 1.3])
        pl.subplot(122)
        lnorm = np.median(longSpec.flam[-1][cloc])
        pl.plot(longSpec.lam, longSpec.flam[-1]/lnorm, "k-")
        pl.plot(cutSpec.lam, cutSpec.flam[-1]/lnorm, "r-")
        pl.fill_between([index['ind_start'][0], index['ind_stop'][0]], 0., [2.0, 2.0], color="blue", alpha=0.2)
        pl.fill_between([index['blue_start'], index['blue_stop']], 0., [2.0, 2.0], color="red", alpha=0.2)
        pl.fill_between([index['red_start'], index['red_stop']], 0., [2.0, 2.0], color="red", alpha=0.2)
        pl.axis([index['blue_start']*0.995, index['red_stop']*1.005, 0.5, 1.3])
        
    # return new shorter and convolved spectrum
    return cutSpec
    

def calcSimpleIndex(spectrum, index, contMethod='mean', disp=None, round_prec=10, \
                    continIncludesPartialPix=True, doPlot=False, verbose=False):
    """
    RH 18/10/2016

    Calculate a simple equivalent width index for absorption lines.

    Inputs:
    =======
    - spectrum: a spectrum class object
    - index: an index class object
    - contMethod: specify whether to use MEAN or MEDIAN averaging along continuum bandpass
    - disp: spectrum dispersion in AA/pix. If None, automatically find it.
    - round_prec: rounding precision for calculating disp
    - continIncludesPartialPix: if bandpass falls within a pixel (which it will), do we include this pixel in
                                the continuum calculation?
    - verbose: print info to screen as running

    Outputs:
    =======
    - index: the measured index
    - indexVAR (optional)

    """

    # sanity check that this is a simpleIndex
    assert index['simpleIndex'], "Sorry, index is not a 'Simple Index', so cannot perform simple calculation of the index."

    # check for uniform spacing
    delta = checkSpacing(spectrum, index, round_prec=round_prec)
    if type(disp)==type(None): disp=delta

    if spectrum.eflam is not None:
        calcVar=True
    else:
        calcVar=False
    

    # if no variance array
    # if not hasattr(spectrum,'eflam'):
    # loop over spectra of different ages
    indVals=[] # init
    if calcVar:
        EindVals = []

    spectrum.unraster() # put into 2D format
    
    for ns in xrange(spectrum.nspec):
        # loop over continuum definitions
        yAv=[] # init
        xAv=[]
        if calcVar:
            VyAv = [] # init variance array
        for j in xrange(index['ncont']):
            # Find first and last data indices within bandpass
            a = np.where(spectrum.lam > index['cont_start'][j])[0][0]
            b = np.where(spectrum.lam < index['cont_stop'][j])[0][-1]
            # Make sure the pixels containing the bandpass edges are included
            if (index['cont_start'][j] - spectrum.lam[a-1]) < (spectrum.lam[a] - index['cont_start'][j]):
                a -= 1
            if (index['cont_stop'][j] - spectrum.lam[b]) > (spectrum.lam[b+1] - index['cont_stop'][j]):
                b += 1

            # Multiplicative factors for start and end pixels
            Cstart_c = (spectrum.lam[a] - index['cont_start'][j] + 0.5*disp)/disp  
            Cend_c = (index['cont_stop'][j] - spectrum.lam[b] + 0.5*disp)/disp 

            # make weighted continuum 
            wCont = np.ones_like(spectrum.flam[ns][a:b+1])
            wCont[0]=Cstart_c
            wCont[-1]=Cend_c
            wCsum = np.sum(wCont)
                        
            # take average
            if contMethod=='mean':
                # use (weighted) mean
                yAv.append(np.sum(spectrum.flam[ns][a:b+1]*wCont)/wCsum)
                xAv.append(np.sum(spectrum.lam[a:b+1]*wCont)/wCsum)

                if calcVar:
                    VyAv.append(np.sum( (spectrum.eflam[ns][a:b+1]*wCont/wCsum)**2.0 ))

            elif contMethod=='median':
                yAv.append(np.median(spectrum.flam[ns][a:b+1]))
                xAv.append(np.median(spectrum.lam[a:b+1]))

                if calcVar:
                    VyAv.append( 1.253**2.0 * np.sum( (spectrum.eflam[ns][a:b+1]*wCont/wCsum)**2.0 ) )
            else:
                raise ValueError("contMethod not understood")

        assert len(yAv)==2, "This should never happen: Can't calculate simple index when More than 2 continuum regions defined."

        # convert to arrays
        yAv = np.array(yAv); xAv = np.array(xAv);
        if calcVar: VyAv=np.array(VyAv)
        
        # calc linear continuum: y=mx+c using simple simultaneous equn solution
        deltaY = yAv[1]-yAv[0]
        deltaX = xAv[1]-xAv[0]
        gradient = (deltaY) / (deltaX)
        intercept    = yAv[0]-gradient*xAv[0]

        if calcVar:
            # calc errors on grad and intercept - see notes 2/11/16
            Vgradient = (np.sum(VyAv)/deltaX**2.)  # no x error
            VCs = VyAv # no x error
            VCinvsum = np.sum(1.0 / VCs)
            alphas = VCs/VCinvsum
            Vintercept = np.sum(alphas**2. * VCs)

        #######################################
        # Calculate the index
        assert index['nfeat']==1, "Cannot calculate simple index for index with more than one feature definition."
        # init
        ind = 0.; ind_var_tot = 0.
        ind_var_1 = 0.; ind_var_2 = 0.
        
        # Find first and last data indices within bandpass
        a = np.where(spectrum.lam > index['ind_start'])[0][0]
        b = np.where(spectrum.lam < index['ind_stop'])[0][-1]
        # Determine which pixel is closest to start and end of bandpass
        # RH doc: determine if a,b pixels currently overlap with bandpass edges or not
        # if not, extend a,b so that they do. 
        if (index['ind_start'] - spectrum.lam[a-1]) < (spectrum.lam[a] - index['ind_start']):
            a -= 1
        if (index['ind_stop'] - spectrum.lam[b]) > (spectrum.lam[b+1] - index['ind_stop']):
            b += 1
        # Multiplicative factors for start and end pixels
        Cstart_c = (spectrum.lam[a] - index['ind_start'] + 0.5*disp)/disp  
        Cend_c = (index['ind_stop'] - spectrum.lam[b] + 0.5*disp)/disp 

        Fi = gradient*spectrum.lam[a:b+1] + intercept
        contNormSubSpec = 1.0 - spectrum.flam[ns][a:b+1] / Fi
        contNormSubSpec[0] *= Cstart_c
        contNormSubSpec[-1] *= Cend_c
        indVal = disp*np.sum(contNormSubSpec)
        indVals.append(indVal)

        if calcVar:
            # calc error on index, folding in uncertainty on continuum (ignoring covariances)
            VFi = spectrum.lam[a:b+1] * Vgradient + Vintercept
            VIi = indVal**2. * ( (spectrum.eflam[ns][a:b+1]/spectrum.flam[ns][a:b+1])**2.0 + (VFi/Fi**2.) )
            EindVal = np.sqrt(np.sum(VIi))
            EindVals.append(EindVal)

    spectrum.ind = np.array(indVals)
    spectrum.unrasteredAttrs.append("ind")
    if calcVar:
        spectrum.eind = np.array(EindVals)
        spectrum.unrasteredAttrs.append("eind")
        
    spectrum.reraster() # put into multi-D format
    

    # determine what to return
    rlist = np.copy(spectrum.ind)
    if calcVar:
        rlist=[rlist]
        rlist.append(spectrum.eind)

    return rlist


def checkSpacing(spectrum, index, round_prec=10):
    """
    Check for uniform wavelength spacing. Code taken from SZ's index calculator
    """
    # check uniform spacing of data in wavelegth
    data_loc_start = np.where(spectrum.lam > np.min(index['cont_start']))[0]
    data_start = data_loc_start[data_loc_start.argmin()]#spectrum.lam[data_loc_start].argmin()
    data_loc_stop = np.where(spectrum.lam < np.max(index['cont_stop']))[0]
    data_stop = data_loc_stop[data_loc_stop.argmax()]
    deltas = np.round(spectrum.lam[data_start+1:data_stop]-spectrum.lam[data_start:data_stop-1],round_prec)
    delta = np.median(deltas)
    
    assert np.all(deltas==delta), "Wavelength spacing not uniform"
    return delta


def calcCenarroIndex(spectrum, index, disp=None, round_prec=10, verbose=False):
    """
    Function to calculate absorption line index strength using inverse variance weighting

    -Code originally written by S. Zieleniewski, adopting
    the Cenarro 2001 method of variance weighted fitting.
    -Later adopted by R. Houghton. 

    Inputs:
        spectrum: a spectrum class object 
        index: an indlib class from indexTools.py
        disp: float - spectral dispersion [A/pixel], NOT anythign to do with spectral resolution!
        round_prec: precision for rounding wavelengths.


    Outputs:
        ind: value of chosen index
        ind_var (optional): variance on index measurement if spectrum class has variance defined

    """
    # check for unifrm spacing
    delta = checkSpacing(spectrum, index, round_prec=round_prec)
    if type(disp)==type(None): disp=delta

    spectrum.unraster() # put into 2D format
    
    # init
    indices = np.zeros(spectrum.nspec,dtype=float)
    index_vars = np.zeros_like(indices)
    # loop over spectra and calc indices
    for spec in xrange(spectrum.nspec):

        SED = np.column_stack((spectrum.lam, spectrum.flam[spec]))
        if spectrum.eflam is not None: #hasattr(spectrum, 'eflam'):
            var_SED = np.column_stack((spectrum.lam, spectrum.eflam[spec]**2.0)) # make variance array
        else:
            var_SED = None

        #Method for calculating pseudo-continuum using equations from
        #Cenarro et al (2001) Appendix A.

        # RH doc: fit the continuum according to Eqs A11-A19
        # If no variance array, no varaince weighting (i.e. use VAR=1.0 for everything)
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

                # RH doc: make a copy of the whole pixels in bandpass
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
        # RH doc: doesn't look like fractional pixels included in continuum fit

        # If variance array exists, use it:
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

        # Calculate the index
        # init
        ind = 0.; ind_var_tot = 0.
        ind_var_1 = 0.; ind_var_2 = 0.
        for j in xrange(index['nfeat']):
            #Find first and last data indices within bandpass
            a = np.where(SED[:,0] > index['ind_start'][j])[0][0]
            b = np.where(SED[:,0] < index['ind_stop'][j])[0][-1]
            #Determine which pixel is closest to start and end of bandpass
            # RH doc: determine if a,b pixels currently overlap with bandpass edges or not
            # if not, extend a,b so that they do. 
            if (index['ind_start'][j] - SED[a-1,0]) < (SED[a,0] - index['ind_start'][j]):
                a -= 1
            if (index['ind_stop'][j] - SED[b,0]) > (SED[b+1,0] - index['ind_stop'][j]):
                b += 1
            # Multiplicative factors for start and end pixels
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


    spectrum.ind = np.array(indices)
    spectrum.unrasteredAttrs.append("ind")
    if var_SED is not None:
        spectrum.eind = np.sqrt(np.array(index_vars))
        spectrum.unrasteredAttrs.append("eind")
        
    spectrum.reraster()

    rlist = spectrum.ind
    if var_SED is not None:
        rlist = [spectrum.ind,spectrum.eind]

    return rlist 

    
def normaliseSpec(lam, flam, eflam=None, polyOrder=3, indLib=None, maxErrVal=999.0, returnPoly=False):
    """
    Normalise spectrum continuum to one, avoiding abs lines (defined by indLib)

    Inputs:
      lam     - wavelength array (1D)
      flam    - flux array (1D)
      [eflam] - flux error array (optional; 1D) 
      [polyOrder]  - polynomial order to fit (default value set)
      [indLib]     - index librray of absorption features to omit when fitting polynomial (optional)
      [returnPoly] - return the polynomial as well as normed flux (and normed flux error if passed eflam)

    Returns:
      nflam     - continuum normalised flux array (1D)
      [neflam]  - continuum normalised flux error array (if eflam given; 1D)
      [polyFit] - polynomial fitted to continuum (if returnPoly set; spectrum, not coefs)
    
    """
    
    # scale wave to be between -1 and 1
    x = 2.0*(lam-lam.min())/(lam.max()-lam.min()) - 1.0

    nflam  = []
    pfits  = []

    # loop over spectra, normalising

    # init
    xfit = np.copy(x)
    wfit = np.copy(lam)
    yfit = np.copy(flam)
    if eflam is not None: efit = np.copy(eflam)

    # remove nans and err<0.0
    if eflam is not None:
        loc=np.where((~np.isfinite(flam)) | (~np.isfinite(eflam)) | (eflam<=0.0))
    else:
        loc=np.where(~np.isfinite(flam))  
    good = np.ones_like(lam, dtype=np.bool)
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
    if polyOrder>=0:
        if eflam is not None:
            coef = cheby.chebfit(xfit, yfit, deg=polyOrder, w=1.0/efit**2.0) # variance weight poly fit
        else:
            coef = cheby.chebfit(xfit, yfit, deg=polyOrder)
        polyFit = cheby.chebval(x,coef)
    elif polyOrder==-1:
        # use median - for use in stacking, to maintain continuum shape
        polyFit = np.ones_like(x) * np.median(yfit)
    else:
        raise ValueError, "polyOrder not understood: -1 (median), or >=0 for true polynomial fit"
        
    # normalise 
    nflam=flam/polyFit

    # return extra arrays
    if eflam is not None:
        neflam = eflam/polyFit

    rlist = [nflam]
    if eflam is not None: rlist.append(neflam)
    if returnPoly: rlist.append(polyFit)

    return rlist
        
