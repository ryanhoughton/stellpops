import pylab as pl
import numpy as np
import atpy as ap
from stellpops import specTools as st
import pdb
import copy
import warnings as warn


# define latex pretty printing labels for indices
prettyPrint = {'CN_1':r'CN$_1$', \
               'CN_2':r'CN$_2$', \
               'Ca4227':r'Ca$_{4227}$', \
               'G4300':r'G4300', \
               'Fe4383':r'Fe$_{4383}$', \
               'Ca4455':r'Ca$_{4455}$', \
               'Fe4531':r'Fe$_{4531}$', \
               'Fe4668':r'C2$_{4668}$', \
               'H_beta':r'H$\beta$', \
               'Fe5015':r'Fe$_{5015}$', \
               'Mg_1':r'Mg$_1$', \
               'Mg_2':r'Mg$_2$', \
               'Mg_b':r'Mg$_b$', \
               'Fe5270':r'Fe$_{5270}$',\
               'Fe5335':r'Fe$_{5335}$',\
               'Fe5406':r'Fe$_{5406}$', \
               'Fe5782':r'Fe$_{5782}$', \
               'Na_D':r'Na$_{D}$',\
               'TiO_1':r'TiO$_1$',\
               'TiO_2':r'TiO$_2$',\
               'Hdelta_A':r'H$\delta_A$', \
               'Hdelta_F':r'H$\delta_F$', \
               'Hgamma_A':r'H$\gamma_A$', \
               'Hgamma_F':r'H$\gamma_F$'}


# from Worthey and Ottaviani, 1997, Table 8, Wavelength in AA and FWHM in AA.
lickResolutions=np.array([[4000,4400,4900,5400,6000],[11.5, 9.2, 8.4, 8.4, 9.8]])
lickResPolyFit=[1.23827561e-12,  -2.53065476e-08,   1.94910696e-04, -6.70562662e-01,   8.77800000e+02]

class ind(dict):
    """
    RH 9/6/16

    The core class for an index object

    """

    def __init__(self, ind_start=None, ind_stop=None, blue_start=None, blue_stop=None, red_start=None, red_stop=None, \
                 resol=None, name=None, wavesyst=None, verbose=False):
        """
        Initalise index class, used to calculate equivalent widths using two (blue and red) continuum bandpasses
        and a (index) feature bandpass.

        The actual calculation of the index is left to a different tool (see specTools.py)

        Inputs:
        - ind_start : defines the (blue) start of the feature bandpass
        - ind_stop  : defines the (red) end of the feature bandpass
        - blue_start: defines the (blue) start of the blue continuum bandpass
        - blue_stop : defines the (red) end of the blue continuum bandpass
        - red_start : defines the (blue) start of the red continuum bandpass
        - resol     : the resolution (FWHM) in AA that the index should be measured at. This is normally set to None
                      (resolution of spectrum is not changed) unless the index is to be calculated on the Lick/IDS system.
        - name      : the name of the index
        - verbose   : print stuff to screen, which can get annoying.
    

        NOTE THAT CONT_START AND CONT_STOP ARE DEFINED INTERNALLY (used for index calc). 
        
        """

        assert np.any(map(lambda a: a is not None, locals())), "One of the required inputs was not defined"

        cont_start=[]
        cont_start.extend(np.atleast_1d(blue_start).tolist())
        cont_start.extend(np.atleast_1d(red_start).tolist())
        cont_stop=[]
        cont_stop.extend(np.atleast_1d(blue_stop).tolist())
        cont_stop.extend(np.atleast_1d(red_stop).tolist())
        ind = {'ind_start':np.atleast_1d(ind_start).tolist(), 'ind_stop':np.atleast_1d(ind_stop).tolist()}
        cont = {'cont_start':cont_start, 'cont_stop':cont_stop, \
                'blue_start':blue_start, 'blue_stop':blue_stop, 'red_start':red_start, 'red_stop':red_stop}
        nfeat = {'nfeat':1}
        ncont = {'ncont':2}
        resol = {'resol':resol}
        ID = {'name':name}
        t = {'simpleIndex': True}
        
        
        self.update(ind)
        self.update(cont)
        self.update(ncont)
        self.update(nfeat)
        self.update(resol)
        self.update(ID)
        self.update(t)

        if wavesyst is not None:
            if (wavesyst=="vac" or wavesyst=="VAC" or wavesyst=="Vac"):
                self.update({"wavesyst":"vac"})
            elif (wavesyst=="air" or wavesyst=="AIR" or wavesyst=="Air"):
                self.update({"wavesyst":"air"})
            else:
                raise ValueError("wavesyst not understood. Should be air or vac.")
        else:
            warn.warn("You failed to specify if the wavelengths are defined in AIR or VAC units.")
            self.wavesyst=None
            
        if verbose: print "Added index "+name

    ## allow self.var instead of self['var']
    #def __getattribute__(self,item):
    #    return self[item]

    def copy(self):
        """
        RH 27/10/2016
        Return a copy of the index (to allow edits without changing original)
        """
        newInd = ind(ind_start=self['ind_start'], ind_stop=self['ind_stop'], \
                     blue_start=self['blue_start'], blue_stop=self['blue_stop'], \
                     red_start=self['red_start'], red_stop=self['red_stop'], \
                     resol=self['resol'], name=self['name'])
        return newInd

    def augmentIndex(self,ind_start=None, ind_stop=None, \
                     red_start=None, red_stop=None, blue_start=None, blue_stop=None):
        """
        
        Add additional features to make a complex index

        MUST SPECIFY RED AND BLUE CONTINUUA SEPARATELY
        
        """
        # add new features
        if (type(ind_start)!=type(None)) & (type(ind_stop)!=type(None)):
            assert len(ind_start)==len(ind_stop), 'Non-equal number of start and stop positions for new feature'
            self['ind_start'].extend([ind_start])
            self['ind_stop'].extend([ind_stop])
            self['nfeat']+=1
        # add new blue cont
        if (type(blue_start)!=type(None)) & (type(blue_stop)!=type(None)):
            assert len(blue_start)==len(blue_stop), 'Non-equal number of start/stop positions for new blue continuum'
            self['blue_start'].extend([blue_start])
            self['blue_stop'].extend([blue_stop])
            self['cont_start'].extend([blue_start])
            self['cont_stop'].extend([blue_stop])
            self['ncont']+=1
        # add new red cont
        if (type(red_start)!=type(None)) & (type(red_stop)!=type(None)):
            assert len(red_start)==len(red_stop), 'Non-equal number of start/stop positions for new red continuum'
            self['red_start'].extend([red_start])
            self['red_stop'].extend([red_stop])
            self['cont_start'].extend([red_start])
            self['cont_stop'].extend([red_stop])
            self['ncont']+=1
        # label as NOT SIMPLE
        self['simpleIndex']=False

    def redshiftIndex(self, redshift, updateName=False):
        """
        Redshift an index's bandpasses, and return new index in the observer frame.
        """

        if updateName:
            name = self['name']+"_z"+str(redshift)
        else:
            name = self['name']
        newind=ind(name=name, \
                   red_start=list(np.atleast_1d(self['red_start'])*(1.0+redshift)), \
                   red_stop=list(np.atleast_1d(self['red_stop'])*(1.0+redshift)), \
                   blue_start=list(np.atleast_1d(self['blue_start'])*(1.0+redshift)), \
                   blue_stop=list(np.atleast_1d(self['blue_stop'])*(1.0+redshift)), \
                   ind_start=list(np.atleast_1d(self['ind_start'])*(1.0+redshift)), \
                   ind_stop=list(np.atleast_1d(self['ind_stop'])*(1.0+redshift)))
        return newind

    def plotIndexBands(self, scale=1.0, alpha=0.1, contCol='blue', indCol='red', \
                       ymin=-100.0, ymax=100.0, autoy=False, noContinuua=False, \
                       addLabel=True, justLine=True, linePos=0.9, labelPos=0.95, \
                       labelSize=20.0, showHorizLine=False, usePrettyPrint=True):
        """
        Overplot the regions used for continuum and feature bandpasses
        
        """

        if autoy:
            gca=pl.gca()
            ymin, ymax = gca.get_ylim()

        if justLine:
            for nf in xrange(self['nfeat']):
                for istart, istop in zip(self['ind_start'], self['ind_stop']):
                    meanPos = 0.5*(istart+istop)
                    pl.plot([meanPos*scale, meanPos*scale], [ymin,ymax*linePos], "k:")
        else:
            for nf in xrange(self['nfeat']):
                pl.fill_between([self['ind_start'][nf]*scale, self['ind_stop'][nf]*scale], \
                                [ymin,ymin], [ymax,ymax], color=indCol, alpha=alpha)
            if not noContinuua:
                for nc in xrange(self['ncont']):
                    pl.fill_between([self['cont_start'][nc]*scale, self['cont_stop'][nc]*scale], \
                                    [ymin,ymin], [ymax,ymax], color=contCol, alpha=alpha)
        if addLabel:
            for start,stop in zip(self['ind_start'], self['ind_stop']):
                if usePrettyPrint:
                    label=prettyPrint[self['name']]
                else:
                    label=self['name']
                xpos = np.mean([start,stop])*scale
                ypos = labelPos*ymax
                pl.text(xpos, ypos, label, \
                        horizontalalignment='center', fontsize=labelSize)
                if showHorizLine:
                    pl.plot([start*scale,stop*scale], [ymax*labelPos*0.95, ymax*labelPos*0.95], "k-")



class indlib():
    """
    RH 14/3/16

    The index library class, using dictionaries for both index names, and then sub dictionaries
    for index bands (blue, red, feature)
    
    """

    def __init__(self, name=None, ind_start=None, ind_stop=None, blue_start=None, \
                 blue_stop=None, red_start=None, red_stop=None, resol=None, table=None, \
                 dicts=None, wavesyst=None, verbose=False):
        """
        Set up the index class in two ways depending on input:
        - do single index with individual components specified
        - do multiple indices using a table

        """

        self.nindex=0
        self.core_index_tags=['name', 'ind_start', 'ind_stop', 'cont_start', 'cont_stop', 'ncont', 'nfeat']
        self.names=[]
        
        if table is not None:
            if verbose: print "Adding index elements via table"
            self.table = table
            self.add_simple_indices_via_table(wavesyst=wavesyst, verbose=verbose)
        elif dicts is not None:
            if verbose: print "Adding index elements via list of dictionaries"
            self.add_simple_index_via_dicts(dicts, wavesyst=wavesyst, verbose=verbose)
        else:
            if verbose: print "Adding index element manually via individual elements"
            self.add_simple_index(name=name,ind_start=ind_start, ind_stop=ind_stop, \
                                  blue_start=blue_start, blue_stop=blue_stop, \
                                  red_start=red_start, red_stop=red_stop, wavesyst=wavesyst, \
                                  verbose=verbose)

    # allow self.var instead of self['var']
    def __getitem__(self,item):
        return getattr(self,item)

    def add_simple_index(self, name=None, ind_start=None, ind_stop=None, blue_start=None, \
                  blue_stop=None, red_start=None, red_stop=None, resol=None, wavesyst=None, verbose=False):
        """
        Add a single index 'manually'

        """
        # sanity check
        assert np.any(map(lambda a: a is not None, locals())), "One of the required inputs was not defined"

        newind = ind(ind_start=ind_start, ind_stop=ind_stop, blue_start=blue_start, \
                     blue_stop=blue_stop, red_start=red_start, red_stop=red_stop, resol=resol, \
                     name=name, wavesyst=wavesyst)
        setattr(self,name,newind)
        self.nindex+=1
        self.names.append(name)
        if verbose: print "Added index "+name, all

    def augment_index(self, name=None, ind_start=None, ind_stop=None, blue_start=None, \
                  blue_stop=None, red_start=None, red_stop=None, wavesyst=None, verbose=False):
        """
        Add aditional continuua / features to an existing index
        """

        # sanity checks
        assert np.any(map(lambda a: a is not None, locals())), "One of the required inputs was not defined"
        assert hasattr(self, name), "This index is not defined, so we cannot augment it"

        #ind = getattr(self, name)
        #ind['ind_start'].extend([ind_start])
        #ind['ind_stop'].extend([ind_stop])
        #ind['cont_start'].extend([ind_start])
        #ind['cont_stop'].extend([ind_stop])
        #ind['nfeat']+=1
        setattr(self,name,ind)
        if verbose: print "Augmented index "+name+" to ", ind

    def add_simple_indices_via_table(self, wavesyst=None, verbose=False):
        """
        Add multiple indices through a table
                
        """

        # ERROR this doesn't work as None has no length, breaking atpy
        # add resol if not in table
        #if self.table.keys().count('resol')==0:
        #    resol = [None for x in xrange(len(self.table))]
        #    self.table.add_column('resol', resol)

        nind = len(self.table.name)
        current_index=None
        for ni in xrange(nind):
            if self.table.name[ni] != '-':
                self.add_simple_index(name=self.table.name[ni], ind_start=self.table.ind_start[ni], ind_stop=self.table.ind_stop[ni], \
                                      blue_start=self.table.blue_start[ni], blue_stop=self.table.blue_stop[ni], \
                                      red_start=self.table.red_start[ni], red_stop=self.table.red_stop[ni], \
                                      wavesyst=wavesyst, verbose=verbose) # resol=self.table.resol[ni], 
                current_index = self.table.name[ni]
            else:
                assert current_index is not None, "No index preceeding index with name=='-'"
                self.augment_index(name=current_index, ind_start=self.table.ind_start[ni], ind_stop=self.table.ind_stop[ni], \
                                   blue_start=self.table.blue_start[ni], blue_stop=self.table.blue_stop[ni], \
                                   red_start=self.table.red_start[ni], red_stop=self.table.red_stop[ni], \
                                   wavesyst=wavesyst, verbose=verbose) # resol=self.table.resol[ni], 

    def add_simple_index_via_dict(self, index_dict, wavesyst=None, verbose=False):
        """
        Add a single index via (unspecified) dict class

        This is likely surplus
        Assume no resolution in dict
        """

        # make sure dict has all the required attributes, and they're not None
        for cb in self.core_index_tags:
            assert cb in index_dict, "Failed to find attribute "+cb+" in dictionary"
            assert index_dict.get(cb) is not None, "Dictionary has None for attribute "+cb
            
        setattr(self, index_dict.get('name'), index_dict)
        self[index_dict.get('name')]['wavesyst']=wavesyst
        self.names.append(index_dict.get('name'))
        
    def add_simple_index_via_dicts(self, dicts, wavesyst=None, verbose=False):
        """
        By providing a list of dicts, add indicies: wrapper for add_simple_index_via_dict
        """
        for d in dicts:
            self.add_simple_index_via_dict(d, wavesyst=wavesyst, verbose=verbose)

    def keys(self):
        """
        Return the index key words (ind names and extra stuff about lib)

        DO NOT USE THIS TO GET INDEX LIST - USE SELF.NAMES
        """
        return self.__dict__.keys()

    def d(self):
        """
        Return class as dictionary
        """
        return self.__dict__

    def redshiftLib(self, redshift, overwrite=False, augment=False):
        """
        Redshift the index libary bandpasses, and return new lib.
        Wrapper for 'redshiftIndex'
        
        """

        newinds=[]
        # loop over indices, calcing new shifted ones
        names = self.names
        for name in names:
            oldind = getattr(self,name)
            newind = oldind.redshiftIndex(redshift)
            if overwrite:
                newind['name']=oldind['name']
                setattr(self,name,newind)
            else:
                newinds.append(newind)
        # return new lib if not overwriting
        if not overwrite:
            newlib = indlib(dicts=newinds)
            return newlib

    def plotIndexBands(self, scale=1.0, alpha=0.1, contCol='blue', indCol='red', \
                       xmin=0.0, xmax=1e9, ymin=-100.0, ymax=100.0, autoxy=False, noContinuua=False, \
                       addLabel=False, justLine=True, linePos=[0.75,0.8,0.85,0.9], \
                       labelPos=[0.8,0.85,0.9,0.95], stagger=0.1, nstagger=2, \
                       labelSize=20.0, indexList=None):
        """        Overplot the index bands, that are in the current plot range
        
        """

        if indexList==None:
            indexList = self.names

        if autoxy:
            gca = pl.gca()
            xmin, xmax = gca.get_xlim()
            ymin, ymax = gca.get_ylim()

        nlabel = len(labelPos)
        nline = len(linePos)
        assert nlabel==nline, "#labels!=#lines"
        nindex=0
        for n, name in enumerate(indexList):
            if np.any(self[name]['blue_start']*scale>xmin) & np.any(self[name]['red_stop']*scale<xmax):
                labpos = labelPos[nindex % nlabel]
                linpos = linePos[nindex % nlabel]
                self[name].plotIndexBands(scale=scale, alpha=alpha, contCol='blue', indCol='red', \
                                          ymin=ymin, ymax=ymax, autoy=False, justLine=justLine, \
                                          addLabel=addLabel, labelPos=labpos, linePos=linpos, \
                                          labelSize=labelSize, noContinuua=noContinuua)
                nindex+=1
    def subset(self, indexList):
        """
        RH 22/8/2016

        Make a new index library based on the list given; essentially a subset of this library.

        - IndexList: a list of strings referring to specific keys in this library
        
        """

        listOfIndices = []
        for ind in indexList:
            listOfIndices.append(self[ind])
        newLib=indlib(dicts=listOfIndices)
        return newLib
        

def loadLickIndicesAir(filename="/home/houghton/z/data/stellar_pops/lickIndicesAir.txt", atLickRes=False, verbose=False):
    """
    Load the Lick indices from Worthey's website (http://astro.wsu.edu/worthey/html/index.table.html, air wavelengths).
    These are a compilation from Trager et al. (1998) and Worthey & Ottaviani (1997). They are the same as used in TMJ10 models. 

    """

    tab = ap.Table(filename, type='ascii')
    inds = indlib(table=tab, wavesyst="air", verbose=verbose)
    if atLickRes: inds = addLickRes2IndLib(inds)
    return inds

def loadLickIndicesVac(filename="/home/houghton/z/data/stellar_pops/lickIndicesAir.txt", atLickRes=False, verbose=False):
    """
    Like above but convert to vacuum wavelengths
    """
    table = ap.Table(filename, type='ascii')
    nline, ncol = table.shape
    for nli in xrange(nline):
        for nci in xrange(1,ncol-2):
            table[nli][nci] = st.air2vac(table[nli][nci],verbose=verbose)
    inds = indlib(table=table, wavesyst="vac", verbose=verbose)
    # inset Lick resolutions - probably should be done at air wavelengths but surely tiny tiny error
    if atLickRes: inds = addLickRes2IndLib(inds)
    return inds

def getLickIndices(verbose=False):
    tab = loadLickIndicesVac()
    return indlib(table=tab,verbose=verbose)

def loadCvD12IndicesVac(filename="/home/houghton/z/data/stellar_pops/CvDIndicesVac.txt", verbose=False):
    """
    Load the indicies presented in Table 1 of Conroy & van Dokkum 2012a into indLib

    Vacuum wavelengths

    """
    tab = ap.Table(filename, type='ascii')
    lib = indlib(table=tab, wavesyst="vac", verbose=verbose)
    return lib


def calcMeanFe(Fe1, Fe2, Fe3=None, eFe1=None, eFe2=None, eFe3=None):
    """
    Wrapper for choosing between mean2Fe and mean3Fe
    """
    if Fe3 is None:
        res = calcMean2Fe(Fe1, Fe2, eFe5270=eFe1, eFe5335=eFe2)
    else:
        res = calcMean3Fe(Fe1, Fe2, Fe3, eFe5270=eFe1, eFe5335=eFe2, eFe5406=eFe3)
    return res

def calcMean2Fe(Fe5270, Fe5335, eFe5270=None, eFe5335=None):
    meanFe = 0.5*(Fe5270+Fe5335)
    if (eFe5270 is not None) & (eFe5335 is not None):
        emeanFe = np.sqrt(0.25*eFe5270**2.0 + 0.25*eFe5335**2.0) 
        rlist = [meanFe, emeanFe]
    else:
        rlist = meanFe
    return rlist

def calcMean3Fe(Fe5270, Fe5335, Fe5406, eFe5270=None, eFe5335=None, eFe5406=None):
    mean3Fe = (Fe5270+Fe5335+Fe5406)/3.
    if (eFe5270 is not None) & (eFe5335 is not None) & (eFe5406 is not None):
        emean3Fe = np.sqrt((1./6.)*eFe5270**2.0 + (1./6.)*eFe5335**2.0 + (1./6.)*eFe5406**2.0) 
        rlist = [mean3Fe, emean3Fe]
    else:
        rlist = mean3Fe
    return rlist


def calcMgFe(Mgb, Fe1, Fe2, Fe3=None, eMgb=None, eFe1=None, eFe2=None, eFe3=None):
    """
    Wrapper for choosing between MgFe2 and MgFe3
    """
    if Fe3 is None:
        res = calcMgFe2(Mgb, Fe1, Fe2, eMgb=eMgb, eFe5270=eFe1, eFe5335=eFe2)
    else:
        res = calcMgFe3(Mgb, Fe1, Fe2, Fe3, eMgb=eMgb, eFe5270=eFe1, eFe5335=eFe2, eFe5406=eFe3)
    return res

def calcMgFe2(Mgb, Fe5270, Fe5335, eMgb=None, eFe5270=None, eFe5335=None):
    """
    Return [MgFe] from Gonzalez 1993

    And error if values passed
    """
    if (eFe5270 is not None) & (eFe5335 is not None):
        avFe, eavFe = calcMean2Fe(Fe5270, Fe5335, eFe5270=eFe5270, eFe5335=eFe5335)
    else:
        avFe = calcMean2Fe(Fe5270, Fe5335)
    a=0.5
    b=0.5
    MgFe = np.sqrt(Mgb * avFe)
    if (eMgb is not None) & (eFe5270 is not None) & (eFe5335 is not None):
        # see RH note book 28/9/16
        eMgFe = (np.sqrt(0.25*(eMgb/Mgb)**2.0 + 0.25*(eavFe/avFe)**2.0)) * MgFe # relative error calc, scaled up by value
        rlist = [MgFe, eMgFe]
    else:
        rlist = MgFe
        
    return rlist

def calcMgFe3(Mgb, Fe5270, Fe5335, Fe5406, eMgb=None, eFe5270=None, eFe5335=None, eFe5406=None):
    """
    Return modified [MgFe] akin to  Gonzalez 1993

    And error if values passed
    """
    if (eFe5270 is not None) & (eFe5335 is not None) & (eFe5406 is not None):
        av3Fe, eav3Fe = calcMean3Fe(Fe5270, Fe5335, Fe5406, eFe5270=eFe5270, eFe5335=eFe5335, eFe5406=eFe5406)
    else:
        av3Fe = calcMean3Fe(Fe5270, Fe5335, Fe5406)
    a=0.5
    b=0.5
    MgFe3 = np.sqrt(Mgb * av3Fe)
    if (eMgb is not None) & (eFe5270 is not None) & (eFe5335 is not None) & (eFe5406 is not None):
        # see RH note book 28/9/16
        eMgFe3 = (np.sqrt(0.25*(eMgb/Mgb)**2.0 + 0.25*(eavFe/avFe)**2.0)) * MgFe # relative error calc, scaled up by value
        rlist = [MgFe3, eMgFe3]
    else:
        rlist = MgFe3
        
    return rlist

def calcMgFePrime(Mgb, Fe5270, Fe5335, eMgb=None, eFe5270=None, eFe5335=None):
    """
    Return the [MgFe]' index from TMB03

    And error if values passed
    """
    a = 0.72
    b = 0.28
    avFe = (a*Fe5270 + b*Fe5335)
    MgFeP = np.sign(Mgb) * np.sign(avFe) * np.sqrt(np.fabs(Mgb * avFe))

    if (eMgb is not None) & (eFe5270 is not None) & (eFe5335 is not None):
        eMgFeP = 0.5 / (MgFeP) * (eMgb/2.0*avFe + Mgb/2.0*(a*eFe5270+b*eFe5335))
        raise ValueError("Error propagation needs correcting, see above and note book 28/9/16")
        rlist = [MgFeP, eMgFeP]
    else:
        rlist = MgFeP

    return rlist


def fitLickRes(data=lickResolutions, order=4, npoints=100, doPlot=True, update=True):
    """
    RH 20/10/2016

    Fit the lick resolution as a function of wavelength with a quadratic, using the data from Worthey and Ottaviani, 1997, Table 8.
    Then update the module global lickResPolyFit, if asked

    With order=4, I get: [  1.23827561e-12,  -2.53065476e-08,   1.94910696e-04, -6.70562662e-01,   8.77800000e+02]

    """
    global lickResPolyFit
    xx=np.linspace(data[0].min()*0.9, data[0].max()*1.1, npoints)
    coef = np.polyfit(data[0], data[1], order)
    fit = np.polyval(coef, xx)

    if doPlot:
        pl.plot(data[0], data[1], "ko", label="Data")
        pl.plot(xx, fit, "r-", label="Polynomial fit")

    if update:
        lickResPolyFit=coef

    return coef

def addLickRes2IndLib(lib, coefs=lickResPolyFit):
    """
    RH 20/10/16

    Cycle through index library, updating the resolution keyword to that of the Lick resolution.

    """
    indices=lib.names
    for ind in indices:
        # calc central feature wavelength
        clam = np.mean(np.array(lib[ind]['ind_start'] + lib[ind]['ind_stop']))
        FWHM = np.polyval(lickResPolyFit, clam)
        lib[ind]['resol']=FWHM

    return lib


def calcVelDispIndexCorrection(specs, indexIN, maxVelDisp=500.0, minVelDisp=None, normVelDisp=None, nSample=20, order=4, \
                               ageMinMax=None, returnAll=False, doPlot=True, verbose=True):
    """
    RH 21/10/2016

    Absorption line indices change with the velocity dispersion of the spectra. This routine
    calculates how a given index changes with increasing galaxy dispersion.

    Note that the spectral resolution of the spectra limit how low in velocity dispersion we can go. 

    e.g.

    s2 = it.calcVelDispIndexCorrection(ms[2], lib['H_beta'], 500.0, None, normVelDisp=200.0, nSample=21, ageMinMax=[1,15])
    
    """

    # make a local copy of the index
    index = indexIN.copy()
    index['resol']=None # this is not allowed for main index calc loop

    # get central wavelength for index
    meanLam = np.mean(np.array([index['ind_start'], index['ind_stop']]))

    # check that all spectra on: 1) same resolution, 2) same age grid
    # get resolution of spectrum - FWHM AA
    specResAA = specs[0].calcResolution(meanLam)
    age = specs[0].age
    for s in specs:
        assert specResAA==s.calcResolution(meanLam), "Cannot make use of spectra with different resolutions"
        assert np.all(s.age==age), "Models on different age grids."
    # get effective dispersion
    sigmaSpec = (specResAA / meanLam) / np.sqrt(8.*np.log(2.)) * st.c / 1e3

    # define minVelDisp
    if minVelDisp is None: minVelDisp=sigmaSpec
    # sanity check
    assert maxVelDisp > sigmaSpec, "Cannot convolve to the max velocity dispersion, library resolution too coarse"
    assert minVelDisp >= sigmaSpec, "Cannot convolve to the min velocity dispersion, library resolution too coarse" 

    # make array of disp values
    outVelDisps = np.linspace(minVelDisp, maxVelDisp, nSample) 

    # now cycle through each spectrum
    allIndices = []
    normAllIndices = []
    for spec in specs:
        
        # now loop over disps
        allInds=[]
        for n, sig in enumerate(outVelDisps):
            # cut and convolve
            if sig==sigmaSpec:
                cspec = spec
            else:
                cspec = st.cutAndGaussVelConvolve(spec, index, sig, verbose=False)
                if verbose: print "Done convolution for "+str(n)+" of "+str(nSample)+" values of sigma"
            # calc index
            inds = cspec.calcIndex(index)
            allInds.append(inds)
        ai = np.array(allInds)
        allIndices.append(ai.T) # convert to array

        # calc the normalising values - do this now to stop resolution errors later
        if (normVelDisp is None) & (indexIN['resol'] is None):
            normLoc = 0
            normVals = ai[normLoc,:]
        elif (normVelDisp is not None) & (indexIN['resol'] is None):
            normLoc = np.where(outVelDisps>=normVelDisp)[0][0] # use nearest value (not ideal)
            actualNormDisp = outVelDisps[normLoc]
            if normVelDisp != actualNormDisp: print "WARNING: adopting "+str(actualNormDisp)+" as normVelDisp"
            normVals = ai[normLoc,:]
        elif (normVelDisp is None) & (indexIN['resol'] is not None):
            # normalise to Lick value
            normVals = spec.calcIndex(indexIN)
        elif (normVelDisp is not None) & (indexIN['resol'] is not None):
            raise ValueError("Both normVelDisp and index.resolution are both NOT None. Either normalise to Lick resolution or normalise to fixed sigma")
        else:
            raise ValueError("How did this happen?")
        
        # cut by ages if asked
        if ageMinMax is not None:
            aloc = np.where( (spec.age >= ageMinMax[0]) & (spec.age <= ageMinMax[1]))[0]
            ai = ai[:,aloc]
            age = spec.age[aloc]
            normVals=normVals[aloc]
        else:
            age = spec.age

        # normalise by values at chosen sigma
        nai = ai/np.tile(normVals,(nSample,1))
        normAllIndices.append(nai.T)
    # make into arrays
    allIndices=np.array(allIndices)
    normAllIndices=np.array(normAllIndices)
    # now reformat array to that Z and age are all on one axis - to average over
    nailoc = normAllIndices.shape
    nAI1D = normAllIndices.reshape((nailoc[0]*nailoc[1],nailoc[2])).T
    
    # calculate the actual correction polynomial
    x = (outVelDisps - outVelDisps.min())/(outVelDisps.max()-outVelDisps.min()) # scale from 0 to 1
    x = (x-0.5)/0.5 # scale from -1 to +1
    mean = np.mean(nAI1D, axis=1)
    std  = np.std(nAI1D, axis=1)
    # make weight array, fixing for 1/0.0 errors
    w = std
    w[np.where(std==0.0)] = std[np.where(std!=0.0)].min()*0.01 # 1% of min error
    w = 1.0/w
    pcoef = np.polynomial.chebyshev.chebfit(x, mean, deg=order, w=w)
    polyfit = np.polynomial.chebyshev.chebval(x, pcoef)
    epcoef = np.polynomial.chebyshev.chebfit(x, std, deg=order,w=None)
    epolyfit = np.polynomial.chebyshev.chebval(x, epcoef)

    corrector = dispCorr(minVelDisp, maxVelDisp, pcoef, eChebyCoefs=epcoef, outputDisp=normVelDisp)

    rlist = corrector

    if returnAll:
        extra = [outVelDisps, mean, std, normAllIndices]
        rlist=[corrector]
        rlist.append(extra)
    
    # plot if asked
    if doPlot:
        pl.figure(figsize=(14,8))
        pl.subplot(121)
        pl.imshow(np.mean(normAllIndices,axis=0).T)
        pl.ylabel("Velocity Dispersion (non-linear units)")
        pl.xlabel("Age (non-linear units)")
        pl.subplot(122)
        pl.plot(outVelDisps, nAI1D, alpha=0.3)
        pl.errorbar(outVelDisps, mean, yerr=std, color="red")
        pl.plot(outVelDisps, polyfit, "k-")
        pl.plot(outVelDisps, polyfit+epolyfit, "k:")
        pl.plot(outVelDisps, polyfit-epolyfit, "k:")
        pl.ylabel(r'Correction $\sigma_x/\sigma_0$')
        pl.xlabel(r'$\sigma_x$')
        fig = pl.gcf()
        fig.canvas.set_window_title(r'Index correction with dispersion for '+index['name'])

    return rlist


class dispCorr():
    """
    RH 24/10/16

    Class to calculate, store and apply dispersion corrections for various indices
    """
    def __init__(self, minDisp, maxDisp, chebyCoefs, eChebyCoefs=None, outputDisp=200.0):
        """
        Initalise the class by assigning memory
        """
        self.minS=minDisp
        self.maxS=maxDisp
        self.rangeS=maxDisp-minDisp
        self.coefs=chebyCoefs
        if eChebyCoefs is not None:
            self.ecoefs=eChebyCoefs
        else:
            self.ecoefs=None
        # done

    def correct(self, inputDisp, inputIndex, giveError=True):
        """
        Calculate the new index value given the input and output dispersions
        """
        assert (inputDisp > self.minS) & (inputDisp < self.maxS), \
               "Input velocity dispersion out of allowed range: "+str(self.minS)+" to "+str(self.maxS)
        #assert (outputDisp > self.minS) & (outputDisp < self.maxS), \
        #       "Output velocity dispersion out of allowed range: "+str(self.minS)+" to "+str(self.maxS)
        # calc scaled x
        x = (inputDisp-self.minS)/self.rangeS
        x = (x-0.5)/0.5
        #xout= (outputDisp-self.minS)/self.rangeS
        #xout= (xout-0.5)/0.5
        yval = np.polynomial.chebyshev.chebval(x, self.coefs)
        #y2 = np.polynomial.chebyshev.chebval(xout, self.coef)
        rlist = inputIndex/yval
        if self.ecoefs is not None:
            eyval = np.polynomial.chebyshev.chebval(x, self.ecoefs)
            rlist=[inputIndex/yval]
            rlist.append(eyval)

        return rlist


