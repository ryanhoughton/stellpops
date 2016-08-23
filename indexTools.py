import pylab as pl
import numpy as np
import atpy as ap
from stellpops.specTools import air2vac
import pdb
import copy

# define latex pretty printing labels for indices
prettyPrint = {'CN_1':r'CN$_1$', \
               'CN_2':r'CN$_2$', \
               'Ca4227':r'Ca$_{4227}$', \
               'G4300':r'G4300', \
               'Fe4383':r'Fe$_{4383}$', \
               'Ca4455':r'Ca$_4455$', \
               'Fe4531':r'Fe$_4531$', \
               'Fe4668':r'Fe$_{4668}$', \
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

class ind(dict):
    """
    RH 9/6/16

    The core class for an index object

    """

    def __init__(self, ind_start=None, ind_stop=None, blue_start=None, blue_stop=None, red_start=None, red_stop=None, \
                 name=None, verbose=False):
        """
        Initalise class

        NOTE THAT CONT_START AND CONT_STOP ARE DEFINED INTERNALLY AND USED TO INITALISE (used for index calc)
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
        ID = {'name':name}
        t = {'simpleIndex': True}
        
        self.update(ind)
        self.update(cont)
        self.update(ncont)
        self.update(nfeat)
        self.update(ID)
        self.update(t)
        if verbose: print "Added index "+name

    ## allow self.var instead of self['var']
    #def __getattribute__(self,item):
    #    return self[item]

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
                       addLabel=True, justLine=True, linePos=0.9, labelPos=0.95, labelSize=20.0):
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
                pl.text(np.mean([start,stop])*scale, labelPos*ymax, self['name'], \
                        horizontalalignment='center', fontsize=labelSize)
                pl.plot([start*scale,stop*scale], [ymax*labelPos*0.95, ymax*labelPos*0.95], "k-")



class indlib():
    """
    RH 14/3/16

    The index library class, using dictionaries for both index names, and then sub dictionaries
    for index bands (blue, red, feature)
    
    """

    def __init__(self, name=None, ind_start=None, ind_stop=None, blue_start=None, \
                 blue_stop=None, red_start=None, red_stop=None, table=None, dicts=None, verbose=False):
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
            self.add_simple_indices_via_table(verbose=verbose)
        elif dicts is not None:
            if verbose: print "Adding index elements via list of dictionaries"
            self.add_simple_index_via_dicts(dicts, verbose=verbose)
        else:
            if verbose: print "Adding index element manually via individual elements"
            self.add_simple_index(name=name,ind_start=ind_start, ind_stop=ind_stop, \
                                  blue_start=blue_start, blue_stop=blue_stop, \
                                  red_start=red_start, red_stop=red_stop, verbose=verbose)

    # allow self.var instead of self['var']
    def __getitem__(self,item):
        return getattr(self,item)

    def add_simple_index(self, name=None, ind_start=None, ind_stop=None, blue_start=None, \
                  blue_stop=None, red_start=None, red_stop=None, verbose=False):
        """
        Add a single index 'manually'

        """
        # sanity check
        assert np.any(map(lambda a: a is not None, locals())), "One of the required inputs was not defined"

        newind = ind(ind_start=ind_start, ind_stop=ind_stop, blue_start=blue_start, \
                     blue_stop=blue_stop, red_start=red_start, red_stop=red_stop, name=name)
        setattr(self,name,newind)
        self.nindex+=1
        self.names.append(name)
        if verbose: print "Added index "+name, all

    def augment_index(self, name=None, ind_start=None, ind_stop=None, blue_start=None, \
                  blue_stop=None, red_start=None, red_stop=None, verbose=False):
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

    def add_simple_indices_via_table(self, verbose=False):
        """
        Add multiple indices through a table
        """

        nind = len(self.table.name)
        current_index=None
        for ni in xrange(nind):
            if self.table.name[ni] != '-':
                self.add_simple_index(name=self.table.name[ni], ind_start=self.table.ind_start[ni], ind_stop=self.table.ind_stop[ni], \
                                      blue_start=self.table.blue_start[ni], blue_stop=self.table.blue_stop[ni], \
                                      red_start=self.table.red_start[ni], red_stop=self.table.red_stop[ni], verbose=verbose)
                current_index = self.table.name[ni]
            else:
                assert current_index is not None, "No index preceeding index with name=='-'"
                self.augment_index(name=current_index, ind_start=self.table.ind_start[ni], ind_stop=self.table.ind_stop[ni], \
                                   blue_start=self.table.blue_start[ni], blue_stop=self.table.blue_stop[ni], \
                                   red_start=self.table.red_start[ni], red_stop=self.table.red_stop[ni], verbose=verbose)

    def add_simple_index_via_dict(self, index_dict, verbose=False):
        """
        Add a single index via (unspecified) dict class

        This is likely surplus 
        """

        # make sure dict has all the required attributes, and they're not None
        for cb in self.core_index_tags:
            assert cb in index_dict, "Failed to find attribute "+cb+" in dictionary"
            assert index_dict.get(cb) is not None, "Dictionary has None for attribute "+cb
        setattr(self, index_dict.get('name'), index_dict)
        self.names.append(index_dict.get('name'))
        
    def add_simple_index_via_dicts(self, dicts, verbose=False):
        """
        By providing a list of dicts, add indicies: wrapper for add_simple_index_via_dict
        """
        for d in dicts:
            self.add_simple_index_via_dict(d, verbose=verbose)

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
                       labelPos=[0.8,0.85,0.9,0.95], \
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
        for n, name in enumerate(indexList):
            if np.any(self[name]['blue_start']*scale>xmin) & np.any(self[name]['red_stop']*scale<xmax):
                labpos = labelPos[n % nlabel]
                linpos = linePos[n % nlabel]
                self[name].plotIndexBands(scale=scale, alpha=alpha, contCol='blue', indCol='red', \
                                          ymin=ymin, ymax=ymax, autoy=False, justLine=justLine, \
                                          addLabel=addLabel, labelPos=labpos, linePos=linpos, \
                                          labelSize=labelSize, noContinuua=noContinuua)
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
        
        
        

def loadLickIndicesAir(filename="/home/houghton/z/data/stellar_pops/lickIndicesAir.txt",verbose=False):
    """
    Load the Lick indices from Worthey's website (http://astro.wsu.edu/worthey/html/index.table.html, air wavelengths). These are a compilation from Trager et al. (1998) and Worthey & Ottaviani (1997). They are the same as used in TMJ10 models. 

    """

    tab = ap.Table(filename, type='ascii')
    inds = indlib(table=tab, verbose=verbose)
    return inds

def loadLickIndicesVac(filename="/home/houghton/z/data/stellar_pops/lickIndicesAir.txt",verbose=False):
    """
    Like above but convert to vacuum wavelengths
    """
    table = ap.Table(filename, type='ascii')
    nline, ncol = table.shape
    for nli in xrange(nline):
        for nci in xrange(1,ncol-2):
            table[nli][nci] = air2vac(table[nli][nci],verbose=verbose)

    inds = indlib(table=table, verbose=verbose)
    return inds

def getLickIndices(verbose=False):
    tab = loadLickIndicesVac()
    return indlib(table=tab,verbose=verbose)

def loadCvD12IndicesVac(filename="/home/houghton/z/data/stellar_pops/CvDIndicesVac.txt"):
    """
    Load the indicies presented in Table 1 of Conroy & van Dokkum 2012a

    Vacuum wavelengths

    """
    tab = ap.Table(filename, type='ascii')
    return tab

def getCvD12Indices(verbose=False):
    return indlib(table=loadCvD12IndicesVac(),verbose=verbose)


def calcMeanFe(Fe5270, Fe5335, eFe5270=None, eFe5335=None):
    meanFe = 0.5*(Fe5270+Fe5335)
    if (eFe5270 is not None) & (eFe5335 is not None):
        emeanFe = 0.5*(eFe5270+eFe5335)
        rlist = [meanFe, emeanFe]
    else:
        rlist = meanFe
    return rlist

def calcMgFe(Mgb, Fe5270, Fe5335, eMgb=None, eFe5270=None, eFe5335=None):
    """
    Return [MgFe] from Gonzalez 1993

    And error if values passed
    """
    avFe = calcMeanFe(Fe5270, Fe5335)
    a=0.5
    b=0.5
    MgFe = np.sqrt(Mgb * avFe)
    if (eMgb is not None) & (eFe5270 is not None) & (eFe5335 is not None):
        eMgFe = 0.5 / (MgFe) * (eMgb/2.0*avFe + Mgb/2.0*(a*eFe5270+b*eFe5335))
        rlist = [MgFe, eMgFe]
    else:
        rlist = MgFe
        
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
        rlist = [MgFeP, eMgFeP]
    else:
        rlist = MgFeP

    return rlist



