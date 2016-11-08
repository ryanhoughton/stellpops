import pylab as pl
import numpy as np
from astropy.table import Table
from specTools import air2vac, vac2air
import pdb
import copy

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
               'Hgamma_F':r'H$\gamma_F$', \
               'CaII86': r'CaT', \
               'FeH': r'FeH', \
               'NaIsdss': r'NaI$_{SDSS}$'
               }


# from Worthey and Ottaviani, 1997, Table 8, Wavelength in AA and FWHM in AA.
lickResolutions=[[4000,4400,4900,5400,6000],[11.5, 9.2, 8.4, 8.4, 9.8]]

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

        #Added by SPV:
        If the index doesn't have a feature defintion, we assume that ind_start and ind_stop will be set to negative numbers. In this case,
        set nfeat=0.0
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

        #Check if our wavelengths for the index feature are positive. If not, assume they don't exist and set nfeat to 0.
        if ind_start <0.0 or ind_stop<0.0:
            nfeat = {'nfeat':0}
        else:
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
        import pdb; pdb.set_trace()

        #make everything at least a list
        self['blue_start']=np.atleast_1d(self['blue_start']).tolist()
        self['blue_stop']=np.atleast_1d(self['blue_stop']).tolist()
        self['red_start']=np.atleast_1d(self['red_start']).tolist()
        self['red_stop']=np.atleast_1d(self['red_stop']).tolist()

        # add new features
        if (type(ind_start)!=type(None)) & (type(ind_stop)!=type(None)):
            #assert len(ind_start)==len(ind_stop), 'Non-equal number of start and stop positions for new feature'


            self['ind_start'].extend([ind_start])
            self['ind_stop'].extend([ind_stop])

            self['nfeat']+=1
        # add new blue cont
        if (type(blue_start)!=type(None)) & (type(blue_stop)!=type(None)):
            #assert len(blue_start)==len(blue_stop), 'Non-equal number of start/stop positions for new blue continuum'

            self['blue_start'].extend(blue_start)
            self['blue_stop'].extend(blue_stop)
            self['cont_start'].extend(blue_start)
            self['cont_stop'].extend(blue_stop)

            self['ncont']+=1
        # add new red cont
        if (type(red_start)!=type(None)) & (type(red_stop)!=type(None)):
            #assert len(red_start)==len(red_stop), 'Non-equal number of start/stop positions for new red continuum'

            self['red_start'].extend(red_start)
            self['red_stop'].extend(red_stop)
            self['cont_start'].append(red_start)
            self['cont_stop'].append(red_stop)

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

    def plotIndexBands(self, ax=None, scale=1.0, alpha=0.1, contCol='blue', indCol='red', \
                       ymin=-100.0, ymax=100.0, autoy=False, noContinuua=False, \
                       addLabel=True, justLine=False, linePos=0.9, labelPos=0.95, \
                       labelSize=20.0, showHorizLine=False, usePrettyPrint=True):
        """
        Overplot the regions used for continuum and feature bandpasses
        
        """

        if ax is None:
            fig, ax=pl.subplots()

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
                 blue_stop=None, red_start=None, red_stop=None, table=None, dicts=None, verbose=True):
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
        print "Augmenting {} with ind: {}, {}; blue: {}, {}; red {}, {}".format(name, ind_start, ind_stop, blue_start, blue_stop, red_start, red_stop)
        self[name].augmentIndex(ind_start=ind_start, ind_stop=ind_stop, blue_start=blue_start, blue_stop=blue_stop, red_start=red_start, red_stop=red_stop)


        # # sanity checks
        # assert np.any(map(lambda a: a is not None, locals())), "One of the required inputs was not defined"
        # assert hasattr(self, name), "This index is not defined, so we cannot augment it"

        # # add new features
        # if (type(ind_start)!=type(None)) & (type(ind_stop)!=type(None)):
        #     #assert len(list(ind_start))==len(ind_stop), 'Non-equal number of start and stop positions for new feature'
        #     self['ind_start'].extend([ind_start])
        #     self['ind_stop'].extend([ind_stop])
        #     self['nfeat']+=1
        # # add new blue cont
        # if (type(blue_start)!=type(None)) & (type(blue_stop)!=type(None)):
        #     #assert len(blue_start)==len(blue_stop), 'Non-equal number of start/stop positions for new blue continuum'
        #     self['blue_start'].extend([blue_start])
        #     self['blue_stop'].extend([blue_stop])
        #     self['cont_start'].extend([blue_start])
        #     self['cont_stop'].extend([blue_stop])
        #     self['ncont']+=1
        # # add new red cont
        # if (type(red_start)!=type(None)) & (type(red_stop)!=type(None)):
        #     #assert len(red_start)==len(red_stop), 'Non-equal number of start/stop positions for new red continuum'
        #     self['red_start'].extend([red_start])
        #     self['red_stop'].extend([red_stop])
        #     self['cont_start'].extend([red_start])
        #     self['cont_stop'].extend([red_stop])
        #     self['ncont']+=1
        # # label as NOT SIMPLE
        # self['simpleIndex']=False
        # if verbose: print "Augmented index "+name+" to ", ind

    def add_simple_indices_via_table(self, verbose=False):
        """
        Add multiple indices through a table
        """


        nind = len(self.table['name'])
        current_index=None


        for ni in xrange(nind):
            if self.table['name'][ni] != '-':
                self.add_simple_index(name=self.table['name'][ni], ind_start=self.table['ind_start'][ni], ind_stop=self.table['ind_stop'][ni], \
                                      blue_start=self.table['blue_start'][ni], blue_stop=self.table['blue_stop'][ni], \
                                      red_start=self.table['red_start'][ni], red_stop=self.table['red_stop'][ni], verbose=verbose)
                current_index = self.table['name'][ni]
            else:

                assert current_index is not None, "No index preceeding index with name=='-'"
                self.augment_index(name=current_index, ind_start=self.table['ind_start'][ni], ind_stop=self.table['ind_stop'][ni], \
                                   blue_start=self.table['blue_start'][ni], blue_stop=self.table['blue_stop'][ni], \
                                   red_start=self.table['red_start'][ni], red_stop=self.table['red_stop'][ni], verbose=verbose)

    def add_simple_index_via_dict(self, index_dict, wavesyst=None, verbose=False):
        """
        Add a single index via (unspecified) dict class

        This is likely surplus 
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
        
        
        

def getLickIndicesAir(filename="/Data/stellarpops/index_definitions/lickIndicesAir.txt",verbose=False):
    """
    Load the Lick indices from Worthey's website (http://astro.wsu.edu/worthey/html/index.table.html, air wavelengths).
    These are a compilation from Trager et al. (1998) and Worthey & Ottaviani (1997). They are the same as used in TMJ10 models. 

    """

    tab = Table.read(filename, format='ascii')
    inds = indlib(table=tab, verbose=verbose)
    return inds

def getLickIndicesVac(filename="/Data/stellarpops/index_definitions/lickIndicesAir.txt",verbose=False):
    """
    Like above but convert to vacuum wavelengths
    """
    
    table = Table.read(filename, format='ascii')
    nline=len(table)
    ncol=len(table.colnames)
    for nli in xrange(nline):
        for nci in xrange(1,ncol-2):
            table[nli][nci] = air2vac(table[nli][nci],verbose=verbose)

    inds = indlib(table=table, verbose=verbose)


    return inds


def getCvD12IndicesVac(filename="/Data/stellarpops/index_definitions/CvDIndicesVac.txt", verbose=True):
    """
    Load the indicies presented in Table 1 of Conroy & van Dokkum 2012a

    Vacuum wavelengths

    """
    tab = Table.read(filename, format='ascii')

    inds=indlib(table=tab,verbose=verbose)



    # #Make a combind CaT index from CaII86_1, _2 and _3

    # #Start with a simple index of definition 1
    # inds.add_simple_index(name='CaT', ind_start=inds.CaII86_1['ind_start'], ind_stop=inds.CaII86_1['ind_stop'], blue_start=inds.CaII86_1['blue_start'], \
    #               blue_stop=inds.CaII86_1['blue_stop'], red_start=inds.CaII86_1['red_start'], red_stop=inds.CaII86_1['red_stop'])

    # #Augment the index with definitions 2 and 3
    # inds.augmentIndex(ind_start=inds.CaII86_2['ind_start'], ind_stop=inds.CaII86_2['ind_stop'], blue_start=inds.CaII86_2['blue_start'], \
    #               blue_stop=inds.CaII86_2['blue_stop'], red_start=inds.CaII86_2['red_start'], red_stop=inds.CaII86_2['red_stop'])
    # inds.augmentIndex(ind_start=inds.CaII86_3['ind_start'], ind_stop=inds.CaII86_3['ind_stop'], blue_start=inds.CaII86_3['blue_start'], \
    #               blue_stop=inds.CaII86_3['blue_stop'], red_start=inds.CaII86_3['red_start'], red_stop=inds.CaII86_3['red_stop'])




    return inds

def getCvD12IndicesAir(filename="/Data/stellarpops/index_definitions/CvDIndicesVac.txt", verbose=True):
    """
    Load the indicies presented in Table 1 of Conroy & van Dokkum 2012a

    Vacuum wavelengths

    """
    tab = Table.read(filename, format='ascii')
    nline=len(tab)
    ncol=len(tab.colnames)
    for nli in xrange(nline):
        for nci in xrange(1,ncol):
            tab[nli][nci] = vac2air(tab[nli][nci],verbose=verbose)



    inds=indlib(table=tab,verbose=verbose)



    # self.ca = np.array((8484., 8513., 8522., 8562., 8642., 8682.))
    # self.pa = np.array((8461., 8474., 8577., 8619., 8730., 8772.))
    # self.cacont = np.array((8474.,8484.,8563.,8577.,8619.,
    #                     8642.,8700.,8725.,8776.,8792.))

    inds.add_simple_index('CaT', )
    
    #Start with a simple index of definition 1
    inds.add_simple_index(name='CaT', ind_start=inds.CaII86_1['ind_start'], ind_stop=inds.CaII86_1['ind_stop'], blue_start=inds.CaII86_1['blue_start'], \
                  blue_stop=inds.CaII86_1['blue_stop'], red_start=inds.CaII86_1['red_start'], red_stop=inds.CaII86_1['red_stop'])



    #Augment the index with definitions 2 and 3
    inds.CaT.augmentIndex(ind_start=inds.CaII86_2['ind_start'], ind_stop=inds.CaII86_2['ind_stop'])

    inds.CaT.augmentIndex(ind_start=inds.CaII86_3['ind_start'], ind_stop=inds.CaII86_3['ind_stop'], blue_start=inds.CaII86_3['blue_start'], \
                  blue_stop=inds.CaII86_3['blue_stop'], red_start=inds.CaII86_3['red_start'], red_stop=inds.CaII86_3['red_stop'])

    inds.CaT.augmentIndex(red_start=inds.CaII86_last_cont_band['red_start'], red_stop=inds.CaII86_last_cont_band['red_stop'])


    

    return inds





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


