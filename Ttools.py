"""
Ryan Houghton 24/8/2016

Functions to read in the Thomas 2003/2010 absorption line indices

loadT03Models - loads the 2003 models
loadT10Models - loads the 2010 models
loadModels    - primary function for above 
showMgFePrime - investigate the Z invariance of MgFe'
showZgrid     - plot a classic grid of models
compareHbMgFe - compare the Hb-MgFe grid for the primed and unprimed
                MgFe indices, using different Fe absorption features

"""


import pylab as pl
import numpy as np
import atpy as ap
import pdb
from stellpops import indexTools as it
from stellpops import M11tools as m11
from stellpops.indexTools import drawGrid

# setup the base directory
basedir="/home/houghton/z/data/stellar_pops/"


def loadT03Models(dir=basedir+'TMB03/', file='alpha-models', abundance='', \
                  Fe1='Fe5270', Fe2='Fe5335', rawTable=False, verbose=True):
    """
    Note that these indices are calibrated to the Lick/IDS system (not from flux calibrated spectra).
    """
    tab = loadModels(dir=dir, file=file, Fe1=Fe1, Fe2=Fe2, rawTable=rawTable, verbose=verbose)
    return tab

def loadT10Models(dir=basedir+'TMJ10/', resolution='Lick', resDict={'Lick':'tmj', 'MILES':'tmj_HighRes'}, abundance='', \
                  Fe1='Fe5270', Fe2='Fe5335', rawTable=False, verbose=True):
    """
    These indices are calibrated to normal flux-calibrated spectra, so no need for Lick standard stars, etc. 
    """
    tab = loadModels(dir=dir, file=resDict[resolution], Fe1=Fe1, Fe2=Fe2, rawTable=rawTable, verbose=verbose)
    return tab

def loadModels(dir=basedir+'TMB03/', file='alpha-models', abundance='', suffix='.dat', \
                    Fe1='Fe5270', Fe2='Fe5335', Fe3='Fe5406', rawTable=False, verbose=True):
    """
    RH 23/6/16

    Load the TMB03 or TMJ10 basic Lick Indices for varying alpha/Fe and Fe/H with age.

    Abundance can be '', -calcium, -carbon, -nitrogen
    
    """

    tab = ap.Table(dir+file+abundance+suffix, type='ascii')
    tab2 = ap.Table(dir+file+abundance+suffix, type='ascii') # 2nd table for removing columns
    tab2.remove_columns(['age', '[Z/H]','[alpha/Fe]'])
    
    # add [MgFe], [MgFe]'
    tab2.add_column('MgFeP', it.calcMgFePrime(tab2.Mgb, tab2[Fe1], tab2[Fe2]))
    tab2.add_column('MgFe', it.calcMgFe(tab2.Mgb, tab2[Fe1], tab2[Fe2]))
    if Fe3 is not None: tab2.add_column('MgFe3', it.calcMgFe(tab2.Mgb, tab2[Fe1], tab2[Fe2], tab2[Fe3]))
    tab2.add_column('meanFe', it.calcMeanFe(tab2[Fe1], tab2[Fe2]))
    if Fe3 is not None: tab2.add_column('mean3Fe', it.calcMeanFe(tab2[Fe1], tab2[Fe2], tab2[Fe3]))
    
    if verbose: print 'Available rows are: ', tab2.keys()
    
    if rawTable:
        newtab=tab
    else:
        Zs = list(set(tab['[Z/H]']))
        As = list(set(tab['[alpha/Fe]']))
        Zs.sort()
        As.sort()

        if verbose:
            print "Found [Z/H] values of ", Zs
            print 'Found [alpha/Fe] values of ', As

        newtab={}
        for zi in Zs:
            minitab={}
            for ai in As:
                loc = np.where((tab['[alpha/Fe]']==ai) & (tab['[Z/H]']==zi))[0]
                minitab['A='+str(ai)]=tab2.rows(loc)
            newtab['Z='+str(zi)]=minitab
            if verbose: print 'Length of Z='+str(zi)+', A='+str(ai)+' data is '+str(len(loc))

        # add helper keys
        newtab['Zkeys']=np.array(['Z='+str(x) for x in Zs]) # don't just call tab.keys() here as order gets scrambled compared to Zs, As
        newtab['Akeys']=np.array(['A='+str(y) for y in As])
        
        newtab['Zs']=np.array(Zs)
        newtab['As']=np.array(As)

        newtab['age']=tab.age[loc]

    return newtab

def showMgFePrime(ages=[12.0], iron1='Fe5270', iron2='Fe5335', rescale=False, tab=None):
    """
    Show the Z invariance of [MgFe]'

    Input:
      ages  - make plot(s) for this/these ages
      iron1 - use this iron index as one of the indices in <fe>
      iron2 - use this iron index as the other index in <fe>
      rescale - subtract off the median value of the caluclated MgFe index
                to aid visualising differences
      tab   - quickly load the models by passing the model table here
    """

    # load the T10 models
    if tab is None:
        tab = tt.loadT10Models()

    Zkeys = tab['Zkeys']
    Zs = tab['Zs']
    Akeys = tab['Akeys']
    As = tab['As']

    for age in ages:
        for zki,zi in zip(Zkeys,Zs):
            index=[]
            alpha=[]
            for aki, ai in zip(Akeys,As):
                loc = np.where(tab['age']==age)[0]
                index.append( it.calcMgFePrime( tab[zki][aki].Mgb[loc], tab[zki][aki][iron1][loc], tab[zki][aki][iron2][loc] ) )
                alpha.append(ai)
            if rescale: index = index-np.median(index)
            pl.plot(alpha,index, label="Age="+str(age)+', [Z/Fe]='+str(zi))

    pl.legend(loc=0)
    if rescale:
        pl.ylabel(r'$\Delta$[MgFe]$^\prime$')
    else:
        pl.ylabel("[MgFe]'")
    pl.xlabel('[alpha/Fe]')


def showZgrid(index1, index2, minage=1.0, maxage=15.0, minZ=-1, maxZ=1.0, alpha=0.0, color="k", \
              ageLabels=True, Zlabels=True, labelsize=10, ZforAgeLabel='min', ageForZlabel='min', \
              aha='right', ava='bottom', zha='right', zva='top', showAlpha=False, \
              ageZforAlphaLabel=['max','max'], Fe1='Fe5270', Fe2='Fe5335', tab=None, \
              alphalabs={'A=0.0':r'[$\alpha$/Fe]=0.0', 'A=0.3':r'[$\alpha$/Fe]=0.3', 'A=0.5':r'[$\alpha$/Fe]=0.5'}):
    """
    Plot a classic grid for abundances

    Inputs:
      index1 - the first index on the x-axis
      index2 - the second index on the y-axis
      minage - minimum age to grid
      maxage - maximum age to grid
      minZ   - minimum Z to grid
      maxZ   - maximum Z to grid
      alpha  - set the alpha/Fe value to use for grid
      colour - colour for grid
      ageLabels - add age labels to grid
      Zlabels   - add Z labels to grid
      ZforAgeLabel - where to put the Age labels: min or max
      ageForZlabel - where to put the Z labels: min or max
      aha    - alpha label horizontal alignment
      ava    - alpha label vertical alignment
      zha    - Z label horizontal alignment
      zva    - Z label vertical alignment
      showAlpha - show the Alpha label
      ageZforAlphaLabel - where to position the alpha label: 2 elem list of min or max
      tab    - provide the Thomas models table here to speed things up
      
    """

    # quick load of models
    if tab is None:
        tab = loadT10Models(Fe1=Fe1, Fe2=Fe2)

    # make alpha key:
    alphaKey='A='+str(alpha)

    # determine label positions
    Zkeys = tab['Zkeys']
    Zs = tab['Zs']
    Zloc = np.where((Zs>=minZ) & (Zs<maxZ))[0]
    ages = tab['age']
    ageloc = np.where((ages>=minage) & (ages<=maxage))[0]

    # draw the primary lines of the grid
    array1 = np.zeros((len(ageloc),len(Zloc)))
    array2 = np.zeros_like(array1)
    for n, zki in enumerate(Zkeys[Zloc]):
        array1[:,n] = tab[zki][alphaKey][index1][ageloc]
        array2[:,n] = tab[zki][alphaKey][index2][ageloc]

    # put age on last axis, as is the norm
    array1=array1.T
    array2=array2.T
    
    drawGrid(index1, index2, array1, array2, ages[ageloc], Zs[Zloc], alpha, color=color,
             ageLabels=ageLabels, Zlabels=Zlabels, labelsize=labelsize, \
             ZforAgeLabel=ZforAgeLabel, ageForZlabel=ageForZlabel, \
             aha=aha, ava=ava, zha=zha, zva=zva, showAlpha=showAlpha, \
             ageZforAlphaLabel=ageZforAlphaLabel, ageSuffix="Gyr", Zprefix='[Z/H]', alphaPrefix=r'[$\alpha$/Fe]')


def compareHbMgFe(minage=1.0, maxage=15.0):
    """
    RH 24/8/2016

    Show the Hb - MgFe grid for the [MgFe] and [MgFe]' indices...
    it looks like the origial [MgFe] is better  for some Fe indices
    (less change with alpha/Fe)
    
    """

    # load models
    tab1 = loadT10Models()
    tab2 = loadT10Models(Fe2='Fe5406')

    # make figure
    pl.figure(figsize=(10,10))
    # show original Hb MgFe
    pl.subplot(221)
    showZgrid('MgFe', 'Hb', alpha=0.0, color='blue', ageForZlabel='max', tab=tab1, \
              minage=minage, maxage=maxage)
    showZgrid('MgFe', 'Hb', alpha=0.3, color='red', ageForZlabel='max', Zlabels=False, ageLabels=False, tab=tab1, \
              minage=minage, maxage=maxage)
    pl.title('Original MgFe index')
    # show the new Hb MgFe' 
    pl.subplot(222)
    
    showZgrid('MgFeP', 'Hb', alpha=0.0, color='blue', ageForZlabel='max', tab=tab1, \
              minage=minage, maxage=maxage)
    showZgrid('MgFeP', 'Hb', alpha=0.3, color='red', ageForZlabel='max', Zlabels=False, ageLabels=False, tab=tab1, \
              minage=minage, maxage=maxage)
    pl.title('New Primed Index')
    # show original Hb MgFe but using different Fe
    pl.subplot(223)
    showZgrid('MgFe', 'Hb', alpha=0.0, color='blue', ageForZlabel='max', tab=tab2, \
              minage=minage, maxage=maxage)
    showZgrid('MgFe', 'Hb', alpha=0.3, color='red', ageForZlabel='max', Zlabels=False, ageLabels=False, tab=tab2, \
              minage=minage, maxage=maxage)
    pl.title('Original MgFe index with Fe5406')
    # show new Hb MgFe using different Fe
    pl.subplot(224)
    showZgrid('MgFeP', 'Hb', alpha=0.0, color='blue', ageForZlabel='max', tab=tab2, \
              minage=minage, maxage=maxage)
    showZgrid('MgFeP', 'Hb', alpha=0.3, color='red', ageForZlabel='max', Zlabels=False, ageLabels=False, tab=tab2, \
              minage=minage, maxage=maxage)
    pl.title('New Primed Index with Fe5406')


def compareT03T11():
    """
    RH 18/10/2016

    Compare old and new models. May also help me figure out if T11 is at Lick resolution.

    """

    t03 = loadT03Models()
    t11 = loadT10Models()

    pl.figure()
    pl.plot(t03['Z=0.0']['A=0.0']['MgFeP'], t03['Z=0.0']['A=0.0']['Fe5335'], "k-")
    pl.plot(t11['Z=0.0']['A=0.0']['MgFeP'], t03['Z=0.0']['A=0.0']['Fe5335'], "r-")
    pl.plot(t03['Z=0.0']['A=0.3']['MgFeP'], t03['Z=0.0']['A=0.3']['Fe5335'], "k:")
    pl.plot(t11['Z=0.0']['A=0.3']['MgFeP'], t03['Z=0.0']['A=0.3']['Fe5335'], "r:")
    pl.plot(t03['Z=0.0']['A=0.5']['MgFeP'], t03['Z=0.0']['A=0.5']['Fe5335'], "k--")
    pl.plot(t11['Z=0.0']['A=0.5']['MgFeP'], t03['Z=0.0']['A=0.5']['Fe5335'], "r--")

    pdb.set_trace()


def calcDispCorsForLickWithM11(indices=None, maxVelDisp=500.0, minVelDisp=100.0, nSample=21, ageMinMax=[1.0, 15.0], order=4, doPlot=False):
    """
    RH 27/10/2016

    Calculate dispersion correctors using the Maraston 2011 models for all Lick indices

    """

    # load Lick indices: using AIR wavelengths and LICK resolutions
    libLick = it.loadLickIndicesAir(atLickRes=True)
    
    # load M11 models
    ms = m11.loadM11ssps(lib='MILES')

    # for each index, calculate a corrector
    cors={}
    if indices is None:
        indices=libLick.names
    
    for ind in indices:
        cors[ind] = it.calcVelDispIndexCorrection(ms, libLick[ind], maxVelDisp=maxVelDisp, minVelDisp=minVelDisp, normVelDisp=None, \
                                                      nSample=nSample, ageMinMax=ageMinMax, order=order, doPlot=doPlot)
    return cors
