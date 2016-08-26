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

# setup the base directory
basedir="/home/houghton/z/data/stellar_pops/"


def loadT03Models(dir=basedir+'TMB03/', file='alpha-models', abundance='', \
                  Fe1='Fe5270', Fe2='Fe5335', rawTable=False, verbose=True):
    """
    Note that these indices are calibrated to the Lick/IDS system (not from flux calibrated spectra).
    """
    tab = loadModels(dir=dir, file=file, Fe1=Fe1, Fe2=Fe2, rawTable=rawTable, verbose=verbose)
    return tab

def loadT10Models(dir=basedir+'TMJ10/', file='tmj', abundance='', \
                  Fe1='Fe5270', Fe2='Fe5335', rawTable=False, verbose=True):
    """
    These indices are calibrated to normal flux-calibrated spectra, so no need for Lick standard stars, etc. 
    """
    tab = loadModels(dir=dir, file=file, Fe1=Fe1, Fe2=Fe2, rawTable=rawTable, verbose=verbose)
    return tab

def loadModels(dir=basedir+'TMB03/', file='alpha-models', abundance='', suffix='.dat', \
                    Fe1='Fe5270', Fe2='Fe5335', rawTable=False, verbose=True):
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
    tab2.add_column('meanFe', it.calcMeanFe(tab2[Fe1], tab2[Fe2]))
    
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
            print 'Found [alpha/H] values of ', As

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


def showZgrid(index1, index2, minage=1.0, maxage=15.0, minZ=-1, maxZ=1.0, alpha='A=0.0', color="k", \
              ageLabels=True, Zlabels=True, labelsize=10, ZforAgeLabel='min', ageForZlabel='min', \
              aha='right', ava='bottom', zha='right', zva='top', showAlpha=False, \
              ageZforAlphaLabel=['max','max'], Fe1='Fe5270', Fe2='Fe5335', tab=None):
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
        tab = tt.loadT10Models(Fe1=Fe1, Fe2=Fe2)

    # determine label positions
    Zkeys = tab['Zkeys']
    Zs = tab['Zs']
    Zloc = np.where((Zs>=minZ) & (Zs<maxZ))[0]
    ages = tab['age']
    ageloc = np.where((ages>=minage) & (ages<=maxage))[0]

    if ZforAgeLabel=='min':
        z4agelab = np.min(Zs[Zloc])
    else:
        z4agelab = np.max(Zs[Zloc])
    if ageForZlabel=='min':
        age4zlab=np.min(ages[ageloc])
    else:
        age4zlab=np.max(ages[ageloc])
    age4zlabloc = np.where(ages==age4zlab)[0]

    # draw the primary lines of the grid
    array1 = np.zeros((len(ageloc),len(Zloc)))
    array2 = np.zeros_like(array1)
    for n, zki in enumerate(Zkeys[Zloc]):
        array1[:,n] = tab[zki][alpha][index1][ageloc]
        array2[:,n] = tab[zki][alpha][index2][ageloc]
        pl.plot(tab[zki][alpha][index1][ageloc], tab[zki][alpha][index2][ageloc], "-", color=color)

        if ageLabels & (z4agelab==Zs[Zloc[n]]):
            for l in ageloc:
                pl.text(tab[zki][alpha][index1][l], tab[zki][alpha][index2][l], str(ages[l])+" Gyr", \
                        fontsize=labelsize, horizontalalignment=aha, verticalalignment=ava)

        if Zlabels:
            pl.text(tab[zki][alpha][index1][age4zlabloc], tab[zki][alpha][index2][age4zlabloc], zki, \
                    fontsize=labelsize, horizontalalignment=zha, verticalalignment=zva)
            
    # now draw connectors at fixed Z, rather than fixed age
    at1=array1#.T
    at2=array2#.T
    for v1, v2 in zip(at1,at2):
        pl.plot(v1, v2, ":", color=color)

    # show alpha value if asked
    if showAlpha:
        if ageZforAlphaLabel[0]=='min':
            age4AlphaLab = np.argmin(ages[ageloc])
        else:
            age4AlphaLab = np.argmax(ages[ageloc])
        if ageZforAlphaLabel[1]=='min':
            z4AlphaLab = np.argmin(Zs[Zloc])
        else:
            z4AlphaLab = np.argmax(ages[Zloc])
        
        pl.text(tab[Zkeys[Zloc[z4AlphaLab]]][alpha][index1][ageloc[age4AlphaLab]],\
                tab[Zkeys[Zloc[z4AlphaLab]]][alpha][index2][ageloc[age4AlphaLab]], \
                alpha, color=color)
          
    pl.xlabel(index1)
    pl.ylabel(index2)
    

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
    showZgrid('MgFe', 'Hb', alpha='A=0.0', color='blue', ageForZlabel='max', tab=tab1, \
              minage=minage, maxage=maxage)
    showZgrid('MgFe', 'Hb', alpha='A=0.3', color='red', ageForZlabel='max', Zlabels=False, ageLabels=False, tab=tab1, \
              minage=minage, maxage=maxage)
    pl.title('Original MgFe index')
    # show the new Hb MgFe' 
    pl.subplot(222)
    
    showZgrid('MgFeP', 'Hb', alpha='A=0.0', color='blue', ageForZlabel='max', tab=tab1, \
              minage=minage, maxage=maxage)
    showZgrid('MgFeP', 'Hb', alpha='A=0.3', color='red', ageForZlabel='max', Zlabels=False, ageLabels=False, tab=tab1, \
              minage=minage, maxage=maxage)
    pl.title('New Primed Index')
    # show original Hb MgFe but using different Fe
    pl.subplot(223)
    showZgrid('MgFe', 'Hb', alpha='A=0.0', color='blue', ageForZlabel='max', tab=tab2, \
              minage=minage, maxage=maxage)
    showZgrid('MgFe', 'Hb', alpha='A=0.3', color='red', ageForZlabel='max', Zlabels=False, ageLabels=False, tab=tab2, \
              minage=minage, maxage=maxage)
    pl.title('Original MgFe index with Fe5406')
    # show new Hb MgFe using different Fe
    pl.subplot(224)
    showZgrid('MgFeP', 'Hb', alpha='A=0.0', color='blue', ageForZlabel='max', tab=tab2, \
              minage=minage, maxage=maxage)
    showZgrid('MgFeP', 'Hb', alpha='A=0.3', color='red', ageForZlabel='max', Zlabels=False, ageLabels=False, tab=tab2, \
              minage=minage, maxage=maxage)
    pl.title('New Primed Index with Fe5406')
