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


