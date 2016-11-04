import numpy as np
import pylab as pl
import string as s
import re
import specTools as t
import os
import pdb
from os.path import expanduser
# load BC03
import BC03tools as BC03
# copy over M05
import M05tools as M05

# set globals
L_sun = M05.L_sun
basedir='/home/houghton/z/data/stellar_pops/MS11/'
dirprefix='SSP_M11_'


def loadM11specWithBC03tools(fname, dir=basedir+'M11_BCformat/'):
    """
    RH 19/8/2016

    The M11 SSPs (MILES and STELIB, Chabrier IMS) also come in BC03 format.
    So use the BC03tools to load them

    """

    spec = BC03.loadBC03spec(dir+fname, IMF=None, loadMassFile=False)
    return spec

def loadM11specWithM05tools(fname, lib='MILES', dir=basedir+dirprefix, \
                            massfile=None, agecol=0, zcol=1, lamcol=2, fluxcol=3, \
                            angstscale=1.0, fluxscale=1.0, skip=0, \
                            Zdict={"10m4":0.00025, "0001":0.001, "001":0.01, \
                                   "002":0.02, "004":0.04, "007":0.07}):
    """
    RH 19/8/2016

    Load a single M11 SSPs in the standard Maraston format, using the M05 tools

    """

    spec = M05.loadM05spec(dir+lib.upper()+"/"+fname, massfile=massfile, \
                           agecol=agecol, zcol=zcol, lamcol=lamcol, \
                           fluxcol=fluxcol, angstscale=angstscale, fluxscale=fluxscale, skip=skip, \
                           Zdict=Zdict)

    return spec


def loadM11ssps(sedpath=basedir, dir=dirprefix, lib="MILES", \
                imf="salpeter", glob="00{1,2,4}", \
                IMFdict={"salpeter":"ss", "kroupa":"kr", "chabrier":"cha"}, \
                Zdict={"10m4":0.00025, "0001":0.001, "001":0.01, \
                "002":0.02, "004":0.04, "007":0.07}, minAge=0.1, maxAge=None, \
                RESdict={'Pickles':[500.,None], 'STELIB':[None, (3.1,3.4)], 'MILES':[None, 2.54], \
                         'ELODIE':[None, 0.55], 'MARCS':[20000,None]}, verbose=False):
    """

    RH 19/8/2016

    Load multiple M11 SSPs in the standard Maraston format, using the M05 tools.

    Note that:
       - the M11 SSPs come with Chabrier IMF as well as Salpeter and Kroupa.
       - the M11 SSPS are available for different stellar libraries
         (MILES, STELIB, ELODIE, MARCS, Pickles), specified by the 'lib' option.

    Note: due to the change in the way the HB morph is used in M11, not all SSPs have a .bhb or .rhb
    suffix. So I've removed the 'morph' option and instead, you should explicity call it in the glob, e.g.
    e.g. glob='00{01.rhb,1,2,4}' for RHB Z=0.001, 0.01, 0.02, 0.04
    or   glob='{10m4.rhb, 0001.rhb, 001, 002, 004}' for ALL Zs, rhb morph.

    As the HB morph only plays a role at low Z

    NOTE: at least for the MILES modes, the ages for each Z are different. Specifying minAge=0.1 solves this as this is the
    minimum age for the most metal rich models and thereafter, sampling of ages is uniform.

    """

    # get the file name
    filespath=sedpath+dirprefix+lib+"/"+"ssp_M11_"+lib+"."+IMFdict[imf]+"z"+glob
    files, nfiles = t.findfiles(filespath)

    # load the first to get the array size
    if verbose:
        print "Reading data from ..."
        for fname in files: print fname

    
    massfile=expanduser(sedpath)+"/"+"stellarmass."+imf
    if not os.path.isfile(massfile): massfile=None

    specs = []
    #specs.append(M05.loadM05spec(files[0], massfile=massfile, resolution=RESdict[lib]))
    #if verbose: print "Read "+files[0]+" (Z="+str(specs[-1].Z)+")"

    # load in all the others
    for fname in files:
        specs.append(M05.loadM05spec(fname, massfile=massfile, resolution=RESdict[lib], minAge=minAge, maxAge=maxAge))
        if verbose: print "Read "+fname+" (Z="+str(specs[-1].Z)+")"

    return specs
