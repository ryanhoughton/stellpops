import numpy as np 
import matplotlib.pyplot as plt 
from Miles.handle_miles_models import load_eMILES_spec
import glob

basedir='/Data/stellar_pops/Miles'
files=glob.glob("{}/*".format(basedir))

spectra=[]

for file in files:

    spec=load_eMILES_spec(file)
    disp=spec.lam[1]-spec.lam[0]

    print "Z={}, [NaFe]={}".format(spec.Z, spec.NaFe)
    spec.irindex(disp, 'FeH', vac_or_air='air', verbose=True)
    spec.irindex(disp, 'NaIsdss', vac_or_air='air', verbose=True)
    print "\n"
    spectra.append(spec)





