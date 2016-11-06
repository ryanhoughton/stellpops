import numpy as np 
import specTools as s
from astropy.io import fits
import os
import scipy.interpolate as si

import glob

basedir='/mnt/sda11/Miles/models'


def _load_single_eMILES_spec(filename, basedir='/mnt/sda11/Miles/models'):

    """Load an e-Miles file and return a spectrum class, with attributes like Z, age, etc accessible

    Inputs:
    -- A filename. Can be full path or just the file. 

    Outputs:
    -- A spectrum object
    """

    model_name=os.path.basename(filename)


    userdict={}

    if filename.split('.')[-1]=='fits':
        hdu=fits.open('{}'.format(filename))
        header=hdu[0].header
        data=hdu[0].data


        lamdas=header['CRVAL1']+(np.arange(header['NAXIS1'], dtype=float)+1.0-header['CRPIX1'])*header['CDELT1']
        #lamdas=np.array([(x+1-header['CRPIX1'])*header['CDELT1']+header['CRVAL1'] for x in xrange(len(data))])

    else:
        lamdas, data=np.genfromtxt('{}'.format(filename), unpack=True, skip_header=61)

    if 'bi' in model_name:
        userdict['IMF_type']='bimodal'
        IMF_pos=model_name.index("bi")+2

        userdict['IMF']=model_name[IMF_pos:IMF_pos+4]

    if 'NaFe' in model_name:
        nafe_pos=model_name.index('NaFe')+4
        userdict['NaFe']=model_name[nafe_pos:nafe_pos+3]
    else:
        userdict['NaFe']=0.0



    z_pos=model_name.index('Z')+1
    userdict['Z']=model_name[z_pos:z_pos+5]

    age_pos=model_name.index('T')+1
    userdict['age']=model_name[age_pos:age_pos+6]

    spec=s.spectrum(lamdas, data, userdict=userdict)
    return spec

def load_eMILES_spectra(basedir='/mnt/sda11/Miles/models'):

    #Solar metallicity

    IMFs=['bi0.30', 'bi0.50', 'bi0.80', 'bi1.00', 'bi1.30', 'bi1.50', 'bi1.80', 'bi2.00']
    Zs=['m2.32', 'm1.71', 'm1.31', 'm0.71', 'm0.40', 'p0.00', 'p0.22']
    specs={}

    #Work out the lamda array and get the length of the spectrum
    hdu=fits.open('{}'.format("{}/E{}Z{}T17.7828_iPp0.00_baseFe.fits".format(basedir, 'bi0.30', 'p0.00')))
    header=hdu[0].header
    data=hdu[0].data    
    lamdas=header['CRVAL1']+(np.arange(header['NAXIS1'], dtype=float)+1.0-header['CRPIX1'])*header['CDELT1']



    for IMF in IMFs:
        for i, Z in enumerate(Zs):

            #Get all the files with Z=Z and IMF=IMF
            files=glob.glob("{}/E{}Z{}T??.????_iPp0.00_baseFe.fits".format(basedir, IMF, Z))

            #Sort so the ages are in order
            files.sort(key=lambda x:float(x.split('T')[-1][:7]))

            #Empty arrays to store all the ages and flams
            ages=np.empty(len(files))
            flams=np.empty([len(Zs), len(files), len(data)])
            for j, file in enumerate(files):

                #Base MILES are fits files, [Na/Fe] enhanced are ascii
                if file.split('.')[-1]=='fits':
                    hdu=fits.open('{}'.format(file))
                    data=hdu[0].data                 

                else:
                    lamdas, data=np.genfromtxt('{}'.format(file), unpack=True, skip_header=61)


                age_pos=file.index('T')+1
                ages[j]=float(file[age_pos:age_pos+6])


                flams[i, j, :]=data

        sp=s.spectrum(lamspec=flams, lam=lamdas, age=ages, Z=Zs, IMF=IMF, model='eMILES', wavesyst="air")
        specs[IMF]=sp


    return specs










def convolve_eMILES_spectra(spec, sigma, verbose=False):

    """Convolve an e-MILES spectrum up to a given sigma. This is made slightly tricky by the fact that the e-MILES spectra are at a fixed FWHM of lambda=2.5A
    below 8950A, but at sigma=60 km/s above that. We split the spectra at 8950A, convolve both appropriately, then join back together.

    Inputs:
    -- spec: an e-MILES spectrum object
    -- sigma: a velocity dispersion to convolve to
    -- verbose: a boolean

    Outputs:
    -- an e-Miles spectrum object

    Info:
    -- Wavelength should be in Angstroms


    ###### TO DO ######

    Sort out what happens at 8950 angstroms- at the moment, the two convolutions do funny things around the join
    """

    #Do this to keep the userdict of the spectrum, instead of making a new spectrum and losing that dictionary
    import copy
    final_spec=copy.deepcopy(spec)

    lamdas=np.array(spec.lam)
    data=np.array(spec.flam)


    assert (lamdas[0]<8950.0) & (lamdas[-1]>8950.0), "Are the wavelengths in Angstroms? Looks like the wavelength array doesn't contain 8950 A"

    lamda_mask=lamdas<8950.0

    

    final_flams=[]

    


    low_lam_spec=s.spectrum(lamdas[lamda_mask], data[lamda_mask])
    high_lam_spec=s.spectrum(lamdas[~lamda_mask], data[~lamda_mask])

    #e-MILES spectral res above 8950A is 60km/s
    high_lam_specsig=60.0
    assert sigma>high_lam_specsig, "Sigma must be greater than the spectral resolution of the models"
    convsig=np.sqrt(sigma**2 - high_lam_specsig**2)

    high_lam_spec.gaussVelConvolve(0.0, convsig, verbose=verbose)
    # loglams=high_lam_spec.loglam      
    # interp=si.interp1d(loglams, high_lam_spec.conflam[0], fill_value='extrapolate')
    # linlams=high_lam_spec.lam
    # high_lam_flam=interp(linlams)


    import ipdb; ipdb.set_trace()

    #high_lam_spec.interp_log_to_lin(verbose=verbose)

    #e-MILES spectral res below 8950A is 2.5A
    
    #For now, we only care about the spectra around NaI 8190
    #At this wavelength, with a FWHM of 2.5A, delta v is 39km/s

    low_lam_specsig=39.0
    assert sigma>low_lam_specsig, "Sigma must be greater than the spectral resolution of the models"
    convsig=np.sqrt(sigma**2 - high_lam_specsig**2)

    low_lam_spec.gaussVelConvolve(0.0, convsig, verbose=verbose)
    # loglams=low_lam_spec.loglam  
    # interp=si.interp1d(loglams, low_lam_spec.conflam[0], fill_value='extrapolate')
    # linlams=low_lam_spec.lam
    # low_lam_flam=interp(linlams)
    #low_lam_spec.interp_log_to_lin(verbose=verbose)


    final_flams.append(np.concatenate([low_lam_flam, high_lam_flam]))


    final_spec.flam=np.array(final_flams)
    setattr(final_spec, 'sigma', sigma)



    return final_spec





    
if __name__=='__main__':
    #Used for testing things

    import glob

    specs=load_eMILES_spectra()

    spec=specs['bi1.30']
    import ipdb; ipdb.set_trace()

    spec.gaussVelConvolve(0.0, 200.0)



    """

    for i, file in enumerate(files):

        print i
        
        spec=load_eMILES_spec(file)
        spec2=convolve_eMILES_spectra(spec, 200.0)

        import ipdb; ipdb.set_trace()
    """









