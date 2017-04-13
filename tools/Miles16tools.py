import numpy as np 
import specTools as s
from astropy.io import fits
import os
import scipy.interpolate as si
import scipy.constants as const

import pandas as pd

import glob

#basedir='/Data/stellarpops/Miles/base_set'


def get_all_specs(imf_type, verbose=True):

    specs=load_eMILES_spectra(NaFe=0.0, imf_type=imf_type, verbose=verbose)

    NaFep03_specs=load_eMILES_spectra(NaFe=0.3, verbose=verbose)
    NaFep06_specs=load_eMILES_spectra(NaFe=0.6, verbose=verbose)
    NaFep09_specs=load_eMILES_spectra(NaFe=0.9, verbose=verbose)

    return specs, NaFep03_specs, NaFep06_specs, NaFep09_specs

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

def pandas_read_file(filename):
    dataframe=pd.read_csv(filename, skiprows=61, delim_whitespace=True, names=['lambda', 'data'])
    lamdas, data=dataframe['lambda'].as_matrix(), dataframe['data'].as_matrix()
    return lamdas, data


def load_eMILES_spectra(basedir='/Data/stellarpops/Miles', NaFe=0.0, Zs=None, verbose=True, imf_type='bi'):


    """
    Load all the eMILES spectra and save into a dictionary. The dictionary keys are IMF values (e.g. 'bi1.30'), and the values are spectrum classes
    of spectra with that IMF and a variety of ages. So to get the flams of all the spectra with a Chabrier IMF, you'd do

    spectra['bi1.30'].flam

    """

    import os
    basedir=os.path.expanduser(basedir)
    if NaFe==0.0:
        if imf_type=='bi':
            folder='bi_base_set'
        elif imf_type=='uni':
            folder='uni_base_set'
    elif NaFe==0.3:
        folder='NaFep03/{}'.format(imf_type)
    elif NaFe==0.6:
        folder='NaFep06/{}'.format(imf_type)
    elif NaFe==0.9:
        folder='NaFep09/{}'.format(imf_type)
    elif NaFe==1.2:
        folder='NaFep12/{}'.format(imf_type)
    else:
        raise NameError('NaFe abundance not understood')

    #Solar metallicity

    if Zs==None:
        #Full e-Miles models have 7 metallicities
        #Na enhanced ones only have 3
        if NaFe==0.0:
            Zs=['m2.32', 'm1.71', 'm1.31', 'm0.71', 'm0.40', 'p0.00', 'p0.22']        
        else:
            Zs=['m0.40', 'p0.00', 'p0.22']


    if imf_type=='bi':
        IMFs=['bi0.30', 'bi0.80', 'bi1.00', 'bi1.30', 'bi1.50', 'bi1.80', 'bi2.00', 'bi2.30', 'bi2.50', 'bi2.80', 'bi3.00', 'bi3.30']
    else:
        IMFs=['un0.30', 'un0.80', 'un1.00', 'un1.30', 'un1.50', 'un1.80', 'un2.00', 'un2.30', 'un2.50', 'un2.80', 'un3.00', 'un3.30']

    #Empty dictionary to fill
    specs={}

    #Work out the lamda array and get the length of the spectrum
    hdu=fits.open('{}'.format("{}/uni_base_set/E{}Z{}T17.7828_iPp0.00_baseFe.fits".format(basedir, 'un0.30', 'p0.00')))
    header=hdu[0].header
    data=hdu[0].data    
    lamdas=header['CRVAL1']+(np.arange(header['NAXIS1'], dtype=float)+1.0-header['CRPIX1'])*header['CDELT1']


    if verbose:
        print "Loading eMILES spectra from {}/{}".format(basedir, folder)
    for IMF in IMFs:
        #Empty array to store the flams
        #Get the number of ages by globbing for the first metallicity
        #import pdb; pdb.set_trace()
        tmp=glob.glob("{}/{}/E{}Z{}T??.????_iPp0.00_baseFe*".format(basedir, folder, IMF, Zs[-1]))

        #assert len(tmp)==12, 'Something may be wrong- expecting 12 spectra in with ages > 5 gyr, got {}'.format(len(tmp))
        flams=np.empty([len(Zs), len(tmp), len(data)])
        for i, Z in enumerate(Zs):

            #Get all the files with Z=Z and IMF=IMF
            files=glob.glob("{}/{}/E{}Z{}T??.????_iPp0.00_baseFe*".format(basedir, folder, IMF, Z))

            #Sort so the ages are in order
            files.sort(key=lambda x:float(x.split('T')[-1][:7]))

            #Empty arrays to store all the ages
            ages=np.empty(len(files))
            
            for j, file in enumerate(files):

                #Base MILES are fits files, [Na/Fe] enhanced are ascii
                if file.split('.')[-1]=='fits':
                    hdu=fits.open('{}'.format(file))
                    data=hdu[0].data
                    hdu.close()                

                else:
                    lamdas, data=pandas_read_file('{}'.format(file))

                #import pdb; pdb.set_trace()
                age_pos=file.index('T')+1
                ages[j]=float(file[age_pos:age_pos+6])


                flams[i, j, :]=data


        ages=np.repeat(ages, len(Zs)).reshape(len(files), len(Zs)).T
        Z_values=np.repeat(Zs, len(files)).reshape(len(Zs), len(files))

        userdict={}
        userdict['NaFe']=NaFe

        sp=s.spectrum(lamspec=flams, lam=lamdas, age=ages, Z=Z_values, IMF=IMF, model='eMILES', wavesyst="air", userdict=userdict)
        specs[IMF]=sp
        if verbose:
            print "Loaded {} SSPs with {} IMF".format(len(files)*len(Zs), IMF)


    return specs



def cut_and_measure_index(spec, index, out_sigma, index_type='Cenarro', model_sigma=None, n_sig=10.0):

    """
    Use specTools to cut a long spectrum down to size and measure an index
    """

    #If the index red stop is less than 8950, the model resolution is FWHM of 2.5A

    if model_sigma is None: 
        if np.atleast_1d(np.array(index['red_stop']))[-1]<8950.0:

            if index['nfeat']>0.0:
                model_sigma=const.c*2.5/(np.sqrt(8.*np.log(2.0))*index['ind_start'][0]*1000.0)
            else:
                model_sigma=const.c*2.5/(np.sqrt(8.*np.log(2.0))*index['blue_stop']*1000.0)

            assert out_sigma > model_sigma, 'Cant convolve to a resolution below the model resolution'
            conv_sigma=np.sqrt(out_sigma**2 - model_sigma**2)
            






        #If it's above 8950, the mdoel resolution is 60km/s
        else:
            model_sigma=60.0
            assert out_sigma > model_sigma, 'Cant convolve to a resolution below the model resolution'
            conv_sigma=np.sqrt(out_sigma**2 - model_sigma**2)
    else:
        assert out_sigma > model_sigma, 'Cant convolve to a resolution below the model resolution'
        conv_sigma=np.sqrt(out_sigma**2 - model_sigma**2)


    cutspec=s.cutAndGaussVelConvolve(spec, index, conv_sigma, verbose=False, n_sig=n_sig)
    if index_type=='Cenarro':
        indvals=s.calcCenarroIndex(cutspec, index)
    elif index_type=='Simple':
        indvals=s.calcSimpleIndex(cutspec, index)
    else:
        raise TypeError('Index Type not understood')


    return indvals


def get_np_indices_for_params(IMFs=['bi0.30', 'bi0.50', 'bi0.80', 'bi1.00', 'bi1.30', 'bi1.50', 'bi1.80', 'bi2.00', 'bi2.30', 'bi2.50', 'bi2.80', 'bi3.00', 'bi3.30', 'bi3.50'], \
    Zs=['m2.32', 'm1.71', 'm1.31', 'm0.71', 'm0.40', 'p0.00', 'p0.22'], \
    ages=[5.011,   5.623,   6.309,   7.079,   7.943,   8.912,  10.0,  11.22,   12.589, 14.125,  15.848,  17.782], verbose=False):

    #We want a dictionary where we can get the index for any IMF, Z or age
    #e.g dictionary['bi0.30']=0, dict['m2.32']=0, dict['p0.22']=6, etc


    param_dict={}


    imf_dict=dict(enumerate(IMFs))
    imf_dict=dict((v,k) for k,v in imf_dict.iteritems())    

    Z_dict=dict(enumerate(Zs))
    Z_dict=dict((v,k) for k,v in Z_dict.iteritems())

    age_dict=dict(enumerate(ages))
    age_dict=dict((v,k) for k,v in age_dict.iteritems())

    
    for d in [imf_dict, Z_dict, age_dict]:
        for k, v in d.iteritems():
            if verbose:  
                print k, v
            param_dict[k]=v

    return param_dict


def alpha_enhanced_spec(specs, index, out_sigma, cvd_dir='/Data/stellarpops/CvD1.2', alpha=0.3, verbose=True):

    import CvD12tools as cvd

    cvd13 = cvd.loadCvD12spec('{}/t13.5_solar.ssp'.format(cvd_dir))
    cvdafe = cvd.loadCvD12spec('{}/t13.5_afe+0.2.ssp'.format(cvd_dir))

    CvDlam=cvd13.lam
    alphafac = (cvdafe.flam[3]/cvd13.flam[3]-1)*((10**(alpha)-1.0)/(10**(0.2)-1.0))

    alphafac=s.vac2air(alphafac)

    if verbose:
        print "Created alphafac: shape is {}".format(alphafac.shape)

    #CvD spectra have a spectral resolution of ~2000
    #This is a sigma of ~63km/s

    #Miles models have a FWHM of 2.5A below 8950A, 60km/s above it.
    #60km/s is a resolving power of 2121: R=c/sigma*sqrt(8ln2)=2121
    #At NaI, FWHM of 2.5A is a sigma of 38.9km/s, or R=3276

    
    #Need to discuss- but will not convolve either model above 8950A for the moment
    #Below, I'll convolve MILES up to a sigma of 63 km/s
    #should test whether this makes a difference!

    if index['red_stop']<8950.0:

        if verbose:
            print "Index is below 8950A. Need to convolve MILES up to CvD resolution"

        #Convolve the Miles models up to the CvD resolution
        CvD_sigma=63.0
        if index['nfeat']>0.0:
            model_sigma=const.c*2.5/(np.sqrt(8.*np.log(2.0))*index['ind_start'][0]*1000.0)            
        else:
            model_sigma=const.c*2.5/(np.sqrt(8.*np.log(2.0))*index['cont_stop'][0]*1000.0)

        assert CvD_sigma > model_sigma, 'Cant convolve to a resolution below the model resolution'
        conv_sigma=np.sqrt(CvD_sigma**2-model_sigma**2)

        MILES_at_CvD_res=s.cutAndGaussVelConvolve(specs, index, conv_sigma, n_sig=30.0)

        if verbose:

            print "Made MILES spec at CvD resolution. Shape is {}".format(MILES_at_CvD_res.flam.shape)


    else:
        MILES_at_CvD_res=specs.copy()


    #Cut Miles specs to finish at the same lamda as the CvD specs
    MILES_at_CvD_res.clipSpectralRange(MILES_at_CvD_res.lam[0], CvDlam[-1])

    if verbose:
        print "Clipped the MILES spectra to end at 2.4um. Shape is {}".format(MILES_at_CvD_res.flam.shape)

    #Clip CvDspecs to start and stop at the same lamda as the (clipped or not clipped) MILES ones
    lamda_mask=np.where((CvDlam>MILES_at_CvD_res.lam[0]) & (CvDlam<MILES_at_CvD_res.lam[-1]))[0]
    alphafac=alphafac[lamda_mask]


    alphafac_interp=si.interp1d(CvDlam[lamda_mask], alphafac, fill_value='extrapolate')
    correction=alphafac_interp(MILES_at_CvD_res.lam)   

    newspec=s.spectrum(lam=MILES_at_CvD_res.lam, lamspec=np.exp(np.log(MILES_at_CvD_res.flam)+correction), wavesyst='air')
    # if index['name']=='TiO89' and specs.IMF=='bi1.30':
    #     import matplotlib.pyplot as plt
    #     plt.clf()
    #     plt.figure()
    #     import pdb; pdb.set_trace()
    """
    for i, flam in enumerate(MILES_at_CvD_res.flam.reshape(-1, MILES_at_CvD_res.flam.shape[-1])):

        assert len(flam)==len(correction), "Lengths of arrays aren't equal!"
        # import pdb; pdb.set_trace()
        # import matplotlib.pyplot as plt
        # plt.figure()
        # plt.plot(MILES_at_CvD_res.lam, flam, c='k')
        flam=np.exp(np.log(flam)+correction)
        # plt.plot(MILES_at_CvD_res.lam, flam, c='r')
        # plt.show()
        # import pdb; pdb.set_trace()


        if verbose:
            print "Applied Correction to spec {} of {}".format(i, MILES_at_CvD_res.flam.reshape(-1, MILES_at_CvD_res.flam.shape[-1]).shape[0])
    """

    return newspec

def alpha_corrected_index(specs, index, out_sigma, cvd_dir='/Data/stellarpops/CvD1.2', alpha=0.3, verbose=True):

    newspec=alpha_enhanced_spec(specs, index, out_sigma, cvd_dir='/Data/stellarpops/CvD1.2', alpha=alpha, verbose=True)

    ind_vals=cut_and_measure_index(newspec, index, out_sigma, index_type='Cenarro')

    return ind_vals

    

def fe_enhanced_spec(specs, index, out_sigma, cvd_dir='/Data/stellarpops/CvD1.2', Fe=0.3, verbose=True):

    import CvD12tools as cvd

    cvd13 = cvd.loadCvD12spec('{}/t13.5_solar.ssp'.format(cvd_dir))
    cvdvar = cvd.loadCvD12spec('{}/t13.5_varelem.ssp'.format(cvd_dir))

    CvDlam=cvd13.lam
    fefac = (cvdvar.flam[4]/cvd13.flam[3]-1)*((10**(Fe)-1.0)/(10**(0.3)-1.0))

    fefac=s.vac2air(fefac)

    if verbose:
        print "Created fefac: shape is {}".format(fefac.shape)

    #CvD spectra have a spectral resolution of ~2000
    #This is a sigma of ~63km/s

    #Miles models have a FWHM of 2.5A below 8950A, 60km/s above it.
    #60km/s is a resolving power of 2121: R=c/sigma*sqrt(8ln2)=2121
    #At NaI, FWHM of 2.5A is a sigma of 38.9km/s, or R=3276

    
    #Need to discuss- but will not convolve either model above 8950A for the moment
    #Below, I'll convolve MILES up to a sigma of 63 km/s
    #should test whether this makes a difference!

    if index['red_stop']<8950.0:

        if verbose:
            print "Index is below 8950A. Need to convolve MILES up to CvD resolution"

        #Convolve the Miles models up to the CvD resolution
        CvD_sigma=63.0
        if index['nfeat']>0.0:
            model_sigma=const.c*2.5/(np.sqrt(8.*np.log(2.0))*index['ind_start'][0]*1000.0)            
        else:
            model_sigma=const.c*2.5/(np.sqrt(8.*np.log(2.0))*index['cont_stop'][0]*1000.0)

        assert CvD_sigma > model_sigma, 'Cant convolve to a resolution below the model resolution'
        conv_sigma=np.sqrt(CvD_sigma**2-model_sigma**2)

        MILES_at_CvD_res=s.cutAndGaussVelConvolve(specs, index, conv_sigma, n_sig=30.0)

        if verbose:

            print "Made MILES spec at CvD resolution. Shape is {}".format(MILES_at_CvD_res.flam.shape)


    else:
        MILES_at_CvD_res=specs.copy()


    #Cut Miles specs to finish at the same lamda as the CvD specs
    MILES_at_CvD_res.clipSpectralRange(MILES_at_CvD_res.lam[0], CvDlam[-1])

    if verbose:
        print "Clipped the MILES spectra to end at 2.4um. Shape is {}".format(MILES_at_CvD_res.flam.shape)

    #Clip CvDspecs to start and stop at the same lamda as the (clipped or not clipped) MILES ones
    lamda_mask=np.where((CvDlam>MILES_at_CvD_res.lam[0]) & (CvDlam<MILES_at_CvD_res.lam[-1]))[0]
    fefac=fefac[lamda_mask]

    #interpolate MILES and CvD onto the same wavelength grid
    fefac_interp=si.interp1d(CvDlam[lamda_mask], fefac, fill_value='extrapolate')
    correction=fefac_interp(MILES_at_CvD_res.lam) 

    if verbose:
        print "Interpolated Miles and CvD to lie on the same wavelength array"  

    newspec=s.spectrum(lam=MILES_at_CvD_res.lam, lamspec=np.exp(np.log(MILES_at_CvD_res.flam)+correction), wavesyst='air')
    if verbose:
        print "Made new spectrum"
    # if index['name']=='TiO89' and specs.IMF=='bi1.30':
    #     import matplotlib.pyplot as plt
    #     plt.clf()
    #     plt.figure()
    #     import pdb; pdb.set_trace()
    """
    for i, flam in enumerate(MILES_at_CvD_res.flam.reshape(-1, MILES_at_CvD_res.flam.shape[-1])):

        assert len(flam)==len(correction), "Lengths of arrays aren't equal!"
        # import pdb; pdb.set_trace()
        # import matplotlib.pyplot as plt
        # plt.figure()
        # plt.plot(MILES_at_CvD_res.lam, flam, c='k')
        flam=np.exp(np.log(flam)+correction)
        # plt.plot(MILES_at_CvD_res.lam, flam, c='r')
        # plt.show()
        # import pdb; pdb.set_trace()


        if verbose:
            print "Applied Correction to spec {} of {}".format(i, MILES_at_CvD_res.flam.reshape(-1, MILES_at_CvD_res.flam.shape[-1]).shape[0])
    """

    return newspec

def element_enhanced_spec(specs, index, out_sigma, element, enhancement=0.3, cvd_dir='/Data/stellarpops/CvD1.2', verbose=True):

    import CvD12tools as cvd

    if element=='Na':
        element_index=1


    cvd13 = cvd.loadCvD12spec('{}/t13.5_solar.ssp'.format(cvd_dir))
    cvdvar = cvd.loadCvD12spec('{}/t13.5_varelem.ssp'.format(cvd_dir))

    CvDlam=cvd13.lam
    fefac = (cvdvar.flam[element_index]/cvd13.flam[3]-1)*((10**(enhancement)-1.0)/(10**(0.3)-1.0))

    fefac=s.vac2air(fefac)

    if verbose:
        print "Created fefac: shape is {}".format(fefac.shape)

    #CvD spectra have a spectral resolution of ~2000
    #This is a sigma of ~63km/s

    #Miles models have a FWHM of 2.5A below 8950A, 60km/s above it.
    #60km/s is a resolving power of 2121: R=c/sigma*sqrt(8ln2)=2121
    #At NaI, FWHM of 2.5A is a sigma of 38.9km/s, or R=3276

    
    #Need to discuss- but will not convolve either model above 8950A for the moment
    #Below, I'll convolve MILES up to a sigma of 63 km/s
    #should test whether this makes a difference!

    if index['red_stop']<8950.0:

        if verbose:
            print "Index is below 8950A. Need to convolve MILES up to CvD resolution"

        #Convolve the Miles models up to the CvD resolution
        CvD_sigma=63.0
        if index['nfeat']>0.0:
            model_sigma=const.c*2.5/(np.sqrt(8.*np.log(2.0))*index['ind_start'][0]*1000.0)            
        else:
            model_sigma=const.c*2.5/(np.sqrt(8.*np.log(2.0))*index['cont_stop'][0]*1000.0)

        assert CvD_sigma > model_sigma, 'Cant convolve to a resolution below the model resolution'
        conv_sigma=np.sqrt(CvD_sigma**2-model_sigma**2)

        MILES_at_CvD_res=s.cutAndGaussVelConvolve(specs, index, conv_sigma, n_sig=30.0)

        if verbose:

            print "Made MILES spec at CvD resolution. Shape is {}".format(MILES_at_CvD_res.flam.shape)


    else:
        MILES_at_CvD_res=specs.copy()


    #Cut Miles specs to finish at the same lamda as the CvD specs
    MILES_at_CvD_res.clipSpectralRange(MILES_at_CvD_res.lam[0], CvDlam[-1])

    if verbose:
        print "Clipped the MILES spectra to end at 2.4um. Shape is {}".format(MILES_at_CvD_res.flam.shape)

    #Clip CvDspecs to start and stop at the same lamda as the (clipped or not clipped) MILES ones
    lamda_mask=np.where((CvDlam>MILES_at_CvD_res.lam[0]) & (CvDlam<MILES_at_CvD_res.lam[-1]))[0]
    fefac=fefac[lamda_mask]

    #interpolate MILES and CvD onto the same wavelength grid
    fefac_interp=si.interp1d(CvDlam[lamda_mask], fefac, fill_value='extrapolate')
    correction=fefac_interp(MILES_at_CvD_res.lam) 

    if verbose:
        print "Interpolated Miles and CvD to lie on the same wavelength array"  

    newspec=s.spectrum(lam=MILES_at_CvD_res.lam, lamspec=np.exp(np.log(MILES_at_CvD_res.flam)+correction), wavesyst='air')
    if verbose:
        print "Made new spectrum"
    # if index['name']=='TiO89' and specs.IMF=='bi1.30':
    #     import matplotlib.pyplot as plt
    #     plt.clf()
    #     plt.figure()
    #     import pdb; pdb.set_trace()
    """
    for i, flam in enumerate(MILES_at_CvD_res.flam.reshape(-1, MILES_at_CvD_res.flam.shape[-1])):

        assert len(flam)==len(correction), "Lengths of arrays aren't equal!"
        # import pdb; pdb.set_trace()
        # import matplotlib.pyplot as plt
        # plt.figure()
        # plt.plot(MILES_at_CvD_res.lam, flam, c='k')
        flam=np.exp(np.log(flam)+correction)
        # plt.plot(MILES_at_CvD_res.lam, flam, c='r')
        # plt.show()
        # import pdb; pdb.set_trace()


        if verbose:
            print "Applied Correction to spec {} of {}".format(i, MILES_at_CvD_res.flam.reshape(-1, MILES_at_CvD_res.flam.shape[-1]).shape[0])
    """

    return newspec

def fe_corrected_index(specs, index, out_sigma, cvd_dir='/Data/stellarpops/CvD1.2', Fe=0.3, verbose=True):

    newspec=fe_enhanced_spec(specs, index, out_sigma, cvd_dir='/Data/stellarpops/CvD1.2', Fe=Fe, verbose=True)

    ind_vals=cut_and_measure_index(newspec, index, out_sigma, index_type='Cenarro')

    return ind_vals


def load_mass_table(dir_name='/Data/stellarpops/Miles/Mass_files'):

    """Load the mass file as a table"""


    table_name='ssp_mass_Padova00_BI_baseFe_v10.txt'


    from astropy.io import ascii

    table=ascii.read('{}/{}'.format(dir_name, table_name))

    return table


def get_masses(ages, Zs, IMF, remnants=True):

    """Get a mass for each SSP in our set of spectra, for a given set of ages, IMFs and Zs

    Ages and Zs should be 1D arrays
    IMF should be a single value: either a float or a string like 'bi1.30' """

    Zs=np.array(Zs)
    ages=np.array(ages)

    #Test if our Zs and IMFs are arrays of strings, and if so convert them to floats
    if type(IMF) ==str:
        IMF=float(IMF.strip("bi"))
    if Zs.dtype.kind=="S":
        Zs=np.array([float(Z.replace("m", "-").replace("p", "+")) for Z in Zs])

    table=load_mass_table()

    #This rounding thing is a bit of a hack, but will work for our purposes as of Nov 2016.
    #The issue is that comparing floats between two lists is tricky
    age_mask=np.in1d(np.round(table['Age'], 1), np.round(ages, 1))
    Z_mask=np.in1d(np.round(table['[M/H]'], 2), np.round(Zs, 2))
    IMF_mask=(table['slope']==IMF)

    final_mask=age_mask & Z_mask & IMF_mask

    if remnants:
        masses=np.array(table['M(*+remn)'])[final_mask].reshape(len(Zs), len(ages))
    else:
        masses=np.array(table['M(*)'])[final_mask].reshape(len(Zs), len(ages))
    return masses

    #age_mask=np.where(table['Age'])




def load_Miles_M_L_table(IMF_type='bi', dir_name='/Data/stellarpops/Miles/Mass_files'):

    """Load the file with M/Ls as a table"""

    if IMF_type=='bi':
        table_name='ssp_phot_Padova00_BI_v10.txt'
    elif IMF_type=='uni' or IMF_type =='un':
        table_name='ssp_phot_Padova00_UN_v10.txt'


    from astropy.io import ascii

    table=ascii.read('{}/{}'.format(dir_name, table_name)) 

    return table

def get_M_L(IMF_type, band, age=14.1254, Z=0.00):

    """
    For a given IMF type, phtometric band, age and Z, get the M/L and IMF slope array ready for plotting
    """
    table=load_Miles_M_L_table(IMF_type)

    mask=np.where((table['Age']==age) & (table['[M/H]']==Z))

    assert len(mask[0])>0, "Seem to have a wrong age or metallicity: the mask has 0 elements"

    return table['slope'][mask], table['(M/L){}'.format(band)][mask]






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


    #import ipdb; ipdb.set_trace()

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











