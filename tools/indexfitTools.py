import numpy as np 
from stellarpops.tools import specTools as ST
#from stellarpops.tools import Miles16tools as M16
from stellarpops.tools import indexTools as IT
from stellarpops.tools import CD12tools as CT
import scipy.interpolate as si

import index_fit_settings as settings


def _CvD_measure_all_indices(base_spectra, var_elem_spectra, indices, sigma, corr_type='mult', outfile_string=None):

    """Take a dictionary of spectra and measure the index for each IMF, age and metallicity.
    Also, add elemental response functions to get index valus at varying [X/H]
    """
    
    imf_names=sorted(base_spectra.keys())
    ages=var_elem_spectra['Solar'].age[:, 0]
    Zs=var_elem_spectra['Solar'].Z[0, :]

    n_ages=len(ages)
    n_Zs=len(Zs)
    n_imfs=len(imf_names)

    #Load the element lists and steps from the settings file
    elem_steps, Na_elem_steps, positive_only_elem_steps=settings.get_steps()
    normal_elems, Na_elem, positive_only_elems=settings.get_elem_lists()    

    elem_param_dict=CT.CD16_get_np_indices_for_elems()
    base_param_dict=CT.CD16_get_np_indices_for_params()



    index_names=[ind['name'] for ind in indices]

    delta_index=np.empty((len(indices), len(normal_elems), len(elem_steps), n_ages, n_Zs, n_imfs))
    index_values=np.empty((len(indices), len(normal_elems), len(elem_steps), n_ages, n_Zs, n_imfs))

    delta_index_Na=np.empty((len(indices), 1, len(Na_elem_steps), n_ages, n_Zs, n_imfs))
    delta_index_positive=np.empty((len(indices), len(positive_only_elems), len(positive_only_elem_steps), n_ages, n_Zs, n_imfs))

    

    ###################################################################################################################################################################################################
    #Measure the indices with 0 enhancement as well as those with 'normal' enhancements

    print '\n\nMeasuring the indices\n\n'
    for ind, index in enumerate(indices):
        print 'Measuring {}'.format(index['name'])

        for i, elem in enumerate(normal_elems):

            for j, imf in enumerate(imf_names):

                for k, enhancement in enumerate(elem_steps):
                    if enhancement >=0.0:
                        e='{}+'.format(elem)
                    else:
                        e='{}-'.format(elem)
                        enhancement=np.abs(enhancement)

                    oldspec=CT.get_element_enhanced_spec(base_spectra[imf], e, enhancement=0.0, varelem_dict=var_elem_spectra, elem_param_dict=elem_param_dict, base_param_dict=base_param_dict)
                    newspec=CT.get_element_enhanced_spec(base_spectra[imf], e, enhancement=enhancement, varelem_dict=var_elem_spectra, elem_param_dict=elem_param_dict, base_param_dict=base_param_dict)

                    index_values[ind, i, k, :, :, j]=CT.CvD_cut_and_measure_index(oldspec, index, sigma, index_type='Cenarro', model_sigma=None)
                    if corr_type=='add':

                        
                        delta_index[ind, i, k, :, :, j]=CT.CvD_cut_and_measure_index(newspec, index, sigma, index_type='Cenarro', model_sigma=None)-index_values[ind, i, k, :, :, j]

                    else:
                        delta_index[ind, i, k, :, :, j]=(CT.CvD_cut_and_measure_index(newspec, index, sigma, index_type='Cenarro', model_sigma=None)-index_values[ind, i, k, :, :, j])/index_values[ind, i, k, :, :, j]

    ###################################################################################################################################################################################################


    ###################################################################################################################################################################################################
    #Measure Sodium Enhancement

    print 'Measuring Na from -0.45 dex to +1.0 dex'
    for ind, index in enumerate(indices):        

        for i, elem in enumerate(Na_elem):

            for j, imf in enumerate(imf_names):

                for k, enhancement in enumerate(Na_elem_steps):
                    if 0.0 <= enhancement <0.45:
                        e='{}+'.format(elem)
                        base_enhancement=0.3

                    elif 0.45 <= enhancement <0.75:
                        e='{}+0.6'.format(elem)
                        base_enhancement=0.6

                    elif 0.75 <= enhancement <1.0:
                        e='{}+0.9'.format(elem)
                        base_enhancement=0.9

                    else:
                        e='{}-'.format(elem)
                        enhancement=np.abs(enhancement)
                        base_enhancement=0.3

                    oldspec=CT.get_element_enhanced_spec(base_spectra[imf], e, enhancement=0.0, base_enhancement=base_enhancement, varelem_dict=var_elem_spectra, elem_param_dict=elem_param_dict, base_param_dict=base_param_dict)
                    newspec=CT.get_element_enhanced_spec(base_spectra[imf], e, enhancement=enhancement, base_enhancement=base_enhancement, varelem_dict=var_elem_spectra, elem_param_dict=elem_param_dict, base_param_dict=base_param_dict)

                    old_value=CT.CvD_cut_and_measure_index(oldspec, index, sigma, index_type='Cenarro', model_sigma=None)
                    if corr_type=='add':

                        
                        delta_index_Na[ind, i, k, :, :, j]=CT.CvD_cut_and_measure_index(newspec, index, sigma, index_type='Cenarro', model_sigma=None)-old_value

                    else:
                        delta_index_Na[ind, i, k, :, :, j]=(CT.CvD_cut_and_measure_index(newspec, index, sigma, index_type='Cenarro', model_sigma=None)-old_value)/old_value

    ###################################################################################################################################################################################################


    ###################################################################################################################################################################################################
    #Measure the elemnts which we only increase

    
    for ind, index in enumerate(indices):
        print 'Measuring {}'.format(index['name'])        

        for i, elem in enumerate(positive_only_elems):

            for j, imf in enumerate(imf_names):

                for k, enhancement in enumerate(positive_only_elem_steps):
                    
                    

                    oldspec=CT.get_element_enhanced_spec(base_spectra[imf], elem, enhancement=0.0, base_enhancement=base_enhancement, varelem_dict=var_elem_spectra, elem_param_dict=elem_param_dict, base_param_dict=base_param_dict)
                    newspec=CT.get_element_enhanced_spec(base_spectra[imf], elem, enhancement=enhancement, base_enhancement=base_enhancement, varelem_dict=var_elem_spectra, elem_param_dict=elem_param_dict, base_param_dict=base_param_dict)

                    old_value=CT.CvD_cut_and_measure_index(oldspec, index, sigma, index_type='Cenarro', model_sigma=None)
                    if corr_type=='add':

                        
                        delta_index_positive[ind, i, k, :, :, j]=CT.CvD_cut_and_measure_index(newspec, index, sigma, index_type='Cenarro', model_sigma=None)-old_value

                    else:
                        delta_index_positive[ind, i, k, :, :, j]=(CT.CvD_cut_and_measure_index(newspec, index, sigma, index_type='Cenarro', model_sigma=None)-old_value)/old_value

    ###################################################################################################################################################################################################

    #Cut down the index_values file (since all the values along the 1st and second axes are the same anyway)
    index_values=index_values[:, -1, -1, :, :, :]


    #Save the measurements as a fits file
    from astropy.io import fits

    hdu1=fits.PrimaryHDU()

    hdu2=fits.ImageHDU()
    hdu3=fits.ImageHDU()
    hdu4=fits.ImageHDU()

    hdu1.data=index_values
    hdu1.header['COMMENT']='measurements of {} at {} km/s with no element variations'.format(index_names, sigma)

    hdu2.data=delta_index
    hdu2.header['COMMENT']='Index measurements of {} at {} km/s with enhancements in {}'.format(index_names, sigma, normal_elems)

    hdu3.data=delta_index_Na
    hdu3.header['COMMENT']='Index measurements of {} at {} km/s with enhancements in Na'.format(index_names, sigma)

    hdu4.data=delta_index_positive
    hdu4.header['COMMENT']='Index measurements of {} at {} km/s with enhancements in {}'.format(index_names, sigma, positive_only_elems)

    hdulist=fits.HDUList([hdu1, hdu2, hdu3, hdu4])


    if outfile_string is None:
        if corr_type=='add':
            fname='model_measurements_additive_corr_sigma_{}.fits'.format(sigma)
        elif corr_type=='mult':
            fname='model_measurements_multiplicative_corr_sigma_{}.fits'.format(sigma)


    else:
        if corr_type=='add':
            fname='model_measurements_additive_corr_sigma_{}_{}.fits'.format(sigma, outfile_string)
        elif corr_type=='mult':
            fname='model_measurements_multiplicative_corr_sigma_{}.fits'.format(sigma, outfile_string)

    try:
        hdulist.writeto(fname, clobber=True)
    except:
        import pdb; pdb.set_trace()
    return 0

#########################################################################################################

def get_correction(interpolator, inds, elems, abunds, age, Z, imf):

    #The interpolator expects a list of 6 numbers. Meshgrid the two arrays which are of different lengths
    # (the indices and the number of elements to enhance) and then create lists of ages, Zs and IMFs of the correct
    # shapes. Then do the same for the abundances. Stack together and pass to the interpolator object!

    points = np.meshgrid(inds, elems, indexing='ij')
    flat = np.array([m.flatten() for m in points])
    #flat is now an array of points of shape 2, len(indices)*len(elems)
    #make arrays of the other variables of length len(indices)*len(elems)
    ages=np.ones_like(points[0])*age
    Zs=np.ones_like(points[0])*Z
    imfs=np.ones_like(points[0])*imf

    #Get the correct abundance for each element- luckily we can index the abundance array by the integer element array
    abunds=abunds[points[1]]

    #Stack together
    xi=np.vstack((flat, abunds.ravel(), ages.ravel(), Zs.ravel(), imfs.ravel()))
    #Do the interpolation
    out_array = interpolator(xi.T)
    #reshape everything to be (len(indices), len(elements))
    result = out_array.reshape(*points[0].shape)

    return result




def lnlike(theta, data, errors, interpolators, method='linear', corr_type='mult'):
        

    index_interpolator, general_delta_index_interpolator, Na_delta_index_interpolator, positive_delta_index_interpolator=interpolators

    #Get the correct values from theta
    Na_abundance=theta[0]
    general_abundances=theta[1:8]
    positive_abundances=theta[8:17]
    age, Z, imf=theta[17:]

    #Make arrays for each element
    all_inds=np.arange(len(data))
    general_elems=np.arange(len(general_abundances))
    positive_elems=np.arange(len(positive_abundances))

    #Get the corrections from the interpolators
    general_correction=get_correction(general_delta_index_interpolator, all_inds, general_elems, general_abundances, age, Z, imf)
    positive_correction=get_correction(positive_delta_index_interpolator, all_inds, positive_elems, positive_abundances, age, Z, imf)
    Na_correction=np.atleast_2d(Na_delta_index_interpolator((all_inds, Na_abundance, age, Z, imf))).T

    correction=np.hstack((general_correction, positive_correction, Na_correction))
    correction=np.sum(correction, axis=1)

    if corr_type=='add':
        model=index_interpolator((all_inds, age, Z, imf), method=method)+correction
    elif corr_type=='mult':
        model=index_interpolator((all_inds, age, Z, imf), method=method)*(1.0+correction)
    old=index_interpolator((all_inds, age, Z, imf), method=method)

    chisq=(data-model)/errors
    chi2=chisq*chisq
    
    return -0.5*np.sum(chi2), model



def lnprior(theta):

    Na_abundance=theta[0]
    general_abundances=theta[1:8]
    positive_abundances=theta[8:17]
    age, Z, imf=theta[17:]

    

    if np.all(general_abundances>=-0.45) and np.all(general_abundances<=0.45)  and np.all(positive_abundances>=0.0) and np.all(positive_abundances<=0.45) and -0.45 <= Na_abundance <= 1.0 and 1.0 < age < 13.0 and -1.5 < Z < 0.4 and 0.0 < imf <3.5:
        return 0.0
    return -np.inf



def lnprob(theta, data, errors, interpolators, method='linear', corr_type='mult'):


    lp=lnprior(theta)

    if np.isfinite(lp):
        ll, model=lnlike(theta, data, errors, interpolators, method, corr_type)
        return ll+lp, model

    return -np.inf, None





#########################################################################################################


def load_measure_NGC1277_WHT_indices(filename, indices):


    import os
    filename=os.path.expanduser(filename)

    lam, flux, err=np.genfromtxt(filename, unpack=True)

    WHT_spectrum=ST.spectrum(lam=lam, lamspec=flux, wavesyst='air', errlamspec=err)

    data=np.empty(len(indices))
    errors=np.empty(len(indices))

    for i, index in enumerate(indices):
        indval=ST.clip_spec_measure_index(WHT_spectrum, index)

        data[i]=indval[0]
        errors[i]=indval[1]

    return data, errors

#########################################################################################################



#########################################################################################################
def load_measure_NGC1277_LRIS_indices(filename, indices):

    import os
    filename=os.path.expanduser(filename)
    
    lam, flux, error, _, _=np.genfromtxt(filename, unpack=True)

    lamdas=lam/1.017044
    spectrum=ST.spectrum(lam=lamdas, lamspec=flux, wavesyst='vac', errlamspec=error)

    data=np.empty(len(indices))
    errors=np.empty(len(indices))



    for i, index in enumerate(indices):
        indval=ST.clip_spec_measure_index(spectrum, index)

        data[i]=indval[0]
        errors[i]=indval[1]

    #Make the TiO index have a small error. FIXME!
    errors[-2]=0.005

    # #indices- NaI, FeH, MgI, TiO89  ### Hbeta, Mgb, Fe52, Fe53
    # data=np.array([1.03586, 0.37062, 0.54611, 1.07290])#, 1.222, 4.704, 2.799, 2.122])


    return data, errors
#########################################################################################################

#########################################################################################################
def get_all_optical_NIR_inds(Lick_filename="~/z/Data/stellarpops/index_definitions/lickIndicesAir.txt", CvD_filename="~/z/Data/stellarpops/index_definitions/CvDIndicesVac.txt"):

    """Get all the indices from 3900A to 10000A.
    """

    import os
    Lick_filename=os.path.expanduser(Lick_filename)
    CvD_filename=os.path.expanduser(CvD_filename)

    Lick_Inds=IT.getLickIndicesVac(filename=Lick_filename, verbose=False)
    CvD_Inds=IT.getCvD12IndicesVac(filename=CvD_filename, verbose=False)

    indices=[CvD_Inds.CaII39, Lick_Inds.CN_1, Lick_Inds.CN_2, Lick_Inds.Ca4227, Lick_Inds.G4300, Lick_Inds.Fe4383, Lick_Inds.Ca4455, Lick_Inds.Fe4531, Lick_Inds.Fe4668, Lick_Inds.H_beta, Lick_Inds.Fe5015, Lick_Inds.Mg_1, 
Lick_Inds.Mg_2, Lick_Inds.Mg_b, Lick_Inds.Fe5270, Lick_Inds.Fe5335, Lick_Inds.Fe5406, Lick_Inds.Na_D,  Lick_Inds.Fe5709, Lick_Inds.Fe5782, Lick_Inds.TiO_1, Lick_Inds.TiO_2, 
CvD_Inds.NaIsdss, CvD_Inds.CaII86_1, CvD_Inds.CaII86_2, CvD_Inds.CaII86_3, CvD_Inds.MgI88, CvD_Inds.TiO89, CvD_Inds.FeH99]

    return indices
#########################################################################################################

#########################################################################################################
def get_NGC1277_WHT_indices(Lick_filename="~/z/Data/stellarpops/index_definitions/lickIndicesAir.txt", CvD_filename="~/z/Data/stellarpops/index_definitions/CvDIndicesVac.txt"):

    import os
    Lick_filename=os.path.expanduser(Lick_filename)
    CvD_filename=os.path.expanduser(CvD_filename)

    Lick_Inds=IT.getLickIndicesVac(filename=Lick_filename, verbose=False)
    CvD_Inds=IT.getCvD12IndicesVac(filename=CvD_filename, verbose=False)


    #WHT data indices
    indices=[CvD_Inds.CaII39, Lick_Inds.CN_1, Lick_Inds.CN_2, Lick_Inds.Ca4227, Lick_Inds.G4300, Lick_Inds.Fe4383, Lick_Inds.Ca4455, Lick_Inds.Fe4531, Lick_Inds.Fe4668, Lick_Inds.H_beta, Lick_Inds.Fe5015, Lick_Inds.Mg_1, 
Lick_Inds.Mg_2, Lick_Inds.Mg_b, Lick_Inds.Fe5270, Lick_Inds.Fe5335, Lick_Inds.Fe5406, Lick_Inds.Na_D,  Lick_Inds.Fe5709, Lick_Inds.Fe5782, Lick_Inds.TiO_1, Lick_Inds.TiO_2, 
CvD_Inds.NaIsdss, CvD_Inds.CaII86_1, CvD_Inds.CaII86_2, CvD_Inds.CaII86_3, CvD_Inds.MgI88]

    return indices
#########################################################################################################

#########################################################################################################
def get_NGC1277_LRIS_indices(Lick_filename="~/z/Data/stellarpops/index_definitions/lickIndicesAir.txt", CvD_filename="~/z/Data/stellarpops/index_definitions/CvDIndicesVac.txt"):

    import os
    Lick_filename=os.path.expanduser(Lick_filename)
    CvD_filename=os.path.expanduser(CvD_filename)

    Lick_Inds=IT.getLickIndicesVac(filename=Lick_filename, verbose=False)
    CvD_Inds=IT.getCvD12IndicesVac(filename=CvD_filename, verbose=False)

    indices=[CvD_Inds.CaII39, Lick_Inds.CN_1, Lick_Inds.CN_2, Lick_Inds.Ca4227, Lick_Inds.G4300, Lick_Inds.Fe4383, Lick_Inds.Ca4455, Lick_Inds.Fe4531, Lick_Inds.Fe4668, Lick_Inds.H_beta, Lick_Inds.Fe5015, Lick_Inds.Mg_1, 
Lick_Inds.Mg_2, Lick_Inds.Mg_b, Lick_Inds.Fe5270, Lick_Inds.Fe5335, Lick_Inds.Fe5406, CvD_Inds.NaIsdss, CvD_Inds.CaII86_1, CvD_Inds.CaII86_2, CvD_Inds.CaII86_3, CvD_Inds.MgI88, CvD_Inds.TiO89, CvD_Inds.FeH99] 
    #Lick.Hdelta_A, Lick.Hgamma_A, Lick.Hdelta_F, Lick.Hgamma_F

    return indices
 #########################################################################################################









