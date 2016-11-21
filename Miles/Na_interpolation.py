import numpy as np 
import stellarpops.tools.specTools as s
import scipy.interpolate as si
import scipy.constants as const

import stellarpops.tools.Miles16tools as M16

import glob
import stellarpops.tools.indexTools as IT




def get_all_specs():

    specs=M16.load_eMILES_spectra(NaFe=0.0)

    NaFep03_specs=M16.load_eMILES_spectra(NaFe=0.3)
    NaFep06_specs=M16.load_eMILES_spectra(NaFe=0.6)
    NaFep09_specs=M16.load_eMILES_spectra(NaFe=0.9)

    return [specs, NaFep03_specs, NaFep06_specs, NaFep09_specs]

def make_interpolating_function(all_spectra, alpha, Fe, out_sigma=200.0):

    specs, NaFep03_specs, NaFep06_specs, NaFep09_specs=all_spectra

    #Get the Ages, Zs and IMFs we interpolate with
    all_ages=NaFep03_specs['bi1.30'].age[0, :]
    all_Zs=NaFep03_specs['bi1.30'].Z[:, 0]
    IMFs=sorted(NaFep03_specs.keys())

    full_array_indices_dict=M16.get_np_indices_for_params(IMFs=sorted(specs.keys()), ages=specs['bi1.30'].age[0, :], Zs=specs['bi1.30'].Z[:, 0])

    Z_points=[-0.4, 0.0, 0.22]
    Na_points=[0.0, 0.3, 0.6, 0.9]
    age_points=all_ages
    Na_points=[0.0, 0.3, 0.6, 0.9]



    
    
    CvDinds=IT.getCvD12IndicesAir(verbose=False)
    #Need to interpolate age, Z, [Na/Fe], alpha




    #This is a bit tricker because the Na enhanced spectra have different params than the base set.
    #Get a dictionary for both
    #Na_array_indices_dict=M16.get_np_indices_for_params(IMFs=sorted(NaFep03_specs.keys()), ages=all_ages, Zs=all_Zs)

    #4 here refers to the number of Sodium enhanced spectra we have- Na/Fe=0.0, 0.3, 0.6, 0.9
    cube=np.empty([len(Na_points), len(all_Zs), len(all_ages)])

    #Loop through the Na enhanced spectra and put the FeH measurements (at enhanced alpha) in the array
    for i, Na_specs in enumerate([NaFep03_specs, NaFep06_specs, NaFep09_specs]):
        if Fe !=0.0:
            fe_spec=M16.fe_enhanced_spec(Na_specs['bi1.30'], CvDinds.NaIsdss, out_sigma, Fe=Fe)
            cube[i+1, :, :]=M16.alpha_corrected_index(fe_spec, CvDinds.NaIsdss, out_sigma, alpha=alpha, verbose=False)
        else:
            cube[i+1,:, :]=M16.alpha_corrected_index(Na_specs['bi1.30'], CvDinds.NaIsdss, out_sigma, alpha=alpha, verbose=False)



    #Index the bigger 'specs' array to only pick out the ages and metallicities we want
    #To do this, select the ages and Zs in common between the base set and Na enhanced model.
    good_ages=[]
    good_Zs=[]

    for age in all_ages:    
        good_ages.append(full_array_indices_dict[age])

    for Z in all_Zs:
        good_Zs.append(full_array_indices_dict[Z])

    #Meshgrid to get the two index arrays we want
    Z_index, age_index=np.meshgrid(np.array(good_Zs), np.array(good_ages))

    #Fill the remaining array elements- those with Na/Fe=0.0

    if Fe !=0.0:
        fe_spec=M16.fe_enhanced_spec(specs['bi1.30'], CvDinds.NaIsdss, out_sigma, Fe=Fe)
        cube[0, :, :]=M16.alpha_corrected_index(fe_spec, CvDinds.NaIsdss, out_sigma, alpha=alpha, verbose=False)[Z_index.T, age_index.T]
    else:
        cube[0, :, :]=M16.alpha_corrected_index(specs['bi1.30'], CvDinds.NaIsdss, out_sigma, alpha=alpha, verbose=False)[Z_index.T, age_index.T]

    interpolating_function = si.RegularGridInterpolator((Na_points, Z_points, age_points), cube, method='linear', bounds_error=False, fill_value=None)

    return interpolating_function, Na_points


def find_Na_enhancement(interpolating_function, Gal_Z, Gal_age, Na_points, data, error, galname):

    Na_vals=np.zeros_like(Na_points)




    Z_points=[-0.4, 0.0, 0.22]
    Na_points=[0.0, 0.3, 0.6, 0.9]

    for i, Na in enumerate(Na_points):


        Na_vals[i]=interpolating_function([Na, Gal_Z, Gal_age])[0]

    #Use a 1D interpolation for the FeH values and the IMF points
    Na_interp=si.interp1d(Na_vals, Na_points, fill_value='extrapolate')
    
    frac_error=error/data
    final_Na=Na_interp(data)

    p_error=Na_interp(data+error)
    m_error=Na_interp(data-error)
    frac_error=Na_interp(data)*(frac_error)
    final_Na_error=frac_error



    print "Final Na enhancement for {} is {} pm {}".format(galname, final_Na, final_Na_error)

    return final_Na, final_Na_error





def Na_enhancement_wrapper(galname, all_spectra, sigma, age, alpha, Z, Na, Na_error):

    print "\n"
    print "Interpolating grid for {}".format(galname)
    print "params: alpha={}, age{}, Z={}, sigma={}".format(alpha, age, Z, sigma)
    interpolating_function, Na_points=make_interpolating_function(all_spectra, alpha=alpha, Fe=0.0, out_sigma=sigma)

    final_Na, final_Na_error=find_Na_enhancement(interpolating_function, Z, age, Na_points, Na, Na_error, galname)



    

    return final_Na, final_Na_error


if __name__=='__main__':


    #all_spectra=get_all_specs()


    NGC1277_Na, NGC1277_Na_err=np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/NGC1277.txt').T[1], np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/NGC1277.txt').T[2]
    IC843_Na, IC843_Na_err=np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/IC843.txt').T[1], np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/IC843.txt').T[2]

    NGC1277_results=np.zeros_like(NGC1277_Na)
    IC843_results=np.zeros_like(IC843_Na)

    for i, (n, e) in enumerate(zip(NGC1277_Na, NGC1277_Na_err)):

        final_Na, err=Na_enhancement_wrapper('NGC 1277', all_spectra, 200.0, 13.0, 0.3, 0.3, n, e)
        NGC1277_results[i]=final_Na

    for i, (n, e) in enumerate(zip(IC843_Na, IC843_Na_err)):
        final_Na, err=Na_enhancement_wrapper('IC843', all_spectra, 200.0, 10.0, 0.3, 0.2, n, e)

        IC843_results[i]=final_Na

    print NGC1277_results, IC843_results











