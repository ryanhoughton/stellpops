import numpy as np 
import stellarpops.tools.specTools as s
import scipy.interpolate as si
import scipy.constants as const

import stellarpops.tools.Miles16tools as M16

import glob
import stellarpops.tools.indexTools as IT
import matplotlib.pyplot as plt


def get_all_specs():

    specs=M16.load_eMILES_spectra(NaFe=0.0)

    NaFep03_specs=M16.load_eMILES_spectra(NaFe=0.3)
    NaFep06_specs=M16.load_eMILES_spectra(NaFe=0.6)
    NaFep09_specs=M16.load_eMILES_spectra(NaFe=0.9)

    return [specs, NaFep03_specs, NaFep06_specs, NaFep09_specs]

def make_interpolating_function(all_spectra, alpha, out_sigma=200.0):

    specs, NaFep03_specs, NaFep06_specs, NaFep09_specs=all_spectra

    #Get the Ages, Zs and IMFs we interpolate with
    all_ages=NaFep03_specs['bi1.30'].age[0, :]
    all_Zs=NaFep03_specs['bi1.30'].Z[:, 0]
    IMFs=sorted(NaFep03_specs.keys())

    full_array_indices_dict=M16.get_np_indices_for_params(IMFs=sorted(specs.keys()), ages=specs['bi1.30'].age[0, :], Zs=specs['bi1.30'].Z[:, 0])


    IMF_points=[float(IMF.strip('bi')) for IMF in IMFs]
    FeH_vals=np.empty_like(IMF_points)

    
    
    CvDinds=IT.getCvD12IndicesAir(verbose=False)
    #Need to interpolate age, Z, [Na/Fe], alpha




    #This is a bit tricker because the Na enhanced spectra have different params than the base set.
    #Get a dictionary for both
    #Na_array_indices_dict=M16.get_np_indices_for_params(IMFs=sorted(NaFep03_specs.keys()), ages=all_ages, Zs=all_Zs)

    #4 here refers to the number of Sodium enhanced spectra we have- Na/Fe=0.0, 0.3, 0.6, 0.9
    Na_IMF_cube=np.empty([len(IMFs), 4, len(all_Zs), len(all_ages)])

    #Loop through the Na enhanced spectra and put the FeH measurements (at enhanced alpha) in the array
    for i, Na_specs in enumerate([NaFep03_specs, NaFep06_specs, NaFep09_specs]):
        for j, IMF in enumerate(IMFs):
            #This is i+1 because we don't want to fill the Na/Fe=0.0 axis yet
            Na_IMF_cube[j, i+1, :, :]=M16.alpha_corrected_index(Na_specs[IMF], CvDinds.FeH99, out_sigma, alpha=alpha, verbose=False)



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

    #Fill the remaining array elements
    for j, IMF in enumerate(IMFs):
        Na_IMF_cube[j, 0, :, :]=M16.alpha_corrected_index(specs[IMF], CvDinds.FeH99, out_sigma, alpha=alpha, verbose=False)[Z_index.T, age_index.T]


    #Now perform basically the reverse process. Given the correct age, Z and alpha values from optical indices, find the FeH values from the models for each IMF


    interpolating_function = si.RegularGridInterpolator((IMF_points, Na_points, Z_points, age_points), Na_IMF_cube)

    return interpolating_function, IMF_points


def find_IMF_slope_from_FeH(interpolating_function, Gal_Na, Gal_Z, Gal_age, IMF_points, data, error):



    FeH_vals=np.zeros_like(IMF_points)

    for i, IMF in enumerate(IMF_points):


        FeH_vals[i]=interpolating_function([IMF, Gal_Na, Gal_Z, Gal_age])[0]

    #Use a 1D interpolation for the FeH values and the IMF points
    FeH_interp=si.interp1d(FeH_vals, IMF_points)
    
    frac_error=error/data

    final_IMF=FeH_interp(data)
    final_IMFs_error=FeH_interp(data)*(frac_error)

    return final_IMF, final_IMFs_error



all_spectra=get_all_specs()


NGC1277_FeH, NGC1277_FeH_err=np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/GLOBAL_NGC1277.txt')[5], np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/GLOBAL_NGC1277.txt')[6]

interpolating_function, IMF_points=make_interpolating_function(all_spectra, alpha=0.3, out_sigma=200.0)

print NGC1277_FeH, NGC1277_FeH_err

print find_IMF_slope_from_FeH(interpolating_function, Gal_Na=0.34, Gal_Z=0.21, Gal_age=13.5, IMF_points=IMF_points, data=NGC1277_FeH, error=NGC1277_FeH_err)



    





