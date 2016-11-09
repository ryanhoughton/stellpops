import numpy as np 
import stellarpops.tools.specTools as s
import scipy.interpolate as si
import scipy.constants as const

import stellarpops.tools.Miles16tools as M16

import glob
import stellarpops.tools.indexTools as IT
import matplotlib.pyplot as plt

import seaborn as sns


def index_index_map(ax, spectra, index1, index2, vary='age', fixed_value='p0.00', alpha=0.0, make_grid=True, cm=plt.get_cmap('magma'), legend_loc='best'):

    """
    Plot an index-index map of two indices as a function of IMF at fixed age and Z. Then loop through another variable
    defined by vary='age' or vary='Z' to make a grid. 
    """

    x_index_vals=[]
    y_index_vals=[]  

    for IMF in sorted(spectra.keys()):

        if alpha==0.0:

            x_index_vals.append(M16.cut_and_measure_index(spectra[IMF], index1, 200.0))
            y_index_vals.append(M16.cut_and_measure_index(spectra[IMF], index2, 200.0))
        else:
            x_index_vals.append(M16.alpha_corrected_index(spectra[IMF], index1, 200.0, alpha=alpha, verbose=False))
            y_index_vals.append(M16.alpha_corrected_index(spectra[IMF], index2, 200.0, alpha=alpha, verbose=False))


    x_index_vals=np.array(x_index_vals)
    y_index_vals=np.array(y_index_vals)

    #Horrible hack to get the right ages and Zs
    #I dont think having an Nd array of ages and Zs is a good idea!
    ages=spectra[IMF].age[0, :]
    Zs=spectra[IMF].Z[:, 0]
    IMFs=sorted(specs.keys())

    marker_sizes=np.linspace(30, 300, len(IMFs))

    if vary=='age':

        assert type(fixed_value)==str, "If we're varying age then the fixed value must be a metallicity!"
        fixed_index=array_indices_dict[fixed_value]

        for i, age in enumerate(ages):

            c=cm(1.0*i/ages.shape[0])




            ax.plot(x_index_vals[:, fixed_index, i], y_index_vals[:, fixed_index, i], label='Age={:.2} Gyr'.format(age), linewidth=2.0, c=c)
            ax.scatter(x_index_vals[:, fixed_index, i], y_index_vals[:, fixed_index, i], marker='o', s=marker_sizes, facecolors=c, linewidth=3.0, zorder=10)

        if make_grid==True:
            for i, IMF in enumerate(IMFs):
                ax.plot(x_index_vals[i, fixed_index, :], y_index_vals[i, fixed_index, :], c='0.75')

        ax.legend(loc=legend_loc, title="Z={}".format(fixed_value))


    elif vary=='Z':
        assert type(fixed_value)!=str, "If we're varying Z then the fixed value must be an age!"
        fixed_index=array_indices_dict[fixed_value]

        for i, Z in enumerate(Zs):
            c=cm(1.0*i/Zs.shape[0])

            ax.plot(x_index_vals[:, i, fixed_index], y_index_vals[:, i, fixed_index], label='Z={}'.format(Z), linewidth=2.0, c=c)
            ax.scatter(x_index_vals[:, i, fixed_index], y_index_vals[:, i, fixed_index], marker='o', s=marker_sizes, facecolors=c, linewidth=3.0, zorder=10)

        if make_grid==True:
            for i, IMF in enumerate(IMFs):
                ax.plot(x_index_vals[i, :, fixed_index], y_index_vals[i, :, fixed_index], c='0.75')

        ax.legend(loc=legend_loc, title="Age={} Gyr".format(fixed_value))


    else:
        raise 'Axis to vary not understood!'

    ax.set_xlabel('{}'.format(index1['name']))
    ax.set_ylabel('{}'.format(index2['name']))

    
    return ax


    
if __name__=='__main__':
    #Used for testing things



    specs=M16.load_eMILES_spectra()
    CvDinds=IT.getCvD12IndicesAir(verbose=False)

    array_indices_dict=M16.get_np_indices_for_params()


    fig, ax1=plt.subplots()
    fig2, ax2=plt.subplots()

    ax1=index_index_map(ax1, specs, CvDinds.NaIsdss, CvDinds.MgI88, vary='Z', fixed_value=14.125)
    ax2=index_index_map(ax2, specs, CvDinds.FeH99, CvDinds.TiO89, vary='Z', fixed_value=14.125)

    ax1=index_index_map(ax1, specs, CvDinds.NaIsdss, CvDinds.MgI88, vary='Z', fixed_value=14.125, alpha=0.3, cm=plt.get_cmap('viridis'))
    ax2=index_index_map(ax2, specs, CvDinds.FeH99, CvDinds.TiO89,  vary='Z', fixed_value=14.125, alpha=0.3, cm=plt.get_cmap('viridis'))
    #ax2=index_index_map(ax2, specs, CvDinds.FeH99, CvDinds.TiO89,  vary='Z', fixed_value=14.125, alpha=0.7, cm=plt.get_cmap('viridis'))
    #axs[2]=index_index_map(axs[2], specs, CvDinds.FeH99, CvDinds.CaII86, vary='Z', fixed_value=10.0)

    plt.show()

