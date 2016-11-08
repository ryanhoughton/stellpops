import numpy as np 
import stellarpops.tools.specTools as s
import scipy.interpolate as si
import scipy.constants as const

import stellarpops.tools.Miles16tools as M16

import glob
import stellarpops.tools.indexTools as IT
import matplotlib.pyplot as plt



def index_index_map(ax, spectra, index1, index2, vary='age', fixed_value='p0.00', make_grid=True, cm=plt.get_cmap('magma'), legend_loc='best'):

    """
    Plot an index-index map of two indices as a function of IMF at fixed age and Z. Then loop through another variable
    defined by vary='age' or vary='Z' to make a grid. 
    """

    x_index_vals=[]
    y_index_vals=[]  

    for IMF in sorted(spectra.keys()):

        x_index_vals.append(M16.cut_and_measure_index(spectra[IMF], index1, 200.0))
        y_index_vals.append(M16.cut_and_measure_index(spectra[IMF], index2, 200.0))

    x_index_vals=np.array(x_index_vals)
    y_index_vals=np.array(y_index_vals)

    ages=spectra[IMF].age
    Zs=spectra[IMF].Z
    IMFs=sorted(specs.keys())

    marker_sizes=np.linspace(30, 300, len(IMFs))

    if vary=='age':

        assert type(fixed_value)==str, "If we're varying age then the fixed value must be a metallicity!"
        fixed_index=array_indices_dict[fixed_value]

        for i, age in enumerate(ages):

            c=cm(1.0*i/len(ages))




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
            c=cm(1.0*i/len(Zs))

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


    fig, axs=plt.subplots(nrows=3, ncols=1)


    axs[0]=index_index_map(axs[0], specs, CvDinds.NaIsdss, CvDinds.MgI88, vary='Z', fixed_value=10.0)
    axs[1]=index_index_map(axs[1], specs, CvDinds.FeH99, CvDinds.TiO89, vary='Z', fixed_value=10.0)
    axs[2]=index_index_map(axs[2], specs, CvDinds.FeH99, CvDinds.CaII86, vary='Z', fixed_value=10.0)

