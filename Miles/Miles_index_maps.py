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
    IMFs=sorted(spectra.keys())

    chab_index=np.where(np.array(IMFs)=='bi1.30')[0]


    array_indices_dict=M16.get_np_indices_for_params(IMFs=sorted(spectra.keys()), ages=ages, Zs=Zs)

    marker_sizes=np.linspace(30, 300, len(IMFs))

    if vary=='age':

        assert type(fixed_value)==str, "If we're varying age then the fixed value must be a metallicity!"
        fixed_index=array_indices_dict[fixed_value]


        for i, age in enumerate(ages):

            c=cm(1.0*i/ages.shape[0])

            label=r'Age={:.2} Gyr'.format(age)

            if alpha != 0.0:
                label=r'{}, $\alpha$={}'.format(label, alpha)
            if spectra[IMF].NaFe!=0.0:
                label=r'{}, [NaFe]={}'.format(label, spectra[IMF].NaFe)


            ax.plot(x_index_vals[:, fixed_index, i], y_index_vals[:, fixed_index, i], label=label, linewidth=2.0, c=c)
            ax.scatter(x_index_vals[:, fixed_index, i], y_index_vals[:, fixed_index, i], marker='o', s=marker_sizes, facecolors=c, linewidth=3.0, zorder=8)

            #Colour Chabrier different
            ax.scatter(x_index_vals[chab_index, fixed_index, i], y_index_vals[chab_index, fixed_index, i], marker='s', s=350, facecolors='b', linewidth=3.0, zorder=8)

        if make_grid==True:
            for i, IMF in enumerate(IMFs):
                ax.plot(x_index_vals[i, fixed_index, :], y_index_vals[i, fixed_index, :], c='0.75')

        ax.legend(loc=legend_loc, title="Z={}".format(fixed_value))


    elif vary=='Z':
        assert type(fixed_value)!=str, "If we're varying Z then the fixed value must be an age!"
        fixed_index=array_indices_dict[fixed_value]

        for i, Z in enumerate(Zs):
            c=cm(1.0*i/Zs.shape[0])


            label='Z={}'.format(Z)

            if alpha != 0.0:
                label=r'{}, $\alpha$={}'.format(label, alpha)
            if spectra[IMF].NaFe!=0.0:
                label=r'{}, [NaFe]={}'.format(label, spectra[IMF].NaFe)

            ax.plot(x_index_vals[:, i, fixed_index], y_index_vals[:, i, fixed_index], label=label, linewidth=2.0, c=c)
            ax.scatter(x_index_vals[:, i, fixed_index], y_index_vals[:, i, fixed_index], marker='o', s=marker_sizes, facecolors=c, linewidth=3.0, zorder=8)

            #Colour Chabrier different
            ax.scatter(x_index_vals[chab_index, i, fixed_index], y_index_vals[chab_index, i, fixed_index], marker='s', s=350, facecolors='b', linewidth=3.0, zorder=8)

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



    specs=M16.load_eMILES_spectra(NaFe=0.0)

    NaFep03_specs=M16.load_eMILES_spectra(NaFe=0.3)

    CvDinds=IT.getCvD12IndicesAir(verbose=False)



    
    def plot(figs, axs, GalName, data_fname, global_fname, marker, colour):

        fig1, fig2, fig3, fig4=figs
        ax1, ax2, ax3, ax4=axs

        fname1='MILES_NaI_MgI_fixed_age_{}.pdf'.format(GalName)
        fname2='MILES_FeH_TiO_fixed_age_{}.pdf'.format(GalName)



        ax1=index_index_map(ax1, specs, CvDinds.NaIsdss, CvDinds.MgI88, vary='Z', fixed_value=14.125)
        ax2=index_index_map(ax2, specs, CvDinds.FeH99, CvDinds.TiO89, vary='Z', fixed_value=14.125)

        ax1=index_index_map(ax1, specs, CvDinds.NaIsdss, CvDinds.MgI88, vary='Z', fixed_value=14.125, alpha=0.3, cm=plt.get_cmap('rainbow'))
        ax2=index_index_map(ax2, specs, CvDinds.FeH99, CvDinds.TiO89, vary='Z', fixed_value=14.125, alpha=0.3, cm=plt.get_cmap('rainbow'))

        ax1=index_index_map(ax1, NaFep03_specs, CvDinds.NaIsdss, CvDinds.MgI88, vary='Z', fixed_value=14.125, cm=plt.get_cmap('viridis'))
        ax2=index_index_map(ax2, NaFep03_specs, CvDinds.FeH99, CvDinds.TiO89,  vary='Z', fixed_value=14.125, cm=plt.get_cmap('viridis'))

        

        fname3='MILES_NaI_MgI_fixed_Z_{}.pdf'.format(GalName)
        fname4='MILES_FeH_TiO_fixed_Z_{}.pdf'.format(GalName)

        ax3=index_index_map(ax3, specs, CvDinds.NaIsdss, CvDinds.MgI88, vary='age', fixed_value='p0.00')
        ax4=index_index_map(ax4, specs, CvDinds.FeH99, CvDinds.TiO89, vary='age', fixed_value='p0.00')

        ax3=index_index_map(ax3, NaFep03_specs, CvDinds.NaIsdss, CvDinds.MgI88, vary='age', fixed_value='p0.00', cm=plt.get_cmap('viridis'))
        ax4=index_index_map(ax4, NaFep03_specs, CvDinds.FeH99, CvDinds.TiO89,  vary='age', fixed_value='p0.00', cm=plt.get_cmap('viridis'))


                
        R, Na, Na_err, CaT, CaT_err, FeH, FeH_err, Mg, Mg_err, TiO, TiO_err=NGC1277_data=np.genfromtxt(data_fname, unpack=True)
        for x, y, xerr, yerr in zip(Na, Mg, Na_err, Mg_err):
            ax1.errorbar(x, y, xerr=xerr, yerr=yerr, mfc=colour, fmt=marker, mew=3.0, mec='k', ms=10, zorder=9)
            ax3.errorbar(x, y, xerr=xerr, yerr=yerr, mfc=colour, fmt=marker, mew=3.0, mec='k', ms=10, zorder=9)

        for x, y, xerr, yerr in zip(FeH, TiO, FeH_err, TiO_err):
            ax2.errorbar(x, y, xerr=xerr, yerr=yerr, mfc=colour, fmt=marker, mew=3.0, mec='k', ms=10, zorder=9)
            ax4.errorbar(x, y, xerr=xerr, yerr=yerr, mfc=colour, fmt=marker, mew=3.0, mec='k', ms=10, zorder=9)



       
        R, Na, Na_err, CaT, CaT_err, FeH, FeH_err, Mg, Mg_err, TiO, TiO_err=NGC1277_data=np.genfromtxt(global_fname, unpack=True)
        ax1.errorbar(Na, Mg, xerr=Na_err, yerr=Mg_err, mfc='w', fmt=marker, mew=3.0, mec='k', ms=10, zorder=10)
        ax3.errorbar(Na, Mg, xerr=Na_err, yerr=Mg_err, mfc='w', fmt=marker, mew=3.0, mec='k', ms=10, zorder=10)
        ax2.errorbar(FeH, TiO, xerr=FeH_err, yerr=TiO_err, mfc='w', fmt=marker, mew=3.0, mec='k', ms=10, zorder=10)
        ax4.errorbar(FeH, TiO, xerr=FeH_err, yerr=TiO_err, mfc='w', fmt=marker, mew=3.0, mec='k', ms=10, zorder=10)


        #######################################

        ax1.set_xlim([0.2, 1.4])
        ax1.set_ylim([0.35, 0.65])

        ax3.set_xlim([0.2, 1.4])
        ax3.set_ylim([0.35, 0.65])

        ax2.set_xlim([0.0, 0.9])
        ax2.set_ylim([1.02, 1.09])

        ax4.set_xlim([0.0, 0.9])
        ax4.set_ylim([1.02, 1.09])

        #######################################

        for fig, fname in zip([fig1, fig2, fig3, fig4], [fname1, fname2, fname3, fname4]):
            fig.savefig(fname)


    fig1, ax1=plt.subplots(figsize=(16.5, 11.6))
    fig2, ax2=plt.subplots(figsize=(16.5, 11.6))
    fig3, ax3=plt.subplots(figsize=(16.5, 11.6))
    fig4, ax4=plt.subplots(figsize=(16.5, 11.6))

    fig5, ax5=plt.subplots(figsize=(16.5, 11.6))
    fig6, ax6=plt.subplots(figsize=(16.5, 11.6))
    fig7, ax7=plt.subplots(figsize=(16.5, 11.6))
    fig8, ax8=plt.subplots(figsize=(16.5, 11.6))

    N='/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/NGC1277.txt'
    I='/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/IC843.txt'

    GN='/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/GLOBAL_NGC1277.txt'
    IN='/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/GLOBAL_NGC1277.txt'


    plot([fig1, fig2, fig3, fig4], [ax1, ax2, ax3, ax4], 'NGC1277', N, GN, 'D', 'b')
    plot([fig5, fig6, fig7, fig8], [ax5, ax6, ax7, ax8], 'IC843', I, GN, 'o', 'r')

    plt.show()

