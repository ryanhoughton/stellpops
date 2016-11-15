import numpy as np 
import stellarpops.tools.specTools as s
import scipy.interpolate as si
import scipy.constants as const

import stellarpops.tools.Miles16tools as M16

import glob
import stellarpops.tools.indexTools as IT
import matplotlib.pyplot as plt

import seaborn as sns


def index_index_map(ax, spectra, index1, index2, param_list, param_name, fixed_value='p0.00', alpha=0.0, make_grid=True, c=plt.get_cmap('magma'), legend_loc='best'):

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
    all_ages=spectra[IMF].age[0, :]
    all_Zs=spectra[IMF].Z[:, 0]
    IMFs=sorted(spectra.keys())

    chab_index=np.where(np.array(IMFs)=='bi1.30')[0]


    array_indices_dict=M16.get_np_indices_for_params(IMFs=sorted(spectra.keys()), ages=all_ages, Zs=all_Zs)

    marker_sizes=np.linspace(30, 300, len(IMFs))


        
    fixed_index=array_indices_dict[fixed_value]

    good_indices=[]
    for cntr, thing in enumerate(param_list):

        i=array_indices_dict[thing]

        good_indices.append(i)

        if type(c)!=tuple:
            c=cm(1.0*cntr/len(param_list))
        


        #import pdb; pdb.set_trace()

        if param_name=='age':
            label=r'Age={0:.2f} Gyr'.format(thing)
            x=x_index_vals[:, fixed_index, i]
            y=y_index_vals[:, fixed_index, i]
            x_chab=x_index_vals[chab_index, fixed_index, i]
            y_chab=y_index_vals[chab_index, fixed_index, i]




        elif param_name=='Z':
            label='Z={}'.format(thing)
            

            x=x_index_vals[:, i, fixed_index]
            y=y_index_vals[:, i, fixed_index]

            x_chab=x_index_vals[chab_index, i, fixed_index]
            y_chab=y_index_vals[chab_index, i, fixed_index]

            

        else:
            raise NameError('Param name "{}" not understood'.format(param_name))


        if alpha != 0.0:
            label=r'{}, [$\alpha$/Fe]=+{}'.format(label, alpha)
        if spectra[IMF].NaFe!=0.0:
            label=r'{}, [Na/Fe]=+{}'.format(label, spectra[IMF].NaFe)


        ax.plot(x, y, label=label, linewidth=2.0, c=c)
        ax.scatter(x, y, marker='o', s=marker_sizes, facecolors=c, linewidth=3.0, zorder=8)

        #Colour Chabrier different
        ax.scatter(x_chab, y_chab, marker='s', s=350, facecolors='green', linewidth=3.0, zorder=8)

    if make_grid==True:
        for i, IMF in enumerate(IMFs):
            if param_name=='age':
                x_grid=x_index_vals[i, fixed_index, good_indices]
                y_grid=y_index_vals[i, fixed_index, good_indices]
            elif param_name=='Z':
                x_grid=x_index_vals[i, good_indices, fixed_index]
                y_grid=y_index_vals[i, good_indices, fixed_index]
            ax.plot(x_grid, y_grid, c='0.75')



    line_legend=ax.legend(loc=legend_loc, title="Age = {} Gyr".format(fixed_value), fontsize=20)
    ax.get_legend().get_title().set_fontsize('20')

    ax.add_artist(line_legend)


    """
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
                ax.plot(, c='0.75')

        ax.legend(loc=legend_loc, title="Age={} Gyr".format(fixed_value))


    else:
        raise 'Axis to vary not understood!'
    """



    
    return ax

def make_alpha_arrow(fig, ax, specs, alpha, index1, index2, arrow_start, **kwargs):


    all_ages=specs['bi1.30'].age[0, :]
    all_Zs=specs['bi1.30'].Z[:, 0]
    IMFs=sorted(specs.keys())

    array_indices_dict=M16.get_np_indices_for_params(IMFs=sorted(specs.keys()), ages=all_ages, Zs=all_Zs)



    arrow_xvals=np.empty(2)
    arrow_yvals=np.empty(2)


    arrow_xvals[0]=M16.alpha_corrected_index(specs['bi1.30'], index1, 200.0, alpha=alpha, verbose=False)[array_indices_dict['p0.00'], array_indices_dict[14.125]]
    arrow_yvals[0]=M16.alpha_corrected_index(specs['bi1.30'], index2, 200.0, alpha=alpha, verbose=False)[array_indices_dict['p0.00'], array_indices_dict[14.125]]

    arrow_xvals[1]=M16.alpha_corrected_index(specs['bi1.30'], index1, 200.0, alpha=0.0, verbose=False)[array_indices_dict['p0.00'], array_indices_dict[14.125]]
    arrow_yvals[1]=M16.alpha_corrected_index(specs['bi1.30'], index2, 200.0, alpha=0.0, verbose=False)[array_indices_dict['p0.00'], array_indices_dict[14.125]]

    dx=arrow_xvals[1]-arrow_xvals[0]
    dy=arrow_yvals[1]-arrow_yvals[0]

    x0, y0=arrow_start



    ax.annotate('', xytext=(x0, y0), xy=(x0+dx, y0+dy), xycoords='data', textcoords='data',arrowprops=dict(facecolor='black', shrink=0.0))
    ax.annotate(r'$\delta[\alpha /Fe]=-0.3$', xytext=(0, -10), xy=(x0+dx, y0+dy), xycoords='data', textcoords='offset pixels',horizontalalignment='center', verticalalignment='top', **kwargs)

    #ax.arrow(x0, y0, dx, dy, length_includes_head=True, **kwargs)

    return ax

def make_Z_arrow(fig, ax, specs, index1, index2, arrow_start, **kwargs):


    all_ages=specs['bi1.30'].age[0, :]
    all_Zs=specs['bi1.30'].Z[:, 0]
    IMFs=sorted(specs.keys())

    array_indices_dict=M16.get_np_indices_for_params(IMFs=sorted(specs.keys()), ages=all_ages, Zs=all_Zs)



    arrow_xvals=np.empty(2)
    arrow_yvals=np.empty(2)


    index_xvals=M16.alpha_corrected_index(specs['bi1.30'], index1, 200.0, alpha=0.0, verbose=False)
    index_yvals=M16.alpha_corrected_index(specs['bi1.30'], index2, 200.0, alpha=0.0, verbose=False)

    x=[ index_xvals[array_indices_dict['m0.40'], array_indices_dict[14.125]], index_xvals[array_indices_dict['p0.00'], array_indices_dict[14.125]], index_xvals[array_indices_dict['p0.22'], array_indices_dict[14.125]] ]
    y=[ index_yvals[array_indices_dict['m0.40'], array_indices_dict[14.125]], index_yvals[array_indices_dict['p0.00'], array_indices_dict[14.125]], index_yvals[array_indices_dict['p0.22'], array_indices_dict[14.125]] ]

    dx=x[0]-x[1]
    dy=y[0]-y[1]

    x0, y0=arrow_start


    ax.annotate('', xytext=(x0, y0), xy=(x0+dx, y0+dy), xycoords='data', textcoords='data',arrowprops=dict(facecolor='black', shrink=0.0))
    ax.annotate(r'$\delta[Z/H]=-0.22$', xytext=(0, -10), xy=(x0+dx, y0+dy), xycoords='data', textcoords='offset pixels', horizontalalignment='center', verticalalignment='top', **kwargs)

    return ax






def _plot_data(fig, ax, x, y, xerr, yerr, R, radial_cbar=True, marker='o'):

    if radial_cbar==True:
        plot=ax.scatter(x, y, s=250, c=R, marker=marker, linewidths=2.0, cmap='inferno', zorder=9, edgecolors='k')
        
        a,b,c =ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='', zorder=8, c='k')

        clb=fig.colorbar(plot, ax=ax)
        #max and min arcseconds from the NGC1277 data
        clb.set_clim(vmin=0.0, vmax=6.4)
        r_colour = clb.to_rgba(R)
        c[0].set_color(r_colour)
        c[1].set_color(r_colour)
        clb.ax.tick_params(labelsize=30.0)
        clb.set_label(r'R ($^{\prime\prime}$)', size=30.0)


    else:
        ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=marker, ms=24, c='w', mec='k', mew=2.0, ecolor='k', zorder=10)
        plot=None

    return ax, plot


    
if __name__=='__main__':



    specs=M16.load_eMILES_spectra(NaFe=0.0)

    NaFep03_specs=M16.load_eMILES_spectra(NaFe=0.3)
    NaFep06_specs=M16.load_eMILES_spectra(NaFe=0.6)
    NaFep09_specs=M16.load_eMILES_spectra(NaFe=0.9)



    #NaFep06_specs=M16.load_eMILES_spectra(NaFe=0.6)

    CvDinds=IT.getCvD12IndicesAir(verbose=False)


    #plot_Zs=['m0.40', 'p0.00', 'p0.22']
    plot_Zs=['p0.22']




    cm=plt.get_cmap('winter')

    

    def Na_Mg_plot(fig, ax, radial_fname, global_fname, GalName, cm=plt.get_cmap('winter')):

        ax=index_index_map(ax, NaFep03_specs, CvDinds.NaIsdss, CvDinds.MgI88, param_list=plot_Zs, param_name='Z', fixed_value=14.125, alpha=0.3, make_grid=True, c=cm(0.3), legend_loc='best')
        ax=index_index_map(ax, NaFep06_specs, CvDinds.NaIsdss, CvDinds.MgI88, param_list=plot_Zs, param_name='Z', fixed_value=14.125, alpha=0.3, make_grid=True, c=cm(0.6), legend_loc='best')
        ax=index_index_map(ax, NaFep09_specs, CvDinds.NaIsdss, CvDinds.MgI88, param_list=plot_Zs, param_name='Z', fixed_value=14.125, alpha=0.3, make_grid=True, c=cm(0.99), legend_loc='best')

        R, Na, Na_err, CaT, CaT_err, FeH, FeH_err, Mg, Mg_err, TiO, TiO_err=NGC1277_data=np.genfromtxt(radial_fname, unpack=True)
        ax, plot=_plot_data(fig, ax, Na, Mg, Na_err, Mg_err, R, radial_cbar=True)

        R, Na, Na_err, CaT, CaT_err, FeH, FeH_err, Mg, Mg_err, TiO, TiO_err=NGC1277_data=np.genfromtxt(global_fname, unpack=True)
        ax, _=_plot_data(fig, ax, Na, Mg, Na_err, Mg_err, R, radial_cbar=False)

        ax=make_alpha_arrow(fig, ax, NaFep03_specs, 0.3, CvDinds.NaIsdss, CvDinds.MgI88, arrow_start=(1.7, 0.32), **{'fontsize':20.0})

        ax=make_Z_arrow(fig, ax, NaFep03_specs, CvDinds.NaIsdss, CvDinds.MgI88, arrow_start=(1.7, 0.32), **{'fontsize':20.0})

        ax.set_xlabel(r'NaI$_{SDSS}$ (\AA)')
        ax.set_ylabel(r'MgI (\AA)')

        leg = ax.legend([plot], [GalName], loc='lower right', frameon=False, numpoints=1, scatterpoints=1, fontsize=20)

        return ax
        
#ax1=make_Z_line(fig1, ax1, NaFep03_specs,CvDinds.NaIsdss, CvDinds.MgI88, **{'ls':'dashed', 'marker':'*', 'ms':[10, 15, 20], 'c':'k', 'zorder':5})
    
    def FeH_TiO_plot(fig, ax, radial_fname, global_fname, GalName, cm=plt.get_cmap('winter')):
    
        #ax2=index_index_map(ax2, NaFep03_specs, CvDinds.FeH99, CvDinds.TiO89,  param_list=plot_Zs, param_name='Z', fixed_value=14.125, alpha=0.0, make_grid=True, c=cm(0.3), legend_loc='best')
        ax=index_index_map(ax, NaFep03_specs, CvDinds.FeH99, CvDinds.TiO89,  param_list=plot_Zs, param_name='Z', fixed_value=14.125, alpha=0.3, make_grid=True, c=cm(0.3), legend_loc='best')

        
        #ax=index_index_map(ax, NaFep06_specs, CvDinds.FeH99, CvDinds.TiO89,  param_list=plot_Zs, param_name='Z', fixed_value=14.125, alpha=0.0, make_grid=True, c=cm(0.6), legend_loc='best')
        ax=index_index_map(ax, NaFep06_specs, CvDinds.FeH99, CvDinds.TiO89,  param_list=plot_Zs, param_name='Z', fixed_value=14.125, alpha=0.3, make_grid=True, c=cm(0.6), legend_loc='best')

        
        #ax=index_index_map(ax, NaFep09_specs, CvDinds.FeH99, CvDinds.TiO89,  param_list=plot_Zs, param_name='Z', fixed_value=14.125, alpha=0.0, make_grid=True, c=cm(0.99), legend_loc='best')
        ax=index_index_map(ax, NaFep09_specs, CvDinds.FeH99, CvDinds.TiO89,  param_list=plot_Zs, param_name='Z', fixed_value=14.125, alpha=0.3, make_grid=True, c=cm(0.99), legend_loc='best')


        
        R, Na, Na_err, CaT, CaT_err, FeH, FeH_err, Mg, Mg_err, TiO, TiO_err=NGC1277_data=np.genfromtxt(radial_fname, unpack=True)
        ax, plot=_plot_data(fig, ax, FeH, TiO, FeH_err, TiO_err, R, radial_cbar=True)


        R, Na, Na_err, CaT, CaT_err, FeH, FeH_err, Mg, Mg_err, TiO, TiO_err=NGC1277_data=np.genfromtxt(global_fname, unpack=True)            
        ax, _=_plot_data(fig, ax, FeH, TiO, FeH_err, TiO_err, R, radial_cbar=False)

        ax=make_alpha_arrow(fig, ax, specs, 0.3, CvDinds.FeH99, CvDinds.TiO89, arrow_start=(0.1, 1.07), **{'fontsize':20.0})

        ax=make_Z_arrow(fig, ax, NaFep03_specs, CvDinds.FeH99, CvDinds.TiO89, arrow_start=(0.1, 1.07), **{'fontsize':20.0})    

        ax.set_xlabel(r'FeH (\AA)')
        ax.set_ylabel(r'TiO')

        leg = ax.legend([plot], [GalName], loc='lower right', frameon=False, numpoints=1, scatterpoints=1, fontsize=20)

  
        return ax

    plt.style.use('publication')

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    fig1, ax1=plt.subplots(figsize=(16.5, 11.6))
    fig2, ax2=plt.subplots(figsize=(16.5, 11.6))

    N='/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/NGC1277.txt'
    GN='/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/GLOBAL_NGC1277.txt'


    ax1=Na_Mg_plot(fig1, ax1, N, GN, GalName=r'NGC~1277')
    ax2=FeH_TiO_plot(fig2, ax2, N, GN, GalName=r'NGC~1277')

    fig1.savefig("NGC1277_NaI_Mg_MILES_index_index_map.pdf", bbox_inches='tight')
    fig2.savefig("NGC1277_FeH_TiO_MILES_index_index_map.pdf", bbox_inches='tight')






    fig3, ax3=plt.subplots(figsize=(16.5, 11.6))
    fig4, ax4=plt.subplots(figsize=(16.5, 11.6))


    I='/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/IC843.txt'
    GI='/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/GLOBAL_IC843.txt'


    ax3=Na_Mg_plot(fig3, ax3, I, GI, GalName=r'IC~843')
    ax4=FeH_TiO_plot(fig4, ax4, I, GI, GalName=r'IC~843')

    fig3.savefig("IC843_NaI_Mg_MILES_index_index_map.pdf", bbox_inches='tight')
    fig4.savefig("IC843_FeH_TiO_MILES_index_index_map.pdf", bbox_inches='tight')


    """
    for ax in [ax1, ax2]:
        ax.tick_params(axis='both', labelsize=tick_label_size)
    """

#    plt.show()

