
import numpy as np
import matplotlib.pyplot as plt
import spectools as s

import scipy.interpolate as si
import scipy.constants as sc

import CvD12tools as cvd
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.colors import colorConverter

import seaborn as sns

import Miles_FeH_interpolation as MilesFeH

####FONTS
tfsize = 28
lfsize = 24
gensize = 20
majtsize = 14
mintsize = 6
twidth = 2
msize = 10
####


sns.set_context('poster')
#For saving the data:


filename="data_IMF_slope_sigma.txt"
with open(filename, 'w') as f:
    pass

def makeErrorBoxes(xdata,ydata,xerror,yerror,fc='r',ec='None',alpha=0.25):

    # Create list for all the error patches
    errorboxes = []

    # Loop over data points; create box from errors at each point
    for x_left,y_left,width,height in zip(xdata,ydata,xerror,yerror):
        
        rect = Rectangle((x_left,y_left),width,height)
        errorboxes.append(rect)

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes,facecolor=fc,alpha=alpha,edgecolor=ec, zorder=1)

    # Add collection to axes
    return pc


def make_arrow(ax, x, y0, y1, c='k'):

    dx=0.0
    dy=y1-y0

    ax.arrow(x, y0, dx, dy, head_width=4.0, head_length=0.05, fc=c, ec='k', length_includes_head=True)

    return 'done'
#------------#
#Plot of IMF v velocity dispersion - comparing our galaxies to literature results
#spiniello = np.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/SZThesis/Chapter_coma_plots/Spiniello2014_data.txt', delimiter=', ')
#ziel = np.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/SZThesis/Chapter_coma_plots/Zieleniewski2015-16_data.txt', delimiter=', ')

#Load up my FeH data for each galaxy
m31 = np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/IMF_Plots/FeH_index_values_mids_m31_1_30-01-15.txt', delimiter=',')
m32 = np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/IMF_Plots/FeH_index_values_mids_m32_30-01-15.txt', delimiter=',')

fehpos = 13
fehvpos = 14
gmp2921 = np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/IMF_Plots/Indices_v_radius_GMP2921_400kms_data_2016-4-13.txt', delimiter=',')
# gmp2921 = np.genfromtxt('GMP2921/GMP2921_global_spec_index_values_400kms_data_2016-4-13.txt', delimiter=',')
#gmp3329 = np.genfromtxt('GMP3329/Indices_v_radius_GMP3329_270kms_data_2016-4-13.txt', delimiter=',')
gmp3329 = np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/IMF_Plots/GMP3329_global_spec_index_values_270kms_data_2016-4-13.txt', delimiter=',')
# gmp4928 = np.genfromtxt('GMP4928/Indices_v_radius_GMP4928_270kms_data_2016-4-13.txt', delimiter=',')
gmp4928 = np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/IMF_Plots/GMP4928_global_spec_index_values_270kms_data_2016-4-13.txt', delimiter=',')
#gmp3367 = np.genfromtxt('GMP3367/GMP3367_global_spec_index_values_200kms_data_2016-4-14.txt', delimiter=',')
gmp3367 = np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/IMF_Plots/GMP3367_global_spec_index_values_185kms_data_2016-6-3.txt', delimiter=',')


#Data from SPV+2017
IC843_FeH, IC843_FeH_err=np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/GLOBAL_IC843.txt')[5], np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/GLOBAL_IC843.txt')[6]
IC843_FeH_var=IC843_FeH_err**2

NGC1277_FeH, NGC1277_FeH_err=np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/GLOBAL_NGC1277.txt')[5], np.genfromtxt('/Volumes/SPV_SWIFT/Science/Resolved_Gals_indices/Index_Index_Plots/GLOBAL_NGC1277.txt')[6]
NGC1277_FeH_var=NGC1277_FeH_err**2


m31feh = np.mean(m31[:2,1])
m31fehv = np.sum(m31[:2,2]**2)/4.

m32feh = np.mean(m32[:2,1])
m32fehv = np.sum(m32[:2,2]**2)/4.

g2921feh = np.mean(gmp2921[:2,fehpos])
g2921fehv = np.sum(gmp2921[:2,fehvpos])/4.
# g2921feh = gmp2921[8]
# g2921fehv = gmp2921[9]**2

# g3329feh = np.mean(gmp3329[:2,fehpos])
# g3329fehv = np.sum(gmp3329[:2,fehvpos])/4.
g3329feh = gmp3329[8]
g3329fehv = gmp3329[9]**2

# g4928feh = np.mean(gmp4928[:1,fehpos])
# g4928fehv = np.sum(gmp4928[:1,fehvpos])/1.
g4928feh = gmp4928[8]
g4928fehv = gmp4928[9]**2

g3367feh = gmp3367[8]
g3367fehv = gmp3367[9]**2

filts = s.loadJCFilters()

gals = ['M32', 'M31', 'NGC4889', 'NGC4874', 'NGC4839', 'NGC4873', 'NGC1277', 'IC843']
fehs = [m32feh, m31feh, g2921feh, g3329feh, g4928feh, g3367feh, NGC1277_FeH, IC843_FeH]
fehvs = [m32fehv, m31fehv, g2921fehv, g3329fehv, g4928fehv, g3367fehv, NGC1277_FeH_var, IC843_FeH_var]

all_spectra=MilesFeH.get_all_specs()
corrections=True

if corrections==True:
    #Alpha, iron and sodium corrections
    ages = np.array([4.0, 10.0, 13.5, 13.5, 13.5, 12.0, 13.5, 10.0])
    alphas = np.array([0.0, 0.2, 0.25, 0.25, 0.25, 0.25, 0.3, 0.3])
    fes = np.array([0.0, 0.1, 0.25, 0.02, 0.07, -0.05, 0.021, 0.021])
    nas = np.array([0.0, 0.3, 0.8, 0.1, 0.4, 0.1, 0.4, 0.2])
    Z_vals=fes+0.93*alphas
    resols = np.array([150.0, 200.0, 400.0, 269.0, 271.0, 185.0, 200.0, 200.0])
    plresols = np.array([65.0, 180.0, 380.0, 267.0, 278.0, 170.0, 410.0, 265.0])
    #savefilename="/Volumes/SPV_SWIFT/LaTeX/RadialGradientsPaper/Plots/IMF_Sigma_plot_WITH_response_functs.pdf"
else:
    # #No alpha, iron and sodium corrections
    ages = np.array([4.0, 10.0, 13.5, 13.5, 13.5, 12.0, 13.5, 10.0])
    alphas = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    fes = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    nas = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    resols = np.array([150.0, 200.0, 400.0, 269.0, 271.0, 185.0, 200.0, 200.0])
    Z_vals=fes+0.93*alphas
    plresols = np.array([65.0, 180.0, 380.0, 267.0, 278.0, 170.0, 410.0, 265.0])
    #savefilename="/Volumes/SPV_SWIFT/LaTeX/RadialGradientsPaper/Plots/IMF_Sigma_plot_WITHOUT_response_functs.pdf"


savefilename="test.pdf"







ind1 = 'FeH'
galimfs = []
galimferrs = []
galm2ls = []
galm2lerrs = []
interpfunctions = []
m2linterpfunctions = []



for i, (age, alpha, Z, fe, na, sigma, feh, feh_var, central_sigma, galname) in enumerate(zip(ages, alphas, Z_vals, fes, nas, resols, fehs, fehvs, plresols, gals)):

    #Hack- set Z=0 and use the Fe correction instead
    fe=0.0
    imf, imferr=MilesFeH.IMF_slope_wrapper(galname, all_spectra, sigma, age, alpha, Z, fe, na, feh, np.sqrt(feh_var))

    galimfs.append(imf)
    galimferrs.append(imferr)

    #galm2ls.append(m2l)
    #galm2lerrs.append(m2l_err)





plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': gensize})
fig, ax=plt.subplots(figsize=(16.5, 11.6))
#plt.subplots_adjust(left=0.14, bottom=0.12)



"""
NGC1277_alpha=make_arrow(ax, 423.0, NGC1277_no_response, NGC1277_with_alpha, c='g')
NGC1277_na=make_arrow(ax, 425.0, NGC1277_no_response, NGC1277_with_na, c='b')
NGC1277_fe=make_arrow(ax, 427.0, NGC1277_no_response, NGC1277_with_fe, c='r')


IC843_alpha=make_arrow(ax, 278.0, IC843_no_response, IC843_with_alpha, c='g')
IC843_na=make_arrow(ax, 280.0, IC843_no_response, IC843_with_na, c='b')
IC843_fe=make_arrow(ax, 282.0, IC843_no_response, IC843_with_fe, c='r')



#Plot the arrows for each response function

"""



#---#


#spin_line, = ax.plot(sigs, s14, 'r-', lw=2.0)
#f13_line, = ax.plot(sigs, f13, 'b-', lw=2.0)
#lb13_line, = ax.plot(sigs, lb13, 'g-', lw=2.0)


ziel_points = ax.errorbar(plresols[:-2], galimfs[:-2], yerr=galimferrs[:-2], fmt='s', ms=15, linewidth=1.0, alpha=1.0, zorder=10, c='k')
SPV_points = ax.errorbar(plresols[-2:], galimfs[-2:], yerr=galimferrs[-2:], fmt='D', ms=15, linewidth=1.0,  alpha=1.0, zorder=10, c='r')#cmap=plt.cm.copper_r, c=ages
#eboxes = ax.add_collection(boxes)


ax.set_xlabel('$\sigma$ [km/s]', fontsize=tfsize)
ax.set_ylabel('IMF slope $x$', fontsize=tfsize)
plt.xticks(size=lfsize)
plt.yticks(size=lfsize)
plt.minorticks_on()
ax.tick_params(axis='both', direction='in', length=majtsize,
              width=1.5, which='major')
ax.tick_params(axis='both', direction='in', length=mintsize,
              width=1.5, which='minor')

ax.set_xlim([50.0, 450.0])
ax.set_ylim([0.3, 3.5])

ax.plot([50.,450.], [1.3,1.3], 'k--', lw=1.5)
#ax.plot([50.0,450.0], [2.35, 2.35], 'k--', lw=1.5)
#ax.annotate('Salpeter', xy=(325,2.35),xycoords='data',xytext=(0,5),textcoords='offset points', color='k', fontsize=gensize)
for i in xrange(len(gals)):
    if gals[i] == 'NGC4873' or gals[i] == 'IC843' or gals[i] == 'NGC4889':
        ax.annotate(gals[i], xy=(plresols[i],galimfs[i]),xycoords='data',xytext=(-10,10),textcoords='offset points', color='k', fontsize=gensize, horizontalalignment='right', verticalalignment='bottom')
    else:
        ax.annotate(gals[i], xy=(plresols[i],galimfs[i]),xycoords='data',xytext=(10,10),textcoords='offset points', color='k', fontsize=gensize, horizontalalignment='left', verticalalignment='bottom')

# leg = ax.legend([spin_line, f13_line, lb13_line, ziel_points, SPV_points],
# ['Spiniello+14 relation', 'Ferreras+13 relation', 'La Barbera+13 relation', 'Zieleniewski+16', "This work"], loc='upper left',
# frameon=False, borderaxespad=0.3, labelspacing=0.5, numpoints=1, scatterpoints=1)
#plt.gca().add_artist(leg2)



#print gals, galm2ls, galm2lerrs
plt.tight_layout()
#plt.savefig(savefilename, bbox_inches='tight')
plt.show()
