from __future__ import print_function

import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits 

from astropy.io import ascii
from scipy import ndimage
import numpy as np
from time import clock


from scipy.interpolate import interp1d

from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
import pdb

#from spectools import *
#from CvD12tools import *
#from python_utils import sky_shift
import matplotlib.pyplot as plt
from sam_python.continuum import fit_continuum



def ppxf_fit_spec(galaxy, noise, templates, fit_range, lamRange2, logLam1, logLam2, velscale, start, sky=None, plot=True, moments=4, degree=-1, mdegree=4, goodpixels=None):


    """Use pPXF to fit kinematics to a spectrum. This function masks the templates so they're the appropriate lengths (a bit longer than the galaxy spectrum) then
    runs pPXF


    inputs:

    galaxy- a full length galaxy spectrum (i.e 4112 pixels for SWIFT)
    noise- noise spectrum, same length as galaxy
    templates- array of templates
    lamRange2- wavelength range of the templates, pre masking
    logLam1- log rebinned wavelength array for galaxy
    logLam2- log rebinned wavelength array for templates
    fit_range- range over which you want to fit the galaxy
    start- starting guess for V and Sigma (or h3, h4, etc if you want higher moments)
    sky- optionally fit a sky spectra at the same time. Set to None if you're not using
  


    outputs:

    pp- the ppxf class
    sky- the input sky spectrum, but masked. If sky=None, then this is None
    galaxy- the input galaxy spectrum, but masked
    noise- the input noise spectrum, but masked



    """
    
    #We need to pad the templates for ~100 angstroms each side of the galaxy spectrum.
    

    lower_lam, upper_lam=np.array(fit_range)

    

    #print(upper_lam) 

    #print(lamRange2)


    pad=300

    #make the template mask
    tmask = np.where( (logLam2>=np.log((lower_lam-pad))) & (logLam2<=np.log((upper_lam+pad))))[0]



    #mask the  wavelength range and the templates
    logLam2=logLam2[tmask]
    templates=templates[tmask, :]


    #make the mask for the galaxy spectrum
    #mask=np.where((logLam1 >= np.log(lower_lam)) & (logLam1 <= np.log(upper_lam)))[0]
    

    #mask the galaxy, sky and variance, plus the wavelengh ranges
    """
    logLam1=logLam1[mask]


    galaxy=galaxy[mask]



    if SKY_FLAG:

        sky=sky[mask, :]

      
    noise=noise[mask]
    """





    #################################################################################

    # The galaxy and the template spectra do not have the same starting wavelength.
    # For this reason an extra velocity shift DV has to be applied to the template
    # to fit the galaxy spectrum. We remove this artificial shift by using the
    # keyword VSYST in the call to PPXF below, so that all velocities are
    # measured with respect to DV. This assume the redshift is negligible.
    # In the case of a high-redshift galaxy one should de-redshift its
    # wavelength to the rest frame before using the line below (see above).
    #
    c = 299792.458
    

    #print("logam2[0] is {}, loglam1[0] is {}, logLam1/(1+z) is {}".format(logLam2[0], logLam1[0],logLam1[0]/(1+z)))
    dv = (logLam2[0]-logLam1[0])*c # km/s



    vel, sigma=start



    z_ppxf = np.exp(vel/c) - 1   # Relation between velocity and redshift in pPXF


    if goodpixels is None:
        goodpixels = util.determine_goodpixels(logLam1, lamRange2, z_ppxf)


    print("#########################################################")
    print("Velocity shift DV is {}".format(dv))
    print("The templates are shape {}".format(np.shape(templates)))
    print("The galaxy is shape {}".format(np.shape(galaxy)))

    print("#########################################################")

    # Here the actual fit starts. The best fit is plotted on the screen.
    # Gas emission lines are excluded from the pPXF fit using the GOODPIXELS keyword.
    #

    t = clock()



    if sky is not None:
        pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels=goodpixels, plot=plot, moments=moments, degree=degree, mdegree=mdegree, vsyst=dv, sky=sky)

    else:
        pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels=goodpixels, plot=plot, moments=moments, degree=degree, mdegree=mdegree, vsyst=dv)



    print("Formal errors:")
    print("     dV    dsigma   dh3      dh4")
    print("".join("%8.2g" % f for f in pp.error*np.sqrt(pp.chi2)))

    print('Elapsed time in PPXF: %.2f s' % (clock() - t))


    return (pp, sky, galaxy, noise)

#------------------------------------------------------------------------------

def load_CvD_templates(FWHM_galaxy, lamRange1, velscale, z=0.0, cvd_dir='/Data/stellarpops/CvD1.2'):


    """
    FWHM_galaxy: this is the FWHM of the galaxy spectrum in Angstoms.
    lamRange1: start and stop wavelength values of the galaxy spectrum
    z: the rough redshift
    """

    import glob
    from sam_python.CvD12tools import loadCvD12spec

    cvd = glob.glob('{}/t*.ssp'.format(cvd_dir))
    cvd.sort()

    #CvD Templates are at resolution 2000, so work out lamda/R for the middle wavelength in your array
    FWHM_tem = np.median(lamRange1*(1+z))/2000

    #FIXME
    #This depends on where in the wavelength range you are. For SWIFT:
    #9300A- 3.44A
    #10120A- 4.55A
    #8700- 4.04
    #To get a more accurate value for a certain lamda, fit a gaussian to a skyline!
    
    #Use Simon's CvDTools function to read in the CvD models and get them into proper units
    cvd_data=loadCvD12spec(cvd[0])
    #They're returned in Ryan's spectrum class. spec.lam is wavelengths, spec.flam is flux in lamda units
    lams=cvd_data.lam
    lamRange2=np.array([lams[0], lams[-1]])
    cdelt=lams[10]-lams[9]

    FWHM_dif = np.sqrt((FWHM_galaxy**2 - FWHM_tem**2).clip(0))
    sigma = FWHM_dif/2.355/cdelt # Sigma difference in pixels
    #Log Rebin one spectrum to get the length of the templates array right
    ssp=cvd_data.flam[0]
    


    sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp, velscale=velscale)

    t_num=63 #There are 63 spectra in the CvD models
    templates = np.empty((len(sspNew),t_num))

    #Do the same for all the models

    for j, filename in enumerate(cvd):
        cvd_data=loadCvD12spec(filename)
        print(filename)
        for ssp in cvd_data.flam:
            #print(ssp)
            ssp = ndimage.gaussian_filter1d(ssp,sigma)
            sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp, velscale=velscale)
            templates[:,j] = sspNew/np.median(sspNew) # Normalizes templates


    return templates, logLam2, lamRange2


def load_CaT_templates(FWHM_galaxy, velscale, cat_dir="FIXME"):

    """
    Load the CaT template library
    """

    
    CaT = glob.glob(cat_dir + 'Cun1.*.fits')
    CaT.sort()

    #CaT Templates are at a FWHM of 1.51 Angstroms
    FWHM_tem = 1.51

    #FIXME
    #This depends on where in the wavelength range you are. For SWIFT:
    #9300A- 3.44A
    #10120A- 4.55A
    #8700- 4.04
    #To get a more accurate value for a certain lamda, fit a gaussian to a skyline!
    
    hdu = fits.open(CaT[0])
    ssp = hdu[0].data
    h2 = hdu[0].header

    lam_temp = h2['CRVAL1'] + h2['CDELT1']*np.arange(h2['NAXIS1'])
    lamRange2 = [np.min(lam_temp), np.max(lam_temp)]
    print("LamRange2 is {}".format(lamRange2))
    sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp, velscale=velscale)
    templates = np.empty((sspNew.size, len(CaT)))
    t_num=len(CaT)
    FWHM_dif = np.sqrt((FWHM_galaxy**2 - FWHM_tem**2))
    sigma = FWHM_dif/2.355/h2['CDELT1'] # Sigma difference in pixels

    print("Sigma is {}".format(sigma))

    for j, fname in enumerate(CaT):
        hdu = fits.open(fname)
        ssp = hdu[0].data
        ssp = util.gaussian_filter1d(ssp, sigma)  # perform convolution with variable sigma
        sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp, velscale=velscale)

        templates[:, j] = sspNew/np.median(sspNew) # Normalizes templates

    return templates, loglam2


def mask_and_log_rebin(lamdas, spec, lamRange):

    """Mask and log rebin a spectrum

    """ 



    low_lam, high_lam=lamRange

    mask=np.where((lamdas>=low_lam) & (lamdas<=high_lam))[0]

    spec=spec[mask]

    spec, logLam1, velscale = util.log_rebin(lamRange, spec)

    return spec, logLam1, velscale







