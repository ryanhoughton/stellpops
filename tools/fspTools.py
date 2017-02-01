import numpy as np 
import scipy.constants as const
import ppxf_util as util
from astropy.io import fits
import scipy.interpolate as si


c = 299792.458
################################################################################################################################################################
def lnlike_CvD(theta, parameters):

    galaxy, noise, velscale, goodpixels, vsyst, interp_funct, logLams=parameters

    # galaxy=galaxy[goodpixels]
    # noise=noise[goodpixels]

    vel, sigma, age, Z, imf=theta

    fit_ranges=np.array([[4000, 4700], [4700, 5500], [8000, 8920], [9630, 10150]])*(np.exp(vel/c))

    # return_models=np.array((len(mask, len(fit_ranges))))
    # return_gals=np.array((len(mask, len(fit_ranges))))

    chisq=0

    for i, fit_range in enumerate(fit_ranges):

        mask=np.where((np.exp(logLams)>fit_range[0]) & (np.exp(logLams)<fit_range[1]))
        g=galaxy[mask]
        n=noise[mask]

        template=make_model_CvD(theta, interp_funct, logLams[mask])

        temp=convolve_template_with_losvd(template, vel, sigma, velscale=velscale, vsyst=vsyst)[:len(g)]




        morder=int(np.ceil(len(g)/10.0))
        poly=fit_legendre_polys(g/temp, morder)

        # return_models[:, i]=temp*poly
        # return_gals[:, i]=g

        import pdb; pdb.set_trace()

        chisq+=((g-temp*poly)/n)**2



    return -0.5*chisq#, return_models, return_gals


def lnprob_CvD(theta, parameters):

    lp=lnprior_CvD(theta)

    if np.isfinite(lp):

        ll=lnlike_CvD(theta, parameters)
        return lp + ll

    return -np.inf

def lnprior_CvD(theta):

    vel, sigma, age, Z, imf=theta

    if 1.0 < age < 13.5 and -1.50 < Z < 0.20 and 0.0 < imf < 3.5 and 0.0 < vel <10000.0 and 0.0 < sigma < 1000.0:
        return 0.0
    return -np.inf



def make_model_CvD(theta, interp_funct, logLams):

    _, _, age, Z, imf=theta

    model=interp_funct((logLams, age, Z, imf))

    return model

################################################################################################################################################################





################################################################################################################################################################

def lnprior(theta):

    vel, sigma, age, Z, imf, NaFe=theta

    if 10.0 < age < 17.782 and -0.40 < Z < 0.22 and 0.3 < imf < 3.3 and 0 < NaFe < 0.9 and 0.0 < vel <10000.0 and 0.0 < sigma < 1000.0:
        return 0.0
    return -np.inf

def lnprob(theta, parameters):

    lp=lnprior(theta)

    if np.isfinite(lp):

        ll, model, gal=lnlike(theta, parameters)
        return lp + ll, model

    return -np.inf, None

def make_model(theta, interp_funct, logLams):

    _, _, age, Z, imf, NaFe=theta

    model=interp_funct((logLams, age, Z, imf, NaFe))

    return model




def lnlike(theta, parameters):

    galaxy, noise, velscale, goodpixels, vsyst, interp_funct, logLams=parameters

    vel, sigma, age, Z, imf, NaFe=theta

    template=make_model(theta, interp_funct, logLams)

    temp=convolve_template_with_losvd(template, vel, sigma, velscale=velscale, vsyst=vsyst)[:len(galaxy)]

    poly=fit_legendre_polys(galaxy/temp, 100)   

    if goodpixels is None:
        goodpixels=np.arange(len(galaxy))


    g=galaxy[goodpixels]
    t=temp[goodpixels]
    n=noise[goodpixels]
    p=poly[goodpixels]

    chisq=((g-t*p)/n)**2

    return -0.5*np.sum(chisq), temp*poly, galaxy

################################################################################################################################################################




################################################################################################################################################################
def fit_legendre_polys(ratio, morder):

    x_vals=np.linspace(-1, 1, len(ratio))
    coeffs=np.polynomial.legendre.legfit(x_vals, ratio, morder)

    polynomial=np.polynomial.legendre.legval(x_vals, coeffs)

    return polynomial


def convolve_template_with_losvd(template, vel=0.0, sigma=0.0, velscale=None, vsyst=0.0):

    t_rfft, npad=_templates_rfft(template)
    losvd_rfft=_losvd_rfft(vel, sigma, npad, velscale, npad, vsyst=vsyst)

    convolved_t=np.fft.irfft(t_rfft*losvd_rfft)

    return convolved_t


def _losvd_rfft(vel, sig, pad, velscale, npad, vsyst=0.0):


    nl = npad//2 + 1

    vel=(vel+vsyst)/velscale
    sig/=velscale

    a, b = [vel, 0.0]/sig

    #print 'Vel is {}, Vsyst is {}, sig/vescale is {}, a is {}'.format(vel, vsyst, sig, a)



    w = np.linspace(0, np.pi*sig, nl)
    #analytic FFT of LOSVD
    losvd_rfft = np.exp(1j*a*w - 0.5*(1 + b**2)*w**2)

    return np.conj(losvd_rfft)

    

def _templates_rfft(templates):
    """
    Pre-compute the FFT (of real input) of all templates

    """
    npad = 2**int(np.ceil(np.log2(templates.shape[0])))
    templates_rfft = np.fft.rfft(templates, n=npad, axis=0)

    return templates_rfft, npad


def gaussian(vel, sigma, velscale, npad, h3h4=None):
    """
    FIXME
    """
    import matplotlib.pyplot as plt 
    sigmaKernel=sigma
    nsig=5.0

    dv = np.ceil(nsig*sigmaKernel/velscale) 
    nv = 2*dv + 1
    v = np.linspace(dv,-dv,nv) 
    w = (v - vel/velscale) / (sigmaKernel/velscale)
    w2= w*w
    if h3h4 != None:
        h3=h3h4[0]
        h4=h3h4[1]
        poly = 1.0 + h3/np.sqrt(3.0)*(w*(2.0*w2-3.0)) + \
               h4/np.sqrt(24.0)*(w2*(4.0*w2-12.0)+3.0)
    else:
        poly = np.ones(nv)

    losvd = np.exp(-0.5*w2)/np.sum(np.exp(-0.5*w2)) * poly 
    
    #import pdb; pdb.set_trace()
    #npad = int(2**np.ceil(np.log2(nf)))
    
    nv=int(nv)
    losvd_final=np.zeros(npad/2)
    #import pdb; pdb.set_trace()
    losvd_final[:nv]=losvd

    losvd_final = np.roll(losvd_final, (2 - nv)//2, axis=0)

    return losvd_final

def plot_chain(samples, labels=None):



    a, b, c=samples.shape

    

    fig, ax=plt.subplots(nrows=c, ncols=1)

    for i in range(a):
        for j in range(c):
            ax.flatten()[j].plot(samples[i, :, j])

    if labels is not None:
        assert len(labels)==c, 'Must have same number of labels as dimensions'
        for j in range(c):
            ax.flatten()[j].set_ylabel(labels[j])


    return ax



def interpolate_Miles_ML(IMF_type='bi', dir_name='/Data/stellarpops/Miles/Mass_files'):

    #Return an interpolate object which has M/L in the V band as a function of age, Z/H and IMF

    if IMF_type=='bi':
        table_name='ssp_phot_Padova00_BI_v10.txt'
    elif IMF_type=='uni' or IMF_type =='un':
        table_name='ssp_phot_Padova00_UN_v10.txt'


    #Use np.genfromtxt as astropy.io.ascii doesn't seem to work...
    table=np.genfromtxt('{}/{}'.format(dir_name, table_name))

    ML_V=table[1:, -13].reshape(50, 7, 12)
    _ages=np.unique(table[1:, 3])
    _Z_Hs=np.unique(table[1:, 2])
    _imfs=np.unique(table[1:, 1])

    interp=si.RegularGridInterpolator((_ages, _Z_Hs, _imfs), ML_V)

    return interp

################################################################################################################################################################


################################################################################################################################################################
#NGC1277 Tools
################################################################################################################################################################

################################################################################################################################################################
def prepare_CvD2_templates(templates_lam_range, velscale, verbose=True):
    import glob

    vcj_models=sorted(glob.glob('/Data/stellarpops/CvD2/vcj_models/VCJ_*.s100'))
    temp_lamdas, x35, x3, x23, kroupa, flat=np.genfromtxt(vcj_models[0], unpack=True)

    n_ages=7
    n_zs=5
    n_imfs=5

    Zs=['m1.5', 'm1.0', 'm0.5', 'p0.0', 'p0.2']
    ages=[1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.5]
    model_imfs_order=['x35', 'x3', 'x23', 'kroupa', 'flat']

    t_mask = ((temp_lamdas > templates_lam_range[0]) & (temp_lamdas <templates_lam_range[1]))


    sspNew, logLam_template, template_velscale = util.log_rebin(templates_lam_range, x35[t_mask], velscale=velscale)
    templates=np.empty((len(sspNew), n_ages, n_zs, n_imfs))

    


    for a, Z in enumerate(Zs):    
        for b, age in enumerate(ages):
            model=glob.glob('/Data/stellarpops/CvD2/vcj_models/VCJ_*{}*{}.ssp.s100'.format(Z, age))[0]
            print 'Loading {}'.format(model)
            data=np.genfromtxt(model)

            for c, counter in enumerate(reversed(range(1, data.shape[-1]))):
                
                sspNew, logLam_template, template_velscale = util.log_rebin(templates_lam_range, data[:, counter][t_mask], velscale=velscale)

                templates[:, b, a, c]=sspNew/np.median(sspNew)

    return templates, logLam_template
################################################################################################################################################################  

################################################################################################################################################################
def prepare_CvD_interpolator(templates_lam_range, velscale, verbose=True):

    templates, logLam_template=prepare_CvD2_templates(templates_lam_range, velscale, verbose=verbose)

    ages=[1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.5]
    Zs=[-1.5, -1.0, -0.5, 0.0, 0.2]
    imfs=[0.0, 1.0, 2.3, 3.0, 3.5]


    linear_interp=si.RegularGridInterpolator(((logLam_template, ages, Zs, imfs)), templates)

    return linear_interp, logLam_template

################################################################################################################################################################

################################################################################################################################################################
def prepare_Miles_templates(templates_lam_range, velscale, template_dictionary=None, NaFe=0.0, imf_type='uni', verbose=True):


    #Prepare the Miles templates for use with pPXF
    #Loading all the Miles templates takes ages, so if you already have them (in an ipython session, say), pass them to this function to save time

    if NaFe!=0.0:
        assert imf_type=='bi', "Can't have NaFe>0.0 if imf type is 'uni'"

    if template_dictionary is None:
        #Load the M16 templates
        from stellarpops.tools import Miles16tools as M16
        template_dictionary=M16.load_eMILES_spectra(basedir='/Data/stellarpops/Miles', NaFe=NaFe, verbose=verbose, imf_type=imf_type)

    #Get the imfs, ages and metallicities of the templates
    imfs=sorted(template_dictionary.keys())
    Zs=template_dictionary[imfs[0]].Z[:, 0]
    ages=template_dictionary[imfs[0]].age[0, :]

    temp_lamdas=template_dictionary[imfs[0]].lam

    t_mask = ((temp_lamdas > templates_lam_range[0]) & (temp_lamdas <templates_lam_range[1]))
    sspNew, logLam_template, template_velscale = util.log_rebin(templates_lam_range, template_dictionary[imfs[0]].flam[0, 0, :][t_mask], velscale=velscale)
    templates=np.empty((len(imfs), template_dictionary[imfs[0]].flam.shape[0], template_dictionary[imfs[0]].flam.shape[1], len(sspNew)))

    #Make the template array
    for i, imf in enumerate(imfs):
        for j in range(len(Zs)):
            for k in range(len(ages)):

                sspNew, logLam_template, _ = util.log_rebin(templates_lam_range, template_dictionary[imf].flam[j, k, :][t_mask], velscale=velscale)
                sspNew/=np.median(sspNew)
                templates[i, j, k, :]=sspNew

        #templates[i, :]=template_dictionary[imfs].flam/(np.repeat(np.median(temp[imfs].flam, axis=2), temp[imfs].flam.shape[-1]).reshape(7, 12, -1))

    #pPXF expects wavelength axis first
    templates=templates.T

    #temp_wave=np.logspace(np.log10(templates_lam_range[0]), np.log10(templates_lam_range[1]), templates.shape[0])

    return templates, logLam_template, template_dictionary

################################################################################################################################################################


def NGC1277_CVD_read_in_data_CvD(file = 'Data/n1277b_cen.dat', c=299792.458):
    ##############################################################################################################################################################

    # Read in the central spectrum from CvD
    #
    lamdas, flux, errors, CvD_weights, inst_res=np.genfromtxt(file, unpack=True)
    #Redshift of NGC1277, used to get the initial velocity
    z=0.017044

    #The CvD data has a chip gap between 5625 and 7071 A. This screws up the fitting if you forget about it (obviously!).
    #This 'gap' array is insterted into the normalised flux and error arrays, before we mask it out in the fitting.

    flux_median=np.median(flux)

    flux/=flux_median
    errors/=flux_median

    chip_gap_start=lamdas[2424]
    chip_gap_stop=lamdas[2425]

    # gap=np.ones(int(chip_gap_stop-chip_gap_start)) 
    # flux_edited=np.insert(flux/flux_median, 2424, gap)
    # errors_edited=np.insert(errors/flux_median, 2424, gap)

    # quick_and_dirty_lamdas=np.linspace(lamdas[0], lamdas[-1], len(flux_edited))



    #Mask to only fit some wavelength ranges?

    # lower=lamdas.min()
    # upper=lamdas.max()
    lower=lamdas.min()
    upper=lamdas.max()

    assert (lower>=lamdas.min()) & (upper<=lamdas.max()), 'Lower and upper limits must be within the wavelength ranges of the data'

    lam_range_gal=np.array([lower, upper])
    mask=np.where((lamdas>lower) & (lamdas<upper))

    flux=flux[mask]
    errors=errors[mask]

    #Log rebin them
    galaxy, logLam, velscale = util.log_rebin(lam_range_gal, flux)
    noise, _, _=util.log_rebin(lam_range_gal, errors)   


    #GoodPixels from pPXF
    goodpixels = np.arange(len(galaxy)) 
    

    return galaxy, noise, velscale, goodpixels, lam_range_gal, logLam

    ################################################################################################################################################################




def NGC1277_SWIFT_read_in_data(file='Data/SPV_NGC1277.dat'):

    lamdas, flux, variance, inst_res=np.genfromtxt(file, unpack=True)

    flux_median=np.median(flux)
    errors=np.sqrt(variance)

    flux/=flux_median
    errors/=flux_median

    lamdas=lamdas*10**4

    lower=lamdas.min()
    upper=lamdas.max()

    lam_range_gal=np.array([lower, upper])

    galaxy, logLam, velscale = util.log_rebin(lam_range_gal, flux)
    noise, _, _=util.log_rebin(lam_range_gal, errors)  

    goodpixels = np.arange(len(galaxy))


    return galaxy, noise, velscale, goodpixels, lam_range_gal, logLam


def NGC1277_CVD_read_in_data_MILES(file = 'Data/n1277b_cen.dat', c=299792.458):
    ##############################################################################################################################################################

    # Read in the central spectrum from CvD
    #
    lamdas, flux, errors, CvD_weights, inst_res=np.genfromtxt(file, unpack=True)
    #Redshift of NGC1277, used to get the initial velocity
    z=0.017044

    #The CvD data has a chip gap between 5625 and 7071 A. This screws up the fitting if you forget about it (obviously!).
    #This 'gap' array is insterted into the normalised flux and error arrays, before we mask it out in the fitting.

    flux_median=np.median(flux)

    chip_gap_start=lamdas[2424]
    chip_gap_stop=lamdas[2425]

    gap=np.ones(int(chip_gap_stop-chip_gap_start)) 
    flux_edited=np.insert(flux/flux_median, 2424, gap)
    errors_edited=np.insert(errors/flux_median, 2424, gap)

    quick_and_dirty_lamdas=np.linspace(lamdas[0], lamdas[-1], len(flux_edited))



    #Mask to only fit some wavelength ranges?

    # lower=lamdas.min()
    # upper=lamdas.max()
    lower=lamdas.min()
    upper=lamdas.max()

    assert (lower>=lamdas.min()) & (upper<=lamdas.max()), 'Lower and upper limits must be within the wavelength ranges of the data'

    lam_range_gal=np.array([lower, upper])
    mask=np.where((quick_and_dirty_lamdas>lower) & (quick_and_dirty_lamdas<upper))

    flux_edited=flux_edited[mask]
    errors_edited=errors_edited[mask]

    #Log rebin them
    galaxy, logLam, velscale = util.log_rebin(lam_range_gal, flux_edited)
    noise, _, _=util.log_rebin(lam_range_gal, errors_edited)   


    #GoodPixels from pPXF
    goodpixels = np.arange(len(galaxy))

    #Points corresponding to the chipgap- remove from goodpixels
    chip_gap_lamdas=np.where((logLam>np.log(chip_gap_start)) & (logLam<np.log(chip_gap_stop)))[0]
    goodpixels = np.array([x for x in goodpixels if x not in chip_gap_lamdas])   
    

    return galaxy, noise, velscale, goodpixels, lam_range_gal, logLam

    ################################################################################################################################################################



def read_in_templates(lam_range_gal, velscale, pad=500.0, imf_type='bi', verbose=True):
    ################################################################################################################################################################
    #Wavelength Range of the templates is much bigger than the galaxy, so we mask the templates to just a bit bigger than the galaxy wavelength range

    lam_range_temp = [lam_range_gal[0]-pad, lam_range_gal[1]+pad]

    #get the NaFe=0.0 templates and clip the very low metallicities  
    templates, logLam_template, template_dictionary=prepare_Miles_templates(lam_range_temp, velscale, imf_type=imf_type, verbose=verbose)
    #templates, logLam_template, template_dictionary=prepare_Miles_templates(lam_range_temp, velscale, template_dictionary=template_dictionary, imf_type='uni')
    imfs=sorted(template_dictionary.keys())
    Zs=template_dictionary[imfs[0]].Z[:, 0]
    ages=template_dictionary[imfs[0]].age[0, :]

    templates=templates[:, -6:, -3:, :]
    Zs=Zs[-3:]

    ages=ages[-6:]

    NaFe3templates, logLam_template, NaFe3template_dictionary=prepare_Miles_templates(lam_range_temp, velscale, NaFe=0.3, imf_type=imf_type, verbose=verbose)
    NaFe3templates=NaFe3templates[:, -6:, -3:, :]

    NaFe6templates, logLam_template, NaFe6template_dictionary=prepare_Miles_templates(lam_range_temp, velscale, NaFe=0.6, imf_type=imf_type, verbose=verbose)
    NaFe6templates=NaFe6templates[:, -6:, -3:, :]

    NaFe9templates, logLam_template, NaFe9template_dictionary=prepare_Miles_templates(lam_range_temp, velscale, NaFe=0.9, imf_type=imf_type, verbose=verbose)
    NaFe9templates=NaFe9templates[:, -6:, -3:, :]
    ################################################################################################################################################################

    all_templates=np.stack((templates, NaFe3templates, NaFe6templates, NaFe9templates), axis=4)
    NaFes=np.array([0.0, 0.3, 0.6, 0.9])
    Zs=np.array([-0.4, 0.0, 0.22])
    imfs=[float(imf.strip('bi')) for imf in imfs]


    linear_interp=si.RegularGridInterpolator(((logLam_template, ages, Zs, imfs, NaFes)), all_templates)

    return linear_interp, logLam_template, lam_range_temp

    ################################################################################################################################################################



def NGC1277_CvD_set_up_emcee_parameters_MILES(file = 'n1277b_cen.dat', verbose=True):

    
     # speed of light in km/s

    galaxy, noise, velscale, goodpixels, lam_range_gal, logLam=NGC1277_CVD_read_in_data_MILES(file=file, c=c)

    linear_interp, logLam_template, lam_range_temp=read_in_templates(lam_range_gal, velscale, pad=500.0, verbose=verbose)

    dv = c*np.log(lam_range_temp[0]/lam_range_gal[0])  # km/s


    return [galaxy, noise, velscale, goodpixels, dv, linear_interp, logLam_template], logLam


def NGC1277_CvD_set_up_emcee_parameters_CvD(file = 'n1277b_cen.dat', verbose=True):

    
    

    galaxy, noise, velscale, goodpixels, lam_range_gal, logLam=NGC1277_CVD_read_in_data_CvD(file=file, c=c)


    pad=500.0
    lam_range_temp = [lam_range_gal[0]-pad, lam_range_gal[1]+pad]
    linear_interp, logLam_template, =prepare_CvD_interpolator(lam_range_temp, velscale, verbose=True)


    dv = c*np.log(lam_range_temp[0]/lam_range_gal[0])  # km/s


    return [galaxy, noise, velscale, goodpixels, dv, linear_interp, logLam_template], logLam


################################################################################