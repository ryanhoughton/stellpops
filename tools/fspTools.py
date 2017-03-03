import numpy as np 
import scipy.constants as const
from ppxf import ppxf_util as util
from astropy.io import fits
import scipy.interpolate as si
from stellarpops.tools import CD12tools as CT
import scipy.constants as const

#Speed of light in km/s
c_light = const.c/1000.0
################################################################################################################################################################
def lnlike_CvD(theta, parameters, plot=False):

    galaxy, noise, velscale, goodpixels, vsyst, interp_funct, correction_interps, logLams, logLam_gal, fit_wavelengths=parameters

    general_interp, na_interp, positive_only_interp=correction_interps

    vel, sigma=theta[0], theta[1]
    Na_abundance=theta[2]
    general_abundances=theta[3:-3]
    positive_abundances=theta[-3]
    Z, imf=theta[-2:]
    
    #Don't fit age- keep it fixed at 13.5 Gyr
    age=13.5

    base_template=make_model_CvD(theta, interp_funct, logLams)


    #If positive abundances has only one element, run it as a usual interpoaltor without going through the get correction function
    #Lists are iterables wheres single floats are not, so this check will pass for floats but not arrays or lists
    if not hasattr(positive_abundances,'__iter__'):
        positive_only_correction=positive_only_interp((positive_abundances, age, Z, logLams))
    else:
        positive_only_correction=get_correction(positive_only_interp, logLams, np.arange(len(positive_abundances)), positive_abundances, age, Z)



    general_correction=get_correction(general_interp, logLams, np.arange(len(general_abundances)), general_abundances, age, Z)   

    na_correction=na_interp((Na_abundance, age, Z, logLams))

    #old_template=template*general_correction*positive_only_correction*na_correction
    template=np.exp(np.log(base_template)+general_correction+positive_only_correction+na_correction)





    temp=convolve_template_with_losvd(template, vel, sigma, velscale=velscale, vsyst=vsyst)[:len(galaxy)]

    chisq=0

    fit_ranges=fit_wavelengths*(np.exp(vel/c_light))

    #Do all the plotting, if required
    if plot==True:
        import matplotlib.pyplot as plt 
        import matplotlib.gridspec as gridspec
        import matplotlib.ticker as ticker   

        gs_1 = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[2, 1])
        gs_2 = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[2, 1])
        gs_3 = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[2, 1])
        gs_4 = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[2, 1])
        

        fig=plt.figure(figsize=(14, 10))
        axs=np.empty((4, 2), dtype='object')
        outer_grid=gridspec.GridSpec(2, 2)

        for i in range(4):
            inner_grid = gridspec.GridSpecFromSubplotSpec(2, 2, width_ratios=[1, 1], height_ratios=[2, 1], subplot_spec=outer_grid[i//2, i%2], hspace=0.0)
            axs[i, 0] = fig.add_subplot(inner_grid[0, :2])
            axs[i, 1] = fig.add_subplot(inner_grid[1, :2], sharex=axs[i, 0])
            plt.setp(axs[i, 0].get_xticklabels(), visible=False)



    for i, fit_range in enumerate(fit_ranges):



        #tmask=np.where((np.exp(logLams)>fit_range[0]) & (np.exp(logLams)<fit_range[1]))
        gmask=np.where((np.exp(logLam_gal)>fit_range[0]) & (np.exp(logLam_gal)<fit_range[1]))
    
        g=galaxy[gmask]
        n=noise[gmask]
        t=temp[gmask]

        morder=int((fit_range[1]-fit_range[0])/100)
        poly=fit_legendre_polys(g/t, morder)


        chisq+=np.sum((((g-t*poly)/n)**2))




        if plot:
            x=np.exp(logLam_gal[gmask])/(np.exp(vel/c_light))
            axs[i, 0].plot(x, g, c='k', linewidth=1.5)
            axs[i, 0].plot(x, poly*t, c='b', linewidth=2.0)
            axs[i, 0].fill_between(x, g-n, g+n, facecolor='k', alpha=0.3)

            axs[i, 1].plot(x, 100*(g-poly*t)/(poly*t), c='k', linewidth=1.5)
            axs[i, 1].axhline(0.0, linestyle='dashed', c='k')

            
            axs[i, 0].set_xlim([x.min(), x.max()])
            axs[i, 1].set_ylim([-5, 5])

            axs[i, 1].set_xlabel('Rest Wavelength (A)')
            axs[i, 0].set_ylabel('Flux (Arbitrary Units)')
            axs[i, 1].set_ylabel('Residuals (%)')

            #Avoid the overlapping labels
            axs[i, 0].yaxis.set_major_locator(ticker.MaxNLocator(prune='lower'))
            axs[i, 1].yaxis.set_major_locator(ticker.MultipleLocator(2))

      
            
            




    if plot:
        return -0.5*chisq, (fig, axs)

    return -0.5*chisq#, return_models, return_gals


def lnprob_CvD(theta, parameters):

    lp=lnprior_CvD(theta)

    if np.isfinite(lp):

        ll=lnlike_CvD(theta, parameters)
        return lp + ll

    return -np.inf

def lnprior_CvD(theta):


    vel, sigma=theta[0], theta[1]
    Na_abundance=theta[2]
    general_abundances=theta[3:-3]
    positive_abundances=theta[-3]
    Z, imf=theta[-2:]
    
    #Don't fit age- keep it fixed at 13.5 Gyr
    age=13.5

    if 0.0 < vel < 7000.0 and 0.0 < sigma < 500.0:

        if np.all(general_abundances>=-0.45) and np.all(general_abundances<=0.45)  and np.all(positive_abundances>=0.0) and np.all(positive_abundances<=0.45) and -0.45 <= Na_abundance <= 1.0 and 1.0 < age <= 14.0 and -0.5 < Z < 0.4 and 0.5 < imf <3.5:
            return 0.0

    return -np.inf




def make_model_CvD(theta, interp_funct, logLams):

    vel, sigma=theta[0], theta[1]
    Na_abundance=theta[2]
    general_abundances=theta[3:-3]
    positive_abundances=theta[-3]
    Z, imf=theta[-2:]
    
    #Don't fit age- keep it fixed at 13.5 Gyr
    age=13.5

    model=interp_funct((logLams, age, Z, imf))

    return model


def get_correction(interpolator, logLams, elems, abunds, age, Z):

    #The interpolator expects a list of 6 numbers. Meshgrid the two arrays which are of different lengths
    # (the indices and the number of elements to enhance) and then create lists of ages, Zs and IMFs of the correct
    # shapes. Then do the same for the abundances. Stack together and pass to the interpolator object!

    points = np.meshgrid(elems, logLams, indexing='ij')
    flat = np.array([m.flatten() for m in points])
    #flat is now an array of points of shape 2, len(indices)*len(elems)
    #make arrays of the other variables of length len(indices)*len(elems)
    ages=np.ones_like(points[0])*age
    Zs=np.ones_like(points[0])*Z
    

    #Get the correct abundance for each element- luckily we can index the abundance array by the integer element array
    abunds=abunds[points[0]]

    # import pdb; pdb.set_trace()

    #Stack together
    xi=np.vstack((flat[0, :], abunds.ravel(), ages.ravel(), Zs.ravel(), flat[1, :]))
    #Do the interpolation
    out_array = interpolator(xi.T)
    #reshape everything to be (len(indices), len(elements))
    result = out_array.reshape(*points[0].shape)

    return np.sum(result, axis=0)



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
    import os
    template_glob=os.path.expanduser('~/z/Data/stellarpops/CvD2/vcj_models/VCJ_*.s100')

    vcj_models=sorted(glob.glob(template_glob))
    temp_lamdas, x35, x3, x23, kroupa, flat=np.genfromtxt(vcj_models[-1], unpack=True)

    n_ages=7
    n_zs=5
    n_imfs=5

    


    Zs=['m1.5', 'm1.0', 'm0.5', 'p0.0', 'p0.2']
    ages=[1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.5]
    model_imfs_order=['x35', 'x3', 'x23', 'kroupa', 'flat']

    t_mask = ((temp_lamdas > templates_lam_range[0]) & (temp_lamdas <templates_lam_range[1]))



    y=x35[t_mask]
    x=temp_lamdas[t_mask]
    #Make a new lamda array, carrying on the delta lamdas of high resolution bit
    new_x=temp_lamdas[t_mask][0]+0.9*(np.arange(np.ceil((temp_lamdas[t_mask][-1]-temp_lamdas[t_mask][0])/0.9))+1)
    interp=si.interp1d(x, y, fill_value='extrapolate')
    out=interp(new_x)

    sspNew, logLam_template, template_velscale = util.log_rebin(templates_lam_range, out, velscale=velscale)
    templates=np.empty((len(sspNew), n_ages, n_zs, n_imfs))



    for a, Z in enumerate(Zs):    
        for b, age in enumerate(ages):
            model=glob.glob(os.path.expanduser('~/z/Data/stellarpops/CvD2/vcj_models/VCJ_*{}*{}.ssp.s100'.format(Z, age)))[0]
            print 'Loading {}'.format(model)
            data=np.genfromtxt(model)

            for c, counter in enumerate(reversed(range(1, data.shape[-1]))):
                
                #Interpolate templates onto a uniform wavelength grid and then log-rebin
                y=data[:, counter][t_mask]   
                x=temp_lamdas[t_mask]
                #Make a new lamda array, carrying on the delta lamdas of high resolution bit
                new_x=temp_lamdas[t_mask][0]+0.9*(np.arange(np.ceil((temp_lamdas[t_mask][-1]-temp_lamdas[t_mask][0])/0.9))+1)

                interp=si.interp1d(x, y, fill_value='extrapolate')
                out=interp(new_x)
                
                sspNew, logLam_template, template_velscale = util.log_rebin(templates_lam_range, out, velscale=velscale)

            

                templates[:, b, a, c]=sspNew/np.median(sspNew)

    
  

    



    return templates, logLam_template


def prepare_CvD2_element_templates(templates_lam_range, velscale, elements, verbose=True):

    import glob

    import os
    template_glob=os.path.expanduser('~/z//Data/stellarpops/CvD2/vcj_models/VCJ_*.s100')

    var_elem_spectra=CT.load_varelem_CvD16ssps(dirname=os.path.expanduser('~/z/Data/stellarpops/CvD2'), folder='atlas_rfn_v3', imf='kroupa')

    ages=var_elem_spectra['Solar'].age[:, 0]
    Zs=var_elem_spectra['Solar'].Z[0, :]
    n_ages=len(ages)
    n_Zs=len(Zs)

    temp_lamdas=var_elem_spectra['Solar'].lam

    t_mask = ((temp_lamdas > templates_lam_range[0]) & (temp_lamdas <templates_lam_range[1]))


    positive_only_elems, Na_elem, normal_elems=elements

    elem_steps=[-0.45, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.45]
    Na_elem_steps=[-0.45, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    positive_only_elem_steps=[0.0, 0.1, 0.2, 0.3, 0.45]


    x=var_elem_spectra['Solar'].lam[t_mask]
    y=var_elem_spectra['Solar'].flam[-1, -1, t_mask]
    #Make a new lamda array, carrying on the delta lamdas of high resolution bit
    new_x=var_elem_spectra['Solar'].lam[t_mask][0]+0.9*(np.arange(np.ceil((var_elem_spectra['Solar'].lam[t_mask][-1]-var_elem_spectra['Solar'].lam[t_mask][0])/0.9))+1)
    interp=si.interp1d(x, y, fill_value='extrapolate')
    data=interp(new_x)


    sspNew, logLam_template, template_velscale = util.log_rebin(templates_lam_range, data, velscale=velscale)

    print 'SSP New is {}'.format(len(sspNew))

    positive_only_templates=np.empty((len(positive_only_elems), len(positive_only_elem_steps), n_ages, n_Zs, len(sspNew)))
    general_templates=np.empty((len(normal_elems), len(elem_steps), n_ages, n_Zs, len(sspNew)))
    
    na_templates=np.empty((len(Na_elem_steps), n_ages, n_Zs, len(sspNew)))

    print 'Making the Positive-Only Correction templates'
    #Do the positve only correction templates:
    for a, elem in enumerate(positive_only_elems):
        print '\t{}'.format(elem)
        for b, step in enumerate(positive_only_elem_steps):
            for c, _ in enumerate(ages):
                for d, _ in enumerate(Zs):

                    if step !=0.0:
                        y=(var_elem_spectra[elem].flam[b, c, t_mask]/var_elem_spectra['Solar'].flam[b, c, t_mask] - 1.0)*((10**(step)-1.0)/(10**(0.3)-1.0))


                    else:
                        y=np.zeros_like(var_elem_spectra['Solar'].flam[c, d, t_mask])

                    x=var_elem_spectra[elem].lam[t_mask]
                    #Make a new lamda array, carrying on the delta lamdas of high resolution bit
                    new_x=var_elem_spectra[elem].lam[t_mask][0]+0.9*(np.arange(np.ceil((var_elem_spectra[elem].lam[t_mask][-1]-var_elem_spectra[elem].lam[t_mask][0])/0.9))+1)
                    interp=si.interp1d(x, y, fill_value='extrapolate')
                    data=interp(new_x)
                            
                    sspNew, logLam_template, template_velscale = util.log_rebin(templates_lam_range, data, velscale=velscale)

                    positive_only_templates[a, b, c, d, :]=sspNew#/np.median(sspNew)


    # import matplotlib.pyplot as plt
    # plt.figure()
    print 'Making the General Correction templates'
    #Do the general templates
    for a, elem in enumerate(normal_elems):
        
        print '\t{}'.format(elem)
        for b, step in enumerate(elem_steps):
            for c, _ in enumerate(ages):
                for d, _ in enumerate(Zs):

                    if step>0.0:
                        e='{}+'.format(elem)
                        gen_step=step
                    elif step<0.0:
                        e='{}-'.format(elem)
                        gen_step=np.abs(step)

                    if step !=0.0:
                        y=(var_elem_spectra[e].flam[c, d, t_mask]/var_elem_spectra['Solar'].flam[c, d, t_mask]-1)*((10**(gen_step)-1.0)/(10**(0.3)-1.0))
                    else:
                        y=np.zeros_like(var_elem_spectra['Solar'].flam[c, d, t_mask])

                    x=var_elem_spectra[e].lam[t_mask]
                    #Make a new lamda array, carrying on the delta lamdas of high resolution bit
                    new_x=var_elem_spectra[e].lam[t_mask][0]+0.9*(np.arange(np.ceil((var_elem_spectra[e].lam[t_mask][-1]-var_elem_spectra[e].lam[t_mask][0])/0.9))+1)
                    interp=si.interp1d(x, y, fill_value='extrapolate')
                    data=interp(new_x)
                            
                    sspNew, logLam_template, template_velscale = util.log_rebin(templates_lam_range, data, velscale=velscale)

                    general_templates[a, b, c, d, :]=sspNew#/np.median(sspNew)

        #     plt.plot(sspNew, label='{} {}'.format(e, step))

        # plt.legend()
        # import pdb; pdb.set_trace()


    #Do the Na templates:
    print 'Making the Na Correction template'
    for a, step in enumerate(Na_elem_steps):
        for b, _ in enumerate(ages):
            for c, _ in enumerate(Zs):

                if step <0.0:
                    e='Na-'
                    base_enhancement=0.3
                    Na_step=np.abs(step)                    
                elif 0.0<=step<0.45:
                    e='Na+'
                    base_enhancement=0.3
                    Na_step=step
                elif 0.45<=step<0.75:
                    e='Na+0.6'
                    base_enhancement=0.6
                    Na_step=step
                elif 0.75<=step<1.0:
                    e='Na+0.9'
                    base_enhancement=0.9
                    Na_step=step
                
                if step !=0.0:
                    y=(var_elem_spectra[e].flam[b, c, t_mask]/var_elem_spectra['Solar'].flam[b, c, t_mask]-1)*((10**(Na_step)-1.0)/(10**(base_enhancement)-1.0))

                else:

                    y=np.zeros_like(var_elem_spectra['Solar'].flam[c, d, t_mask])


                x=var_elem_spectra[e].lam[t_mask]
                #Make a new lamda array, carrying on the delta lamdas of high resolution bit
                new_x=var_elem_spectra[e].lam[t_mask][0]+0.9*(np.arange(np.ceil((var_elem_spectra[e].lam[t_mask][-1]-var_elem_spectra[e].lam[t_mask][0])/0.9))+1)
                interp=si.interp1d(x, y, fill_value='extrapolate')
                data=interp(new_x)
                    
                sspNew, logLam_template, template_velscale = util.log_rebin(templates_lam_range, data, velscale=velscale)
                na_templates[a, b, c, :]=sspNew



    return [general_templates, na_templates, positive_only_templates], logLam_template



################################################################################################################################################################  

################################################################################################################################################################
def prepare_CvD_interpolator(templates_lam_range, velscale, verbose=True):

    templates, logLam_template=prepare_CvD2_templates(templates_lam_range, velscale, verbose=verbose)

    ages=[  1.,   3.,   5.,  7., 9.,  11.0, 13.5]
    Zs=[-1.5, -1.0, -0.5, 0.0, 0.2]
    imfs=[0.0, 1.8, 2.3, 3.0, 3.5]


    linear_interp=si.RegularGridInterpolator(((logLam_template, ages, Zs, imfs)), templates, bounds_error=False, fill_value=None)

    return linear_interp, logLam_template

################################################################################################################################################################

################################################################################################################################################################
def prepare_CvD_correction_interpolators(templates_lam_range, velscale, elements, verbose=True):

    all_corrections, logLam_template=prepare_CvD2_element_templates(templates_lam_range, velscale, elements, verbose=verbose)

    general_templates, na_templates, positive_only_templates=all_corrections

    positive_only_elems, Na_elem, normal_elems=elements



    ages=np.array([  1.,   3.,   5.,   9.,  13.])
    Zs=[-1.5, -1.0, -0.5, 0.0, 0.2]


    elem_steps=[-0.45, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.45]
    Na_elem_steps=[-0.45, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    positive_only_elem_steps=[0.0, 0.1, 0.2, 0.3, 0.45]

    # np.save('na_templates.npy', na_templates)
    # np.save('general_templates.npy', general_templates)

    general_interp=si.RegularGridInterpolator(((np.arange(len(normal_elems)), elem_steps, ages, Zs, logLam_template)), general_templates, bounds_error=False, fill_value=None, method='linear')
    na_interp=si.RegularGridInterpolator(((Na_elem_steps, ages, Zs, logLam_template)), na_templates, bounds_error=False, fill_value=None, method='linear')

    #If we only have one positive element to check, we need to do something different- can't have a dimension with only one element in RegularGridInterpolator apparently.
    if len(positive_only_elems)>1:
        positive_only_interp=si.RegularGridInterpolator(((np.arange(len(positive_only_elems)), positive_only_elem_steps, ages, Zs, logLam_template)), positive_only_templates, bounds_error=False, fill_value=None, method='linear')
    else:
        positive_only_interp=si.RegularGridInterpolator((positive_only_elem_steps, ages, Zs, logLam_template), positive_only_templates[0, :], bounds_error=False, fill_value=None, method='linear')


    #import pdb; pdb.set_trace()

    correction_interps=[general_interp, na_interp, positive_only_interp]

    return correction_interps, logLam_template

################################################################################################################################################################



# ################################################################################################################################################################
# def prepare_Miles_templates(templates_lam_range, velscale, template_dictionary=None, NaFe=0.0, imf_type='uni', verbose=True):


#     #Prepare the Miles templates for use with pPXF
#     #Loading all the Miles templates takes ages, so if you already have them (in an ipython session, say), pass them to this function to save time

#     if NaFe!=0.0:
#         assert imf_type=='bi', "Can't have NaFe>0.0 if imf type is 'uni'"

#     if template_dictionary is None:
#         #Load the M16 templates
#         from stellarpops.tools import Miles16tools as M16
#         template_dictionary=M16.load_eMILES_spectra(basedir='/Data/stellarpops/Miles', NaFe=NaFe, verbose=verbose, imf_type=imf_type)

#     #Get the imfs, ages and metallicities of the templates
#     imfs=sorted(template_dictionary.keys())
#     Zs=template_dictionary[imfs[0]].Z[:, 0]
#     ages=template_dictionary[imfs[0]].age[0, :]

#     temp_lamdas=template_dictionary[imfs[0]].lam

#     t_mask = ((temp_lamdas > templates_lam_range[0]) & (temp_lamdas <templates_lam_range[1]))
#     sspNew, logLam_template, template_velscale = util.log_rebin(templates_lam_range, template_dictionary[imfs[0]].flam[0, 0, :][t_mask], velscale=velscale)
#     templates=np.empty((len(imfs), template_dictionary[imfs[0]].flam.shape[0], template_dictionary[imfs[0]].flam.shape[1], len(sspNew)))

#     #Make the template array
#     for i, imf in enumerate(imfs):
#         for j in range(len(Zs)):
#             for k in range(len(ages)):

#                 sspNew, logLam_template, _ = util.log_rebin(templates_lam_range, template_dictionary[imf].flam[j, k, :][t_mask], velscale=velscale)
#                 sspNew/=np.median(sspNew)
#                 templates[i, j, k, :]=sspNew

#         #templates[i, :]=template_dictionary[imfs].flam/(np.repeat(np.median(temp[imfs].flam, axis=2), temp[imfs].flam.shape[-1]).reshape(7, 12, -1))

#     #pPXF expects wavelength axis first
#     templates=templates.T

#     #temp_wave=np.logspace(np.log10(templates_lam_range[0]), np.log10(templates_lam_range[1]), templates.shape[0])

#     return templates, logLam_template, template_dictionary

# ################################################################################################################################################################


def NGC1277_CVD_read_in_data_CvD(file = '~/z/Data/IMF_Gold_Standard/n1277b_cen.dat', c=299792.458):
    ##############################################################################################################################################################

    # Read in the central spectrum from CvD
    #

    import os
    file=os.path.expanduser(file)
    lamdas, flux, errors, CvD_weights, inst_res=np.genfromtxt(file, unpack=True)
    #Redshift of NGC1277, used to get the initial velocity
    z=0.017044

    #The CvD data has a chip gap between 5625 and 7071 A. This screws up the fitting if you forget about it (obviously!).
    #This 'gap' array is insterted into the normalised flux and error arrays, before we mask it out in the fitting.

    flux_median=np.median(flux)

    # flux/=flux_median
    # errors/=flux_median

    chip_gap_start=lamdas[2424]
    chip_gap_stop=lamdas[2425]

    gap=np.ones(int(chip_gap_stop-chip_gap_start)) 
    flux=np.insert(flux/flux_median, 2424, gap)
    errors=np.insert(errors/flux_median, 2424, gap)

    quick_and_dirty_lamdas=np.linspace(lamdas[0], lamdas[-1], len(flux))
 
    #Mask to only fit some wavelength ranges?

    # lower=lamdas.min()
    # upper=lamdas.max()
    lower=4100
    upper=lamdas.max()

    assert (lower>=lamdas.min()) & (upper<=lamdas.max()), 'Lower and upper limits must be within the wavelength ranges of the data'

    

    lam_range_gal=np.array([lower, upper])
    mask=np.where((quick_and_dirty_lamdas>lower) & (quick_and_dirty_lamdas<upper))

    # print 'Lam Range Gal is {}'.format(lam_range_gal)

    flux=flux[mask]
    errors=errors[mask]


    #Log rebin them
    galaxy, logLam, velscale = util.log_rebin(lam_range_gal, flux)
    noise, _, _=util.log_rebin(lam_range_gal, errors, velscale=velscale)   


    #GoodPixels from pPXF
    goodpixels = np.arange(len(galaxy)) 
    

    return galaxy, noise, velscale, goodpixels, lam_range_gal, logLam

    ################################################################################################################################################################


def NGC1277_CVD_read_in_data_MN(file = '~/z/Data/IMF_Gold_Standard/n1277b_cen.dat', c=299792.458):
    ##############################################################################################################################################################

    # Read in the central spectrum from CvD
    #

    import os
    file=os.path.expanduser(file)
    lamdas, flux, errors=np.genfromtxt(file, unpack=True)
    #Redshift of NGC1277, used to get the initial velocity
    z=0.017044

    #The CvD data has a chip gap between 5625 and 7071 A. This screws up the fitting if you forget about it (obviously!).
    #This 'gap' array is insterted into the normalised flux and error arrays, before we mask it out in the fitting.

    flux_median=np.median(flux)

    # flux/=flux_median
    # errors/=flux_median
 
    lower=4150
    upper=lamdas.max()

    assert (lower>=lamdas.min()) & (upper<=lamdas.max()), 'Lower and upper limits must be within the wavelength ranges of the data'

    

    lam_range_gal=np.array([lower, upper])
    mask=np.where((lamdas>lower) & (lamdas<upper))

    # print 'Lam Range Gal is {}'.format(lam_range_gal)

    flux=flux[mask]
    errors=errors[mask]


    #Log rebin them
    galaxy, logLam, velscale = util.log_rebin(lam_range_gal, flux)
    noise, _, _=util.log_rebin(lam_range_gal, errors, velscale=velscale)   


    #GoodPixels from pPXF
    goodpixels = np.arange(len(galaxy)) 
    

    return galaxy, noise, velscale, goodpixels, lam_range_gal, logLam

    ################################################################################################################################################################


def NGC1277_CVD_read_in_data_SPV(file = '~/z/Data/IMF_Gold_Standard/n1277b_cen.dat', c=299792.458):
    ##############################################################################################################################################################

    # Read in the central spectrum from CvD
    #

    import os
    file=os.path.expanduser(file)
    lamdas, flux, variance, inst_res=np.genfromtxt(file, unpack=True)
    #Redshift of NGC1277, used to get the initial velocity
    z=0.017044
    errors=np.sqrt(variance)
    lamdas*=10**4

    #The CvD data has a chip gap between 5625 and 7071 A. This screws up the fitting if you forget about it (obviously!).
    #This 'gap' array is insterted into the normalised flux and error arrays, before we mask it out in the fitting.

    flux_median=np.median(flux)

    # flux/=flux_median
    # errors/=flux_median
 
    lower=lamdas.min()
    upper=lamdas.max()

    assert (lower>=lamdas.min()) & (upper<=lamdas.max()), 'Lower and upper limits must be within the wavelength ranges of the data'

    

    lam_range_gal=np.array([lower, upper])
    mask=np.where((lamdas>lower) & (lamdas<upper))

    # print 'Lam Range Gal is {}'.format(lam_range_gal)

    flux=flux[mask]
    errors=errors[mask]


    #Log rebin them
    galaxy, logLam, velscale = util.log_rebin(lam_range_gal, flux)
    noise, _, _=util.log_rebin(lam_range_gal, errors, velscale=velscale)   


    #GoodPixels from pPXF
    goodpixels = np.arange(len(galaxy)) 
    

    return galaxy, noise, velscale, goodpixels, lam_range_gal, logLam




# def NGC1277_CVD_read_in_data_MILES(file = 'Data/n1277b_cen.dat', c_light=299792.458):
#     ##############################################################################################################################################################

#     # Read in the central spectrum from CvD
#     #
#     lamdas, flux, errors, CvD_weights, inst_res=np.genfromtxt(file, unpack=True)
#     #Redshift of NGC1277, used to get the initial velocity
#     z=0.017044


#     #The CvD data has a chip gap between 5625 and 7071 A. This screws up the fitting if you forget about it (obviously!).
#     #This 'gap' array is insterted into the normalised flux and error arrays, before we mask it out in the fitting.

#     flux_median=np.median(flux)

#     chip_gap_start=lamdas[2424]
#     chip_gap_stop=lamdas[2425]

#     gap=np.ones(int(chip_gap_stop-chip_gap_start)) 
#     flux_edited=np.insert(flux/flux_median, 2424, gap)
#     errors_edited=np.insert(errors/flux_median, 2424, gap)

#     quick_and_dirty_lamdas=np.linspace(lamdas[0], lamdas[-1], len(flux_edited))



#     #Mask to only fit some wavelength ranges?

#     # lower=lamdas.min()
#     # upper=lamdas.max()
#     lower=lamdas.min()
#     upper=lamdas.max()

#     assert (lower>=lamdas.min()) & (upper<=lamdas.max()), 'Lower and upper limits must be within the wavelength ranges of the data'

#     lam_range_gal=np.array([lower, upper])
#     mask=np.where((quick_and_dirty_lamdas>lower) & (quick_and_dirty_lamdas<upper))

#     flux_edited=flux_edited[mask]
#     errors_edited=errors_edited[mask]

#     #Log rebin them
#     galaxy, logLam, velscale = util.log_rebin(lam_range_gal, flux_edited)
#     noise, _, _=util.log_rebin(lam_range_gal, errors_edited)   


#     #GoodPixels from pPXF
#     goodpixels = np.arange(len(galaxy))

#     #Points corresponding to the chipgap- remove from goodpixels
#     chip_gap_lamdas=np.where((logLam>np.log(chip_gap_start)) & (logLam<np.log(chip_gap_stop)))[0]
#     goodpixels = np.array([x for x in goodpixels if x not in chip_gap_lamdas])   
    

#     return galaxy, noise, velscale, goodpixels, lam_range_gal, logLam

#     ################################################################################################################################################################



# def read_in_templates(lam_range_gal, velscale, pad=500.0, imf_type='bi', verbose=True):
#     ################################################################################################################################################################
#     #Wavelength Range of the templates is much bigger than the galaxy, so we mask the templates to just a bit bigger than the galaxy wavelength range

#     lam_range_temp = [lam_range_gal[0]-pad, lam_range_gal[1]+pad]

#     #get the NaFe=0.0 templates and clip the very low metallicities  
#     templates, logLam_template, template_dictionary=prepare_Miles_templates(lam_range_temp, velscale, imf_type=imf_type, verbose=verbose)
#     #templates, logLam_template, template_dictionary=prepare_Miles_templates(lam_range_temp, velscale, template_dictionary=template_dictionary, imf_type='uni')
#     imfs=sorted(template_dictionary.keys())
#     Zs=template_dictionary[imfs[0]].Z[:, 0]
#     ages=template_dictionary[imfs[0]].age[0, :]

#     templates=templates[:, -6:, -3:, :]
#     Zs=Zs[-3:]

#     ages=ages[-6:]

#     NaFe3templates, logLam_template, NaFe3template_dictionary=prepare_Miles_templates(lam_range_temp, velscale, NaFe=0.3, imf_type=imf_type, verbose=verbose)
#     NaFe3templates=NaFe3templates[:, -6:, -3:, :]

#     NaFe6templates, logLam_template, NaFe6template_dictionary=prepare_Miles_templates(lam_range_temp, velscale, NaFe=0.6, imf_type=imf_type, verbose=verbose)
#     NaFe6templates=NaFe6templates[:, -6:, -3:, :]

#     NaFe9templates, logLam_template, NaFe9template_dictionary=prepare_Miles_templates(lam_range_temp, velscale, NaFe=0.9, imf_type=imf_type, verbose=verbose)
#     NaFe9templates=NaFe9templates[:, -6:, -3:, :]
#     ################################################################################################################################################################

#     all_templates=np.stack((templates, NaFe3templates, NaFe6templates, NaFe9templates), axis=4)
#     NaFes=np.array([0.0, 0.3, 0.6, 0.9])
#     Zs=np.array([-0.4, 0.0, 0.22])
#     imfs=[float(imf.strip('bi')) for imf in imfs]


#     linear_interp=si.RegularGridInterpolator(((logLam_template, ages, Zs, imfs, NaFes)), all_templates)

#     return linear_interp, logLam_template, lam_range_temp

#     ################################################################################################################################################################



# def NGC1277_CvD_set_up_emcee_parameters_MILES(file = 'n1277b_cen.dat', verbose=True):

    
#      # speed of light in km/s

#     galaxy, noise, velscale, goodpixels, lam_range_gal, logLam_gal=NGC1277_CVD_read_in_data_MILES(file=file, c=c_light)

#     linear_interp, logLam_template, lam_range_temp=read_in_templates(lam_range_gal, velscale, pad=500.0, verbose=verbose)

#     dv = c_light*np.log(lam_range_temp[0]/lam_range_gal[0])  # km/s


#     return [galaxy, noise, velscale, goodpixels, dv, linear_interp, logLam_template], logLam_gal


def NGC1277_CvD_set_up_emcee_parameters_CvD(file = '~/z/Data/IMF_Gold_Standard/n1277b_cen.dat', verbose=True):

    #positive_only_elems=['Cr+', 'Ni+', 'Co+', 'Eu+', 'Sr+', 'K+', 'V+', 'Cu+', 'as/Fe+']
    positive_only_elems=['as/Fe+']
    Na_elem=['Na']
    #normal_elems=['Ca', 'Fe', 'C', 'N', 'Ti', 'Mg']
    normal_elems=['Ca', 'Fe', 'C', 'N', 'Ti', 'Mg', 'Si']

    assert len(positive_only_elems)==1, 'Need to change the code if you want more than 1 positive element!'

    elements=(positive_only_elems, Na_elem, normal_elems)

    fit_wavelengths=np.array([[4000, 4700], [4700, 5500], [8000,  9100], [9600, 10150]])
    

    galaxy, noise, velscale, goodpixels, lam_range_gal, logLam_gal=NGC1277_CVD_read_in_data_CvD(file=file, c=c_light)


    pad=500.0
    lam_range_temp = [lam_range_gal[0]-pad, lam_range_gal[1]+pad]
    linear_interp, logLam_template =prepare_CvD_interpolator(lam_range_temp, velscale, verbose=True)
    correction_interps, logLam_template=prepare_CvD_correction_interpolators(lam_range_temp, velscale, elements, verbose=True)


    dv = c_light*np.log(lam_range_temp[0]/lam_range_gal[0])  # km/s

    #ndim is number of elements plus V, Sigma plus age, IMF and Z
    ndim=len(positive_only_elems)+len(Na_elem)+len(normal_elems)+2+3


    return [galaxy, noise, velscale, goodpixels, dv, linear_interp, correction_interps, logLam_template, logLam_gal, fit_wavelengths], logLam_gal, ndim

def NGC1277_CvD_set_up_emcee_parameters_MN(file = '~/z/Data/IMF_Gold_Standard/NGC1277_RAD0.00_PPXF_NEW.cxt', verbose=True):

    #positive_only_elems=['Cr+', 'Ni+', 'Co+', 'Eu+', 'Sr+', 'K+', 'V+', 'Cu+', 'as/Fe+']
    positive_only_elems=['as/Fe+']
    Na_elem=['Na']
    #normal_elems=['Ca', 'Fe', 'C', 'N', 'Ti', 'Mg']
    normal_elems=['Ca', 'Fe', 'C', 'N', 'Ti', 'Mg', 'Si']

    assert len(positive_only_elems)==1, 'Need to change the code if you want more than 1 positive element!'

    elements=(positive_only_elems, Na_elem, normal_elems)

    fit_wavelengths=np.array([[4000, 4700], [4700, 5500], [5500, 6800], [8000,  8700]])
    

    galaxy, noise, velscale, goodpixels, lam_range_gal, logLam_gal=NGC1277_CVD_read_in_data_MN(file=file, c=c_light)


    pad=500.0
    lam_range_temp = [lam_range_gal[0]-pad, lam_range_gal[1]+pad]
    linear_interp, logLam_template =prepare_CvD_interpolator(lam_range_temp, velscale, verbose=True)
    correction_interps, logLam_template=prepare_CvD_correction_interpolators(lam_range_temp, velscale, elements, verbose=True)


    dv = c_light*np.log(lam_range_temp[0]/lam_range_gal[0])  # km/s

    #ndim is number of elements plus V, Sigma plus age, IMF and Z
    ndim=len(positive_only_elems)+len(Na_elem)+len(normal_elems)+2+3


    return [galaxy, noise, velscale, goodpixels, dv, linear_interp, correction_interps, logLam_template, logLam_gal, fit_wavelengths], logLam_gal, ndim


################################################################################


def NGC1277_CvD_set_up_emcee_parameters_SPV(file = '~/z/Data/IMF_Gold_Standard/SPV_NGC1277.dat', verbose=True):

    fit_wavelengths=np.array([[6300, 10412]])

    positive_only_elems=['as/Fe+']#, 'as/Fe+']
    Na_elem=['Na']
    normal_elems=['Ca', 'Fe', 'Ti', 'Mg']

    #Try not fitting age, since we're not very sensitive to it in the SWIFT range.
    population_params=['Z', 'IMF']
    kinematic_parmas=['V', 'Sigma']

    #assert len(positive_only_elems)==1, 'Need to change the code if you want more than 1 positive element!'

    elements=(positive_only_elems, Na_elem, normal_elems)
    

    galaxy, noise, velscale, goodpixels, lam_range_gal, logLam_gal=NGC1277_CVD_read_in_data_SPV(file=file, c=c_light)


    pad=500.0
    lam_range_temp = [lam_range_gal[0]-pad, lam_range_gal[1]+pad]

    
    linear_interp, logLam_template =prepare_CvD_interpolator(lam_range_temp, velscale, verbose=True)
    correction_interps, logLam_template=prepare_CvD_correction_interpolators(lam_range_temp, velscale, elements, verbose=True)


    dv = c_light*np.log(lam_range_temp[0]/lam_range_gal[0])  # km/s

    #ndim is number of elements plus V, Sigma plus age, IMF and Z
    ndim=len(positive_only_elems)+len(Na_elem)+len(normal_elems)+len(kinematic_parmas)+len(population_params)


    return [galaxy, noise, velscale, goodpixels, dv, linear_interp, correction_interps, logLam_template, logLam_gal, fit_wavelengths], logLam_gal, ndim


################################################################################




