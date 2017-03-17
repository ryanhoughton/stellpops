import scipy.optimize as opt
import numpy as np 
import matplotlib.pyplot as plt 

def convolve_variable_width(a, sig, prec=1.):
    '''
    approximate convolution with a kernel that varies along the spectral
        direction, by stretching the data by the inverse of the kernel's
        width at a given position
    N.B.: this is an approximation to the proper operation, which
        involves iterating over each pixel of each template and
        performing ~10^6 convolution operations
    Parameters:
     - a: N-D array; convolution will occur along the final axis
     - sig: 1-D array (must have same length as the final axis of a);
        describes the varying width of the kernel
     - prec: precision argument. When higher, oversampling is more thorough
    '''

    assert (len(sig) == a.shape[-1]), '\tLast dimension of `a` must equal \
        length of `sig` (each element of a must have a convolution width)'

    sig0 = sig.max()  # the "base" width that results in minimal blurring
    # if prec = 1, just use sig0 as base.

    n = np.rint(prec * sig0/sig).astype(int)
    print n.min()
    # print n
    print '\tWarped array length: {}'.format(n.sum())
    # define "warped" array a_w with n[i] instances of a[:,:,i]
    a_w = np.repeat(a, n, axis=-1)
    # now a "warped" array sig_w
    sig_w = np.repeat(sig, n)

    # define start and endpoints for each value
    nl = np.cumsum(np.insert(n, 0, 0))[:-1]
    nr = np.cumsum(n)
    # now the middle of the interval
    nm = np.rint(np.median(np.column_stack((nl, nr)), axis=1)).astype(int)

    # print nm

    # print a_w.shape, sig_w.shape # check against n.sum()

    # now convolve the whole thing with a Gaussian of width sig0
    print '\tCONVOLVE...'
    # account for the increased precision required
    a_w_f = np.empty_like(a_w)
    # have to iterate over the rows and columns, to avoid MemoryError
    c = 0  # counter (don't judge me, it was early in the morning)
    # #for i in range(a_w_f.shape[0]):
    #     #for j in range(a_w_f.shape[1]):
    #         '''c += 1
    #         print '\t\tComputing convolution {} of {}...'.format(
    #             c, a_w_f.shape[0] * a_w_f.shape[1])'''
    #         #a_w_f[i, j, :] = ndimage.gaussian_filter1d(a_w[i, j, :], prec*sig0)
    a_w_f = ndimage.gaussian_filter1d(a_w, prec*sig0)
    # print a_w_f.shape # should be the same as the original shape

    # and downselect the pixels (take approximate middle of each slice)
    # f is a mask that will be built and applied to a_w_f

    # un-warp the newly-convolved array by selecting only the slices
    # in dimension -1 that are in nm
    a_f = a_w_f[nm]

    return a_f


def gauss(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def fit_gaussian_to_skyline(lamdas, skyspec, line, plot=False):


    start, stop=line

    assert start<stop, "Your starting wavelength must be below your stopping one"

    mask=np.where((lamdas>start) & (lamdas < stop))

    assert len(mask[0]) >0, "The mask has no elements!"

    x=lamdas[mask]
    y=skyspec[mask]-skyspec[mask][0]



    popt,pcov = opt.curve_fit(gauss,x,y,p0=[1,np.mean(x),10])

    if plot:
        plt.plot(x, y, c='k')
        plt.plot(x, gauss(x, *popt), c='r')

        plt.show()


    sig=np.abs(popt[-1])

    return sig


def get_MUSE_skylines():

    skylines=np.array([[5560, 5590], [5927, 5940],  [6280, 6291], [6340, 6380], [6572, 6581], [6942, 6955], [7333, 7347], [7395, 7405], [7742, 7753], [7811, 7827], [7985, 7999], [8374, 8390], [8458, 8470], [8497, 8508], [8783, 8797], [8879, 8890], [8950, 8964], [9096, 9110]])
    return skylines













