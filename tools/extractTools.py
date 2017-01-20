import numpy as np
import numpy.ma as ma
from astropy.io import fits

import opt_extract_functs as optex

#import settings

import argparse
#from python_utils import *
import pdb



def optimal_extraction(x, y, binNum, cube, varcube):

    """A function to take the output from spatial binning (e.g Voronoi/Sector binning) and optimally extract a spectrum from each bin. 




    NEEDS WORK
    """

    binNum=binNum.astype(int)
    number_of_bins=len(np.unique(binNum)[np.unique(binNum)!=-1])

    d1, d2, d3=cube.shape


    print "\n\nThere are {} bins".format(number_of_bins)

    for i in range(number_of_bins):


            #Create a mask with True where the bin indices are, false everywhere else
            inds=np.where(binNum==(i))
            x_inds=x[inds].astype(int)
            y_inds=y[inds].astype(int)

            aperture_indices=[y_inds, x_inds]

            #True corresponds to masked
            mask=np.ones_like(cube[0, :, :]).astype(bool)
            mask[y_inds, x_inds]=False     



            aperture_mask = np.repeat(mask[np.newaxis, :, :], d1, axis=0)




            print("Making Spectra {}".format(i+1))
            cntr=1.0
            ##########################################################################################
            
            # plt.imshow(mask)
            # plt.show()

            frame=np.median(cube, axis=0)


            spec=np.nansum(np.nansum(cube*aperture_mask.astype(int), axis=2), axis=1)
            nspec=np.sqrt(np.nansum(np.nansum(var*aperture_mask.astype(int), axis=2), axis=1))


            #import pdb; pdb.set_trace()

            temp=cube*aperture_mask.astype(int)
            temp2=var*aperture_mask.astype(int)
            sig_temp=temp[temp>0.1]
            noise_temp=temp2[temp2>0.1]

            if args.plot:
                import matplotlib.pyplot as plt
                plt.imshow(frame, alpha=0.5)
                plt.imshow(mask, alpha=0.5)
                plt.show()

            # plt.plot(spec)
            # plt.figure()
            # plt.plot(nspec)
            # plt.show()
            #print np.shape(np.where(mask==True)[0])
            #print np.mean(cube*aperture_mask.astype(int))

            SNR=np.mean(np.nansum(np.nansum((cube*aperture_mask.astype(int)), axis=2), axis=1)/np.sqrt(np.nansum(np.nansum((var*aperture_mask.astype(int)), axis=2), axis=1)))
            print "SNR IS {}".format(SNR)
            print "Median Signal is {}".format(np.median(sig_temp))
            print "Median noise is {}".format(np.median(noise_temp))
            
            #pdb.set_trace()
            ##########################################################################################
            print "Extracted Spectrum"
            weights=False
            verbose=False
            w=1
            sigma=50
            while cntr!=0.0:
            #for c in range(5):
                print("Wavelength Iteration {}".format(w))

                if not weights:
                    print("\n\nCONSTRUCTING SPATIAL PROFILE\n")
                    print lamdas.shape
                    
                    mcube=optex.wavelegth_fitter(blank_cube, cube, var, aperture_mask, lamdas, l_order, aperture_indices, plot=False, verbose=verbose)
                else:
                    print("\n\nUSING WEIGHTS FROM FILE {}".format(weights))

                print("\n\nREVISING VARIANCE ESTIMATES\n")

                #pdb.set_trace()


                var=optex.variance(mcube, spec, skycube, 0.0, gain)
                varbad=np.where(~np.isfinite(var))
                var[varbad]=np.inf
                print("\n\nMASKING COSMICS/BAD PIXELS\n")


                aperture_mask, cntr=optex.sigma_clipping(cube, var, mcube, spec, aperture_mask, sigma, aperture_indices, verbose=True)
                var[~aperture_mask]=np.inf
                cube[~aperture_mask]=0



                print("\n\nEXTRACTING OPTIMAL SPECTRUM\n")

                spec=np.sum(np.sum(aperture_mask*cube*mcube/var, axis=2), axis=1)/np.sum(np.sum(aperture_mask*mcube*mcube/var, axis=2), axis=1)
                vspec=np.sum(np.sum(aperture_mask*mcube, axis=2), axis=1)/np.sum(np.sum(aperture_mask*mcube*mcube/var , axis=2), axis=1)
                w+=1    
                print(cntr)


            #import pdb; pdb.set_trace()

            spectra[:, i]=spec
            noisespectra[:, i]=np.sqrt(vspec)


            #pdb.set_trace()
            print("DONE")


     
            for (spec, nspec) in zip(spectra.T, noisespectra.T):

                print "Median S/N per A is {}".format(np.median(spec/nspec))

            if args.plot:
                import matplotlib.pyplot as plt
                plt.figure()
                plt.plot(lamdas, spectra)
                plt.figure()
                cm = plt.get_cmap('inferno')
                for j, (spec, nspec) in enumerate(zip(spectra.T, noisespectra.T)):

                    plt.plot(lamdas, spec/nspec, c=cm(1.0*j/13.0), linewidth=2.0)
                plt.show()

def simple_extraction(x, y, binNum, cube, verbose=False):

    """Extract a spectrum by simply summing all the pixels in a bin at each wavelength value

    Ignore nans by making a global 'nan' mask and using the bitwise and with that and the bin mask
    """

    binNum=binNum.astype(int)
    number_of_bins=len(np.unique(binNum)[np.unique(binNum)!=-1])

    d1, d2, d3=cube.shape

    spectra=np.empty((number_of_bins, d1))

    
    #Mask all nans in the cube
    nan_mask=np.zeros_like(cube).astype(bool)
    nan_values=np.where(~np.isfinite(cube))

    nan_mask[nan_values]=True

    for i in range(number_of_bins):
        if verbose:
            print "Extracting spectrum {} of {}".format(i, number_of_bins)

        #Create a mask with True where the bin indices are, false everywhere else
        inds=np.where(binNum==(i))
        x_inds=x[inds].astype(int)
        y_inds=y[inds].astype(int)

        aperture_indices=[y_inds, x_inds]

        #True corresponds to masked
        mask=np.ones_like(cube[0, :, :]).astype(bool)
        mask[y_inds, x_inds]=False



        aperture_mask = np.repeat(mask[np.newaxis, :, :], d1, axis=0)

        final_mask=np.bitwise_or(aperture_mask,nan_mask)

        masked_cube=ma.array(cube, mask=final_mask)


        spectra[i :]=ma.sum(ma.sum(masked_cube, axis=2), axis=1)

    return spectra, nan_mask


def save_array_as_fits(spectra, filename, header=None, clobber=False, verbose=False):

    """
    Save a 2D array of spectra (of the form [n_spectra, lamda]) as a multi extension fits file
    """

    newHDU=fits.HDUList()




    for j in range(spectra.shape[0]):


        if header is not None:
            header['BinNum']=j

        newHDU.append(fits.ImageHDU(spectra[j, :], header=header))

    newHDU.writeto(filename, clobber=clobber)

    if verbose:
        print "Saved to {}".format(filename)













