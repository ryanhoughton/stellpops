import numpy as np
import numpy.ma as ma
from astropy.io import fits

import opt_extract_functs as optex

#import settings

import argparse
#from python_utils import *
import pdb
import time as T
import datetime



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

def simple_extraction(x, y, binNum, cube, verbose=False, type='sum'):

    """Extract a spectrum by simply summing all the pixels in a bin at each wavelength value

    Ignore nans by making a global 'nan' mask and using the bitwise and with that and the bin mask.

    Can either add spectra by summing or median combining (very simply! No continuum adition or anything)
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

        if type=='sum':
            spectra[i :]=ma.sum(ma.sum(masked_cube, axis=2), axis=1)
        elif type=='median':
            spectra[i :]=ma.median(ma.median(masked_cube, axis=2), axis=1)
        else:
            return NameError('Type of combination not understood')

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



def plot_specs(i):
    fig, axs=plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(20, 14))
    axs[0].plot(sig3spectra[:, i], c='k', linewidth=2.0)
    axs[0].plot(sig500spectra[:, i], c='r')

    axs[1].plot(sig3spectra[:, i]/sig500spectra[:, i], c='k')

    axs[1].set_title('3 Sig Clipping/ 500 sig Clipping')

    return axs



def extract_spectra(cube, varcube, bin_text_file, offset=None, plot=False):

    """
    Given a SWIFT cube, a variance cube and a text file of bins (from sector binning- need to change the code to use voronoi bins), extract all the spectra to numpy arrays.
    Use make_fits_tables to extract and save to a fits table

    Inputs:
        cube: a 3D numpy array of a datacube. Shape is dlam, dx, dy
        varcube: a 3D numpy array of a variance cube. Should be the same shape as varcube
        bin_text_file: a text file with a 1D array of bin numbers. Should be of shape dx*dy. 
    Optional Inputs:
        offset: an integer pixel shift between the bin text file and the cube. Useful when we want to extract spectra from an individual O-S cube but the bin text file
                refers to the combined cube of many O-S observations
        plot: boolean. If true, plot the bins over the cube to make sure that they lie on top of each other (i.e check your 'offset' is sensible)

    Returns:
        lamdas: an array np.arange(6300, 6300+d1) of wavelengths for the SWIFT observation
        spectra: an array of spectra: shape [d1, n_bins]
        noise_spectra: an array of NOISE spectra: shape [d1, n_bins]
    """

    #Get header info for wavelength masks
    d1, d2, d3=np.shape(cube)
    """
    crval=header["CRVAL3"]
    crpix=header["CRPIX3"]
    cdelt=header["CDELT3"]


    lamdas=np.array([(l+1-crpix)*cdelt + crval for l in range(d1)])*10**4
    """
    lamdas=np.arange(6300, 6300+d1)

    # low_lam, high_lam, z=settings.return_values(GalName, index_to_measure)


    # d_lam=high_lam-low_lam
    # lamda_mask=np.where((lamdas>=low_lam) & (lamdas<=high_lam))[0]

    # #Find the indices covered by the lamda mask.
    # low_index=np.where((lamdas>=low_lam))[0][0]
    # high_index=np.where(lamdas<=high_lam)[0][-1]

    # #Clip the cubes and lamda array
    # lamdas=lamdas[low_index:high_index]


    # cube=cube[low_index:high_index,:,:]
    # var=var[low_index:high_index,:,:]
    # skycube=skycube[low_index:high_index,:,:]


    #Get the shape of the cube
    d1, d2, d3=np.shape(cube)

    #Make a blank cube to be used as a mask
    blank_cube=np.zeros_like(cube)



    #Load our bin text file
    y, x, binNum=np.loadtxt(bin_text_file, unpack=1)
    if offset is not None:
        x_offset, y_offset=offset
        x+=int(np.round(x_offset))
        y+=int(np.round(y_offset))

    # if GalName=="NGC1277":
    #     dims=np.loadtxt("/Volumes/SPV_SWIFT/SectorBinning/param_files/{}.txt".format(GalName))
    #     print "Loading parameters from {}.txt".format(GalName)
    #     x1=int(dims[0])
    #     x2=int(dims[1])
    #     y1=int(dims[2])
    #     y2=int(dims[3])

    #     y=y+y1
    #     x=x+x1

    #Change binNum to be int type and find the number of bins.
    #The Number of bins is number of unique numbers in binNum, ignoring any bins which are assigned -1
    binNum=binNum.astype(int)
    number_of_bins=len(np.unique(binNum)[np.unique(binNum)!=-1])

    print "\n\nThere are {} bins".format(number_of_bins)
    #Create empty arrays to be filled with the spectra and noise spectra.
    #Note- Noise spectra, not variance!
    spectra=np.zeros([len(lamdas), number_of_bins])
    noisespectra=np.zeros([len(lamdas), number_of_bins])

    # def find_centre_bin(bins):
    #     centre=np.where(bins==1.0)
    #     return int(np.median(centre[0])), int(np.median(centre[1]))

    # import matplotlib.pyplot as plt 
    # plt.figure()
    # plt.imshow(cube[2000, :, :], vmax=400.0, vmin=0.0)
    # plt.imshow(binNum.reshape(len(np.unique(y)), len(np.unique(x))), interpolation='nearest', alpha=0.6, cmap='viridis', extent=[x.min(), x.max(), y.min(), y.max()])
    # plt.show()
    # import pdb; pdb.set_trace()


    l_order=d1/100

    gain=0.98

    ORIGINAL=np.copy(cube)
    ORIGINALVAR=np.copy(varcube)
    ORIGINALBLANKCUBE=np.copy(blank_cube)



    for i in range(number_of_bins):


        cube=np.copy(ORIGINAL)
        var=np.copy(ORIGINALVAR)
        blank_cube=np.copy(ORIGINALBLANKCUBE)

        inds=np.where(binNum==(i+1))




        x_inds=x[inds].astype(int)
        y_inds=y[inds].astype(int)



        aperture_indices=[y_inds, x_inds]

        mask=np.zeros_like(cube[10, :, :]).astype(bool)
        mask[y_inds, x_inds]=True



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

        if plot:
            import matplotlib.pyplot as plt
            plt.imshow(frame, alpha=0.5)
            plt.imshow(mask, alpha=0.5)
            plt.show()
            import pdb; pdb.set_trace()

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
        sigma=500.0
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


            var=optex.variance(mcube, spec, var, 0.0, gain)
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

    return lamdas, spectra, noisespectra

def make_fits_tables(cubepath, varcubepath, bin_text_file, outfilename, offset=None):

    """
    A wrapper around extract_spectra. Run that function and then save the results to a fits table

    Inputs:
        cubepath: a filepath to a datacube
        varcube path: a path to a variance cube
        bin_text_file: a text file with a bin number for each pixel in the cube
        outfilename: a path to save the fits table to

    Optional Inputs:
        offsets: a 2-element-list/tuple/1D array of x, y offsets to pass to extract_spectra.

    Returns:
        None
    """

    cubeHDU=fits.open(cubepath)
    varcubeHDU=fits.open(varcubepath)
    cube=cubeHDU[0].data
    varcube=varcubeHDU[0].data
    lamdas, spectra, noise_spectra=extract_spectra(cube, varcube, bin_text_file, offset=offset, plot=False)

    n_lamdas, n_spectra=spectra.shape
    newHDUlist = [fits.PrimaryHDU()]

    for j in range(n_spectra):

        #Chop off the first 50 pixels as they're usually bad
        print 'REMOVING THE FIRST 50 PIXELS'

        col1 = fits.Column(name='lamda', format='D', array=lamdas[50:])
        col2 = fits.Column(name='flux', format='D', array=spectra[50:, j])
        col3 = fits.Column(name='noise', format='D', array=noise_spectra[50:, j])

        cols = fits.ColDefs([col1, col2, col3])

        tbhdu = fits.BinTableHDU.from_columns(cols)


        #Write the Cube's header to the table header
        tbhdu.header.extend(cubeHDU[0].header)   


        tbhdu.header['MJD-OBS']=float(cubeHDU[0].header['MJD'])

        #Get the time in seconds
        time=cubeHDU[0].header['UTC-TCS'].split()[1]
        x=T.strptime(time,'%H:%M:%S.%f')
        t_seconds=datetime.timedelta(hours=x.tm_hour,minutes=x.tm_min,seconds=x.tm_sec).total_seconds()

        tbhdu.header['TM-START']=t_seconds

        amend=float(cubeHDU[0].header['AMEND'])
        telalt=np.arccos(1.0/amend)*(180/np.pi)
        tbhdu.header['TEL-ALT']=telalt

        newHDUlist.append(tbhdu)

    fits.HDUList(newHDUlist).writeto(outfilename, clobber=True)







