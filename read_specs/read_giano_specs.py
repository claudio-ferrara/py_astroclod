import numpy as np

from astropy.io import fits


def read_giano_ms1d(filename):

    """
    ---Function read_giano_ms1d---
    ---This function reads GIANO-B data products of "ms1d" type---
    ---ms1d are 2D spectra not corrected for the barycentric hearth radial velocity----
    ---Wavelengths are in NANOMETERS and VACUUM---
    ---Input parameters:
    -----filename - name of fits file to be read
    ---Returns:
    -----lam - array of wavelengths
    -----flx - array of fluxes
    -----snr - array of signal-to-noise ratio
    -----ords - array of order numbers
    ---Last modification - January 24th, 2023
    """

    #----getting fluxes from fits file
    im = fits.getdata(filename)

    #----getting wavelengths, fluxes, snr and order numbers
    for i in range(len(im['WAVE'])):
        srt = np.argsort(im['WAVE'][i])
        im['WAVE'][i] = im['WAVE'][i][srt]
        im['FLUX'][i] = im['FLUX'][i][srt]
        im['SNR'][i] = im['SNR'][i][srt]
    
    lam = im['WAVE']
    flx = im['FLUX']
    snr = im['SNR']
    ords = im['ORDER']

    return lam, flx, snr, ords


def read_giano_s1d(filename):
    
    """
    ---Function read_giano_s1d---
    ---This function reads GIANO-B data products of "s1d" type---
    ---s1d are 1D re-sampled spectra corrected for the barycentric hearth radial velocity----
    ---Wavelengths are in NANOMETERS and VACUUM---
    ---Input parameters:
    -----filename - name of fits file to be read
    ---Returns:
    -----w - array of wavelengths
    -----f - array of fluxes
    ---Last modification - January 24th, 2023
    """

    #----getting fluxes and header from fits file
    f = fits.getdata(filename)
    header = fits.getheader(filename)

    #----creating wavelength arrays----
    w = np.linspace(header['CRVAL1'],header['CRVAL1'] + header['NAXIS1']*header['CDELT1'],header['NAXIS1'],endpoint=True)
    
    return w, f