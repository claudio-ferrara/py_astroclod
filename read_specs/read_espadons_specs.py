import numpy as np

from astropy.io import fits


def cut_esp_vect(lo,fo):

    """
    ---Function cut_esp_vect---
    ---This function generate a series of (l,f) orders from---
    --- ESPADONS "i" data products
    ---Input parameters:
    -----lo - array of unsorted wavelengths
    -----fo - array of (unsorted) fluxes
    ---Returns:
    -----l - 2D array of wavelengths
    -----f - 2D array of fluxes
    ---Last modification - January 24th, 2023
    """

    #----creating "d_lo" shift array----
    svec = lo-np.roll(lo,1)
    #----searching points where wavelength difference i large
    tmp = np.where(np.abs(svec) > 0.01)[0]
    svec = svec[tmp]

    #----composing 2D arrays
    l = []
    f = []
    for i in range(len(tmp)-2):
        l.append(lo[tmp[i]:tmp[i+1]-1])
        f.append(fo[tmp[i]:tmp[i+1]-1])
        

    #-----(l,f) - l and f are lists - they cannot be turned into matrixes
    return l,f


def read_esp_2D(filename):

    """
    ---Function read_esp_2D---
    ---This function reads ESPADONS "i" data products
    ---Input parameters:
    -----filename - name of fits file
    ---Returns:
    -----lam - 2D array of wavelengths
    -----flx - 2D array of fluxes
    ---Last modification - January 24th, 2023
    """


    im = fits.getdata(filename)
    lam,flx = cut_esp_vect(im[0],im[1])
    #w,f = totspectrum_2_esp(lam,flx)
    return lam, flx