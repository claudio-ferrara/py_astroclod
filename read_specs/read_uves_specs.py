import numpy as np

from astropy.io import fits


def read_uves_spec(filename):
    
    """
    ---Function read_uves_spec---
    ---This function reads UVES data products of "ADP" type---
    ---Input parameters:
    -----filename - name of fits file to be read
    ---Returns:
    -----w - array of wavelengths
    -----f - array of fluxes
    -----ef - array of errors on fluxes
    ---Last modification - January 24th, 2023
    """

    #-----reading fits file-----
    data = fits.getdata(filename)

    w = np.transpose(data['WAVE'])
    f = np.transpose(data['FLUX_REDUCED'])
    ef = np.transpose(data['ERR_REDUCED'])

    return w, f, ef

