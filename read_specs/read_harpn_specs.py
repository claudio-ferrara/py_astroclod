import numpy as np
import pandas as pd

from astropy.io import fits

def extract_cards_vals_from_header(h):
        
    """
    ---Function extract_cards_vals_from_header---
    ---This function extracts cards (keywords), header values and comments from an astropy header h---
    ---Input parameters:
    -----h - astropy header of fits file
    ---Returns:
    -----carte - array of header cards
    -----hvals - array of header vals
    -----comments - array of comments
    ---Last modification - January 24th, 2023
    """

    carte = []
    hvals = []
    comments = []
    for cards in h: 
        carte.append(cards) 
        hvals.append(h[cards])
        comments.append(h.comments[cards]) 
        
    return carte, hvals, comments


def read_harpn_s1d(filename):
    
    """
    ---Function read_harpn_s1d---
    ---This function reads HARPS-N data products of "s1d" type---
    ---Input parameters:
    -----filename - name of fits file to be read
    ---Returns:
    -----w - array of wavelengths
    -----f - array of fluxes
    ---Last modification - January 24th, 2023
    """

    #-----reading fluxes and header from fits file-----
    data = fits.getdata(filename)
    header = fits.getheader(filename)
    #-----constructing wavelengths array-----
    wavelength = np.linspace(header['CRVAL1'],header['CRVAL1'] + header['NAXIS1']*header['CDELT1'],header['NAXIS1'],endpoint=True)    
    
    return wavelength, data



def read_harpn_e2ds(filename):

    """
    ---Function read_harpn_e2ds---
    ---This function reads HARPS-N data products of "e2ds" type---
    ---It returns the 2D arrays of wavelength and fluxes---
    ---Input parameters:
    -----filename - name of fits file containing the spectrum
    ---Returns:
    -----lams - 2D array of wavelengths
    -----flxs - 2D array of fluxes
    ---Last modification - January 24th, 2023
    """


    #----reading fluxes from fits file----
    im = fits.getdata(filename)
    #----getting the length in pixes of the orders----
    npix= len(im[0])
    #----getting cards and keywords from the header with apposite function----
    carte,hvals,cmmnts = extract_cards_vals_from_header(fits.getheader(filename))
    
    head = pd.DataFrame({'keywords':carte,'values':hvals,'comments':cmmnts})
    #----getting the degree of the polynomial
    head['indexes'] = head['keywords'].str.find("CAL TH DEG LL")
    deg = head[head['indexes'] >= 0]['values'].to_numpy()[0]
    #----getting and sorting polynomial coefficients----
    head['indexes'] = head['keywords'].str.find("CAL TH COEFF LL")
    coeffs = head[head['indexes'] >= 0]['values'].to_numpy()
    coeff_numbers = head[head['indexes'] >= 0]['keywords']
    coeff_numbers = coeff_numbers.str.slice(start=23).to_numpy(dtype=int)
    srt = np.argsort(coeff_numbers) 
    coeff_numbers = coeff_numbers[srt]
    coeffs = np.flip(coeffs[srt])
    
    #-----constructing wavelength arrays-----
    #x = np.arange(1,npix+1)
    x = np.arange(npix)
    lam = []
    for i in range(0,len(coeffs),deg+1): 
        lam.append(np.polyval(coeffs[i:i+deg+1],x))
    
    lam = np.array(lam)
    flx = np.flip(im,axis=0)

    return lam, flx