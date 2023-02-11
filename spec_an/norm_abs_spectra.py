import numpy as np
import matplotlib.pyplot as plt

def med_rej_out(y,fact):

    """
    ---Function med_rej_out---
    ---This function takes as input a vector y of values and returns
    ---residual values after n-sigma clipping
    -----Input parameters:
    --- y - array of numerical values
    --- fact - number of sigma for sigma clipping
    -----Output parameters:
    --- tmp - array of indexes of "good" values, surviving sigma clipping
    ---
    """

    #-----calculating mean and standard deviation-----
    med = np.mean(y)   
    sig = np.std(y)
    #-----condition for sigma clipping-------
    tmp = np.where((y > (med-fact*sig)) & (y < (med + fact*sig)))[0]
    #------tmp is the array of good value------

    return tmp
    
def med_reject_outliers(y,fact = 3,niter=5):

    """
    ---Function med_reject_outliers-----
    ---This function performs n-sigma clipping of an array y n_iter times
    -----Input parameters:
    --- y - array of numerical values
    --- fact - number of sigma for sigma clipping (default = 3)
    --- niter - number of iteration for sigma clipping (default = 5)
    -----Output parameters:
    --- tmp - array of indexes of "good" values, surviving sigma clipping
    """

    ynew = y

    #-----first iteration of sigma clipping-----
    tmp = med_rej_out(ynew,fact)
    ynew = ynew[tmp]

    #------iterating over niter
    for i in range(1,niter):

        tmp2 = med_rej_out(ynew,fact)
        ynew = ynew[tmp2]
        #-----updating tmp array
        tmp = tmp[tmp2]

    return tmp

def norm_order(wocc,focc,deg=6,fact=1,niter=3,niter_fit=2,ifplot=False):

    """
    ---Function norm_order-----
    ---This functions does continuum normalization of an "order" of a high 
    ---resolution spectrum (about 50 Angstroms)
    ---Normalization is carried out using sigma clipping on the derivative of the spectrum (np.gradient)
    ---and then doing a least squares polynomial fitting using lower sigma clipping to exclude possible residual lines
    ---THIS ROUTINES IS DESIGN FOR ABSORPTION SPECTRA WITH METAL LINES
    ---IT SHOULD NOT BE USED FOR BALMER LINES NORMALIZATION
    -----Input parameters:
    --- wocc - array of wavelengths
    --- focc - array of fluxes
    --- deg - polynomial degree for continuum fit (default = 6)
    --- fact - factor for derivative sigma clipping (default = 1)
    --- niter - number of derivative sigma clipping iterations (default = 3)
    --- niter_fit - number of fit sigma clipping iterations (default = 2)
    --- ifplot - set as "True" if you want to plot the spectrum with the fitted continuum superimposed
    -----Output parameters:
    --- focn - normalized flux (corresponding to wocc wavelengths)
    """
    #------removing nans--------
    tmp = np.where(np.isnan(focc) == False)[0]

    if len(tmp) != 0:
        woc = wocc[tmp]
        foc = focc[tmp]
    else:
        return focc

    #----calculating the derivative using numpy----
    fd = np.gradient(foc)

    #-----performing sigma clipping on the derivative-----
    #-----preliminary identification of continuum---------
    tmp = med_reject_outliers(fd,fact=fact,niter=niter)

    wcont = woc[tmp]
    fcont = foc[tmp]
    xbar = wcont.mean()
    xcont = wcont-xbar
    
    #-----beginning polynomial fitting with outlier rejection of residuals----
    coeffs = np.polynomial.Legendre.fit(xcont,fcont,deg = deg)
    
    fcont_fit = coeffs(xcont)

    #------I am fitting (and rejecting outliers) niter_fit times-----
    for i in range(1,niter_fit):
        resfit = fcont-fcont_fit
        tmp = np.where(resfit > 0)[0]
        xcont = xcont[tmp]
        fcont = fcont[tmp]
        coeffs = np.polynomial.Legendre.fit(xcont,fcont,deg=deg)
        fcont_fit = coeffs(xcont)

    
    foc_cont = coeffs(woc-xbar)
    focc_cont = coeffs(wocc-xbar)

    #-------plotting if ifplot==True----
    if ifplot==True:
        plt.close()
        plt.plot(woc,foc)
        plt.plot(woc,foc_cont,'red')
        plt.show()

    focn = focc/focc_cont

    return focn