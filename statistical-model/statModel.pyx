"""
Cythonised code for the statistical model.

To use with the python scripts, please build this extension first:

$ python setup.py build_ext --inplace

Author: Ross <ross dot warren at pm dot me>
Date:   04/12/2019
"""

import time
from libc.math cimport exp, sqrt, pow, pi
from scipy.integrate import quad
from scipy.optimize import curve_fit
import numpy as np
cimport numpy as np


cdef double gaussianDoS(float E, float N0, float E0, float sigma):
    """Gaussian density of states."""
    cdef double g, Etemp
    Etemp = E - E0
    Etemp = pow(Etemp, 2)
    g = (N0 / (sqrt(2 * pi) * sigma)) * exp(-(Etemp) / (2 * sigma ** 2))
    return g

cdef double pee(float E, float kT, float Ef, float N0, float sigma, float E0):
    """Number density of carriers."""
    cdef double g
    g = gaussianDoS(E, N0, E0, sigma)
    cdef double p
    p = g * (1 - (1 / (1 + (exp((E - Ef) / (kT))))))
    return p

cpdef double ICTCs(float E, float kT, float Ef, float NA, float sigmaB, float EB):
    """Number density of bound holes, with binding energy EB."""
    cdef double gt
    cdef double fT
    gt = gaussianDoS(E, NA, EB, sigmaB)
    fT = gt * (1 - (1 / (1 + (exp((E - Ef) / (kT))))))
    return fT
    
cpdef double ICTCsSEP(float E, float kT, float Ef, float NA, float sigmaB, float EB):
    """Number density of free electrons, with binding energy EB."""
    cdef double gtSEP
    cdef double fTSEP
    gtSEP = gaussianDoS(E, NA, EB, sigmaB)
    fTSEP = gtSEP * ((1 / (1 + (exp((E - Ef) / (kT))))))
    return fTSEP

cdef double ionisedDopants(float E, float kT, float Ef, float sigma, float EA):
    """Number density of thermally activated dopants."""
    cdef double gNa
    cdef double NAionised
    gNa = gaussianDoS(E, 1, EA, sigma)
    NAionised = gNa * (1 / (1 + exp((E - Ef) / kT)))
    return NAionised
    
cdef double ionisedDopantsSEP(float NaNeg, float kT, float Ef, float EB):
    """Proportion of ionised dopants not bound."""
    cdef double NAionisedSEP
    NAionisedSEP = NaNeg * ( 1 / (1 + exp(((EB - Ef) / kT))))
    return NAionisedSEP


    
def shifting(int iterations, Nrange, float kT, float E0, float sigma, float E1, float EB, float EA):
    """Simulate dopant shifting - that is no relative shifts."""
    dataOut = []

    # Generate data over host blend ratio ZnPc:F8ZnPc, ZnPc content = N
    for N in Nrange:
        NA = np.empty(iterations)
        NaNegSEP = np.ones(iterations)
        NaNeg = np.empty(iterations)
        p = np.empty(iterations)
        q = np.empty(iterations)
        eFermi = np.empty(iterations)
        Efs = np.linspace(-1.5, 1.5, iterations)
        i = 0
        for Ef in Efs:
            # proportion of dopants ionised
            NaNeg1, abserr1 = quad(
                ionisedDopants, -5, 5, (kT, Ef, sigma, EA), points=[EA])
            # number of carriers in ZnPc HOMO
            p1, abserr2 = quad(
                pee, -5, 5, (kT, Ef, N, sigma, E0), points=[E0])
            # number of carriers in F8ZnPc HOMO
            p2, abserr3 = quad(
                pee, -5, 5, (kT, Ef, (1 - N), sigma, E1), points=[E1])

            eFermi[i] = Ef 
            p[i] = p1
            q[i] = p2
            NaNeg[i] = NaNeg1
            i += 1

        NA = ( p + q ) / (NaNeg)
        MR = NA

        # Solved for molar ratio over Fermi level, so to find relevant concentration:
        def find_nearest(array, value):
            array = np.asarray(array)
            idx = (np.abs(array - value)).argmin()
            return idx, array[idx]

        # Format dataframe for saving
        dat = [MR, eFermi, p, q, NaNeg]
        i, x = find_nearest(dat[0], 0.05)
        dat = np.transpose(dat)
        dataOut.append([N, E0, E1, dat[i][0], dat[i][1], dat[i][2], dat[i][3], dat[i][4]])

    return dataOut

    
def fixed(int iterations, Nrange, float kT, float E0, float sigma, float E1, float EB, float EA):
    """Simulate dopant level fixed - only the hosts move in energy."""
    dataOut = []

    # For calculating how much the host energy levels shift
    # Data taken from UPS data in DOI: 10.1126/science.aaf0590
    znpc_ratio = [0.37, 0.67, 0.84, 1.00]
    znpc_energy = [6.12, 5.83, 5.66, 5.55]
    f8znpc_ratio = [0.00, 0.37, 0.67, 0.84]
    f8znpc_energy = [6.73, 6.45, 6.17, 6.04]

    def Linear_fit(x, m, c):
        return (x * m) + c

    popt, pcov = curve_fit(Linear_fit, znpc_ratio, znpc_energy)
    perr = np.sqrt(np.diag(pcov))
    x = np.linspace(0, 1, 20)
    f = Linear_fit(x, popt[0], popt[1])
    
    popt8, pcov8 = curve_fit(Linear_fit, f8znpc_ratio, f8znpc_energy)
    perr8 = np.sqrt(np.diag(pcov8))
    f8 = Linear_fit(x, popt8[0], popt8[1])

    # Generate data over host blend ratio ZnPc:F8ZnPc, ZnPc content = N
    for N in Nrange:
        E0 = -1 * Linear_fit(N, popt[0], popt[1]) + Linear_fit(1, popt[0], popt[1])
        E1 = -1 * Linear_fit(N, popt8[0], popt8[1]) + Linear_fit(1, popt[0], popt[1])
        EBshift = E0 + EB       # Binding level remains an equal distance ZnPc.
        NA = np.empty(iterations)
        NaNegSEP = np.ones(iterations)
        NaNeg = np.empty(iterations)
        p = np.empty(iterations)
        q = np.empty(iterations)
        eFermi = np.empty(iterations)
        Efs = np.linspace(-2, 2.0, iterations)
        i = 0
        for Ef in Efs:
            # proportion of dopants ionised
            NaNeg1, abserr1 = quad(
                ionisedDopants, -4, 4, (kT, Ef, sigma, EA), points=[EA])
            # number of carriers in ZnPc HOMO
            p1, abserr2 = quad(
                pee, -4, 4, (kT, Ef, N, sigma, E0), points=[E0])
            # number of carriers in F8ZnPc HOMO
            p2, abserr3 = quad(
                pee, -4, 4, (kT, Ef, (1 - N), sigma, E1), points=[E1])

            eFermi[i] = Ef 
            p[i] = p1
            q[i] = p2
            NaNeg[i] = NaNeg1
            i += 1

        NA = ( p + q ) / (NaNeg)
        MR = NA

        # Solved for molar ratio over Fermi level, so to find relevant concentration:
        def find_nearest(array, value):
            array = np.asarray(array)
            idx = (np.abs(array - value)).argmin()
            return idx, array[idx]

        # Format dataframe for saving
        dat = [MR, eFermi, p, q, NaNeg]
        i, x = find_nearest(dat[0], 0.05)
        dat = np.transpose(dat)
        dataOut.append([N, E0, E1, dat[i][0], dat[i][1], dat[i][2], dat[i][3], dat[i][4]])

    return dataOut
    

def shiftingEB(int iterations, Nrange, float kT, float E0, float sigma, float E1, float EB, float EA):
    """Simulate dopant shifting - that is no relative shifts PLUS output EB parameters"""
    dataOut = []

    # Generate data over host blend ratio ZnPc:F8ZnPc, ZnPc content = N
    for N in Nrange:
        NA = np.empty(iterations)
        NaNegSEP = np.ones(iterations)
        NaNegICTC = np.ones(iterations)
        NaNeg = np.empty(iterations)
        p = np.empty(iterations)
        q = np.empty(iterations)
        eFermi = np.empty(iterations)
        Efs = np.linspace(-1.5, 1.5, iterations)
        i = 0
        for Ef in Efs:
            # proportion of dopants ionised
            NaNeg1, abserr1 = quad(
                ionisedDopants, -5, 5, (kT, Ef, sigma, EA), points=[EA])
            # number of carriers in ZnPc HOMO
            p1, abserr2 = quad(
                pee, -5, 5, (kT, Ef, N, sigma, E0), points=[E0])
            # number of carriers in F8ZnPc HOMO
            p2, abserr3 = quad(
                pee, -5, 5, (kT, Ef, (1 - N), sigma, E1), points=[E1])
            
            # Of the ionised dopants, this many are BOUND in ICTCs
            NaNegICTC1, abserr4 = quad(ICTCs, -4, 4, (kT, Ef, NaNeg1, sigma, EB), points=[EB])
            # Of the ionised dopants, this many are NOT bound in ICTCs
            NaNegSEP1, abserr3 = quad(ICTCsSEP, -4, 4, (kT, Ef, NaNeg1, sigma, EB), points=[EB])

            eFermi[i] = Ef 
            p[i] = p1
            q[i] = p2
            NaNegICTC[i] = NaNegICTC1
            NaNegSEP[i] = NaNegSEP1
            NaNeg[i] = NaNeg1
            i += 1

        NA = ( p + q ) / (NaNeg)
        MR = NA

        # Solved for molar ratio over Fermi level, so to find relevant concentration:
        def find_nearest(array, value):
            array = np.asarray(array)
            idx = (np.abs(array - value)).argmin()
            return idx, array[idx]

        # Format dataframe for saving
        dat = [MR, eFermi, p, q, NaNeg, NaNegICTC, NaNegSEP]
        i, x = find_nearest(dat[0], 0.05)
        dat = np.transpose(dat)
        dataOut.append([N, E0, E1, dat[i][0], dat[i][1], dat[i][2], dat[i][3], dat[i][4], dat[i][5], dat[i][6]])

    return dataOut