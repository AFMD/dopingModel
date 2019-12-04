"""
Simulation code for figure 4c, DOI:

Please see readme on how to run. Tested only on Arch Linux.

Author: Ross <ross dot warren at pm dot me>
Date:   04/12/2019
"""

import sys
import time
import datetime
import os
import pandas as pd
from numpy import linspace
from multiprocessing import Pool

import statModel

date = str(datetime.date.today())

def info(title):
    print(title)
    print('module name:', __name__)
    if hasattr(os, 'getppid'):
        print('parent process:', os.getppid())
    print('process id:', os.getpid())


def statModelShifting(iterations, Nrange, kT, E0, sigma, E1, EB, EA):
    '''Calculation for all levels shifting (relative distances fixed).'''
    info('Calculation for all levels shifting')
    # Run Program
    start = time.time()
    calcRESULTS = statModel.shifting(iterations, Nrange, kT, E0, sigma, E1, EB, EA)
    end = time.time()
    print('Elapsed time: %.2f seconds' % (end - start))
    return calcRESULTS


if __name__ == '__main__':
    info('main line')

    iterations = 6000           # iterations of Fermi level to solve over
    kT = 0.025                  # thermal energy eV
    E0 = 0                      # ZnPc DOS centre
    E1 = -0.34                  # F8ZnPc DOS centre
    sigma = 0.160               # Standard dev of all Gaussians
    EB = 0.5                    # CT binding energy
    Nrange = linspace(0, 1, 51)  # blend ratio interval in MR
    
    pool = Pool(4)      # Create pool of four workers, same number as CPUs
    # Now create the tasks for workes
    EA = 0.40
    test1 = [iterations, Nrange, kT, E0, sigma, E1, EB, EA]
    EA = 0.45
    test2 = [iterations, Nrange, kT, E0, sigma, E1, EB, EA]
    EA = 0.50
    test3 = [iterations, Nrange, kT, E0, sigma, E1, EB, EA]
    EA = 0.55
    test4 = [iterations, Nrange, kT, E0, sigma, E1, EB, EA]

    results = [pool.apply_async(statModelShifting, test1),
               pool.apply_async(statModelShifting, test2),
               pool.apply_async(statModelShifting, test3),
               pool.apply_async(statModelShifting, test4)]

    # For saving the simulation data
    for idx, result in enumerate(results):
        x = result.get()
        if idx == 0:
            df = pd.DataFrame(
                x, columns=['N', 'E1', 'E2', 'MR', 'eFermi', 'p', 'q', 'NaNeg'])
            df.to_csv('data/SHIFTING-DOP-EA-0p40.csv', sep='\t')
        if idx == 1:
            df = pd.DataFrame(
                x, columns=['N', 'E1', 'E2', 'MR', 'eFermi', 'p', 'q', 'NaNeg'])
            df.to_csv('data/SHIFTING-DOP-EA-0p45.csv', sep='\t')
        if idx == 2:
            df = pd.DataFrame(
                x, columns=['N', 'E1', 'E2', 'MR', 'eFermi', 'p', 'q', 'NaNeg'])
            df.to_csv('data/SHIFTING-DOP-EA-0p50.csv', sep='\t')
        if idx == 3:
            df = pd.DataFrame(
                x, columns=['N', 'E1', 'E2', 'MR', 'eFermi', 'p', 'q', 'NaNeg'])
            df.to_csv('data/SHIFTING-DOP-EA-0p55.csv', sep='\t')

    sys.exit()
