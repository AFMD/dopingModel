"""
Simulation code for figure 4b, DOI:

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


def statModelFixed(iterations, Nrange, kT, E0, sigma, E1, EB, EA):
    '''Calculation for fixed dopant level.'''
    info('Calculation for fixed dopant level')
    start = time.time()
    calcRESULTS = statModel.fixed(iterations, Nrange, kT, E0, sigma, E1, EB, EA)
    end = time.time()
    print('Elapsed time: %.2f seconds' % (end - start))
    return calcRESULTS


if __name__ == '__main__':
    info('main line')

    iterations = 6000           # iterations of Fermi level to solve over
    kT = 0.025                  # thermal energy eV
    E0 = 0                      # ZnPc DOS centre
    E1 = -0.34                  # F8ZnPc DOS centre
    EA = 0.45                   # F6-TCNNQ DOS centre
    EB = 0.5                    # CT binding energy
    Nrange = linspace(0, 1, 51)  # blend ratio interval in MR

    pool = Pool(4)      # Create pool of four workers, same number as CPUs
    # Now create the tasks for workes
    sigma = 0.125
    test1 = [iterations, Nrange, kT, E0, sigma, E1, EB, EA]
    sigma = 0.150
    test2 = [iterations, Nrange, kT, E0, sigma, E1, EB, EA]
    sigma = 0.175
    test3 = [iterations, Nrange, kT, E0, sigma, E1, EB, EA]
    sigma = 0.200
    test4 = [iterations, Nrange, kT, E0, sigma, E1, EB, EA]

    results = [pool.apply_async(statModelFixed, test1),
               pool.apply_async(statModelFixed, test2),
               pool.apply_async(statModelFixed, test3),
               pool.apply_async(statModelFixed, test4)]

    for idx, result in enumerate(results):
    # for result in results:
        x = result.get()
        if idx == 0:
            df = pd.DataFrame(
                x, columns=['N', 'E1', 'E2', 'MR', 'eFermi', 'p', 'q', 'NaNeg'])
            df.to_csv('data/FIXED-DOP-SIGMA-125.csv', sep='\t')
        if idx == 1:
            df = pd.DataFrame(
                x, columns=['N', 'E1', 'E2', 'MR', 'eFermi', 'p', 'q', 'NaNeg'])
            df.to_csv('data/FIXED-DOP-SIGMA-150.csv', sep='\t')

        if idx == 2:
            df = pd.DataFrame(
                x, columns=['N', 'E1', 'E2', 'MR', 'eFermi', 'p', 'q', 'NaNeg'])
            df.to_csv('data/FIXED-DOP-SIGMA-175.csv', sep='\t')
        if idx == 3:
            df = pd.DataFrame(
                x, columns=['N', 'E1', 'E2', 'MR', 'eFermi', 'p', 'q', 'NaNeg'])
            df.to_csv('data/FIXED-DOP-SIGMA-200.csv', sep='\t')

    sys.exit()
