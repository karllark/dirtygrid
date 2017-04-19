#!/usr/bin/env python

from __future__ import (absolute_import, print_function, division)

#import tables as pt
#import numpy as np
#from astropy.io import fits
import matplotlib.pyplot as plt

import PhotDG

if __name__ == '__main__':
    # get the DG
    dg = PhotDG.PhotDG()

    # print the DG parameters
    dg.print_parameters()

    # get the photometry for a grid point
    dg.photGet('RV31_BC6-WD01', 'dusty', 'burst', 0.1, 10.0, 1e6, 1.0)

    # plot that SED
    dg.photPlot()
