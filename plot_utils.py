#!/usr/bin/env python
'''
Plotting Utilities
This module contains utilities for:
    - spiketrain visualization
    - quick axes manipulation
    - automatic caption generation and printing.
Requires matplotlib.pylab
'''

import matplotlib.pylab as plt
import numpy as np

__all__ = ['plotRastergram']

def plotRastergram(spiketrains,offset=0,color=[0.1,0.1,0.1]):
    '''
    Plots a rastergram of the spike trains
    Spiketrains should be a list of spiketrains
    '''
    for i,sp in enumerate(spiketrains):
        X = np.tile(sp,(2,1))
        Y = np.transpose(np.tile(np.transpose([i,i+1]),(len(sp),1))+offset)
        plt.plot(X,Y,color=color,lw=1)
