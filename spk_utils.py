#!/usr/bin/env python
'''
Spike utilities.
This module contains functions for:
    - spiketrain statistics
    - spike detection
Joao Couto - 17 Jan 2013
'''

import numpy as np
from scipy.signal import find_peaks_cwt

__all__=['entropy']

def _calculate_edges_with_fixed_nbin(minimum, maximum, nbins=20, 
                                     scale='log', verbose=False):
    ''' Module's internal function to calculate the edges with fixed bin number.
    The 
    The output is a np.array vector.
    '''
    if scale == 'log':
        if verbose:
            print('Returning logarithmic scale edges.')
        return min_edge*(10**(np.arange(1,nbins+1)/float(nbins)))
    else:
        if verbose:
            print('Returning linear scale edges.')
        return np.linspace(min_edge,max_edge,nbins)
    
def entropy(spiketrain, nbins=20, scale='log'):
    ''' Calculates the entropy of a spiketrain.
    For reference see Dorval et al. 2009.
    Order [2] defines the order to use. 
    '''
    isi = np.diff(spiketrain)
    if not len(isi):
        # There is no entropy on an empty spiketrain.
        # NOTE: It is zero by definition.
        return np.nan
    else:
        # Calculate edges
        edges = _calculate_edges_with_fixed_nbin(np.min(isi)*0.8, 
                                                 np.max(isi)*1.2,
                                                 nbins,scale)
        
        return 1
    
def argfindpeaks(data, threshold=-20,deadwindow=30):
    ''' Extracts the indexes of the peaks in the data with threshold crossing.
    Uses a dead window/minimum peak distance.
    '''
    N = len(data)
    ii = np.arange(0,N)
    # threshold crossing
    dx = data > threshold
    idx=ii[dx][np.diff(ii[dx])>1]
    idx = np.append(idx,[ii[dx][-1]])
    # find peaks using the dead window
    index = []
    for ii in idx:
        lower = ii - deadwindow
        upper = ii + deadwindow
        if  lower < 0:
            lower = 0
        if upper > N:
            upper = N
        index.append(lower + np.argmax(data[lower:upper]))
    return np.array(index)
