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

def _calculate_edges_with_fixed_resolution(minedge, maxedge, res=20, 
                                     scale='log', verbose=False):
    ''' Module's internal function to calculate the edges with fixed time bin resolution.
    For the log scale res is the parameter K. The number of bins 
    The output is a np.array vector.
    '''
    res = float(res)
    if scale == 'log':
        K = res*np.log10(maxedge/float(minedge))
        print K
        
        if verbose:
            print('Returning log scale edges (base 10).')
        print np.arange(1,K+1)
        return float(minedge)*10**(np.arange(1,K+1)/res)
    else:
        if verbose:
            print('Returning linear scale edges.')
        return np.arange(minedge,maxedge,res)
    
def entropy(spiketrain, edges=None, res=20, scale='log'):
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
        if edges is None:
            # Calculate edges
            edges = _calculate_edges_with_fixed_resolution(np.min(isi)*0.8, 
                                                     np.max(isi)*1.2,
                                                     res,scale)
        # In Dorval et al. 2009 the probabity is defined 
        # as the # of isis in each bin divided by the number of events,
        # so we can not use the density.
        counts,edges = np.histogram(isi,edges,density=False)
        Pisi = counts/float(len(isi))
        idx = np.nonzero(Pisi)
        H = -np.sum(Pisi[idx]*np.log2(Pisi[idx]))
        return H
    
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
