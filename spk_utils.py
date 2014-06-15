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
import pylab as plt

__all__=['entropy','argfindpeaks']

def _calculate_edges_with_fixed_resolution(minedge, maxedge, res=20, 
                                     scale='log', verbose=False):
    ''' Module's internal function to calculate the edges with fixed time bin resolution.
    For the log scale res is the parameter K. The number of bins 
    The output is a np.array vector.
    '''
    res = float(res)
    if scale == 'log':
        K = res*np.log10(maxedge/float(minedge))
        if verbose:
            print('Returning log scale edges (base 10).')
        return float(minedge)*10**(np.arange(1,K+1)/res)
    else:
        if verbose:
            print('Returning linear scale edges.')
        return np.arange(minedge,maxedge,res)

def _calculate_edges_for_entropy_map(minimum=0.4,maximum=250,parameter=20):
    ''' Calculates the edges using trial and error and with fixed limits.
    '''
    edges = [minimum]
    k = 1
    while edges[len(edges)-1]<maximum+1:
        edges.extend([round(minimum*10**(float(k)/parameter),2)]) #adjust parameter....
        k+=1
        if k>1000:
            print "too many edges......k above 1000... restarting with diferent parameter..."
            k=1
            edges=[minimum]
            parameter-=3
            if parameter<0:
                parameter=1
                print "---> We have a problem...."
    return np.array(edges)

def entropy(spiketrain,edges=None, res=20, scale='log',order=2,input_switch = 'spiketrain',plot=False):
    ''' Calculates the entropy of a spiketrain.
    For reference see Dorval et al. 2009.
    Order [2] defines the order to use. 
    '''
    if input_switch in ['spiketrain']:
        isi = np.diff(spiketrain)
    else:
        isi = spiketrain
    if not len(isi)>1:
        # There is no entropy on an empty spiketrain.
        # NOTE: It is zero by definition.
        return np.nan
    else:
        if edges is None:
            # Calculate edges
            # Use limits 0.4 250 in analysing MRG otherwise let the axes be given by the dataset
            edges = _calculate_edges_with_fixed_resolution(0.4,250,res,scale)
                                                           #np.min(isi)*0.8, 
                                                           #np.max(isi)*1.2,
                                                           #res,scale)
#            edges = _calculate_edges_for_entropy_map(0.4,250,parameter=20)
        # In Dorval et al. 2009 the probabity is defined 
        # as the # of isis in each bin divided by the number of events,
        # so we can not use the density.
        if order == 1:
            counts,edges = np.histogram(isi,edges,density=False)
            Pisi = counts/float(len(isi))
            idx = np.nonzero(Pisi)
            return -np.sum(Pisi[idx]*np.log2(Pisi[idx]))
        if order == 2:
            counts,edges0,edges1 = np.histogram2d(isi[:-1],isi[1:],edges)
            Pisi = counts/float(len(isi)-1)
            # idx = np.nonzero(Pisi)

            x,y = np.where(Pisi)
            if plot:
                ax =  plt.gca()
                linedge = edges.copy()
                for ii in range(0,len(edges0)-1):
                    edges0[ii] = edges0[ii]+(edges0[ii+1]-edges0[ii])/2.0 
                    edges1[ii] = edges1[ii]+(edges1[ii+1]-edges1[ii])/2.0
                ax.scatter(np.log(edges0[x]),np.log(edges1[y]),c=Pisi[x,y], 
                           s=counts[x,y] / np.sum(counts[x, y]),
                           cmap=plt.cm.gray, vmax=1, vmin=0)
                t = [1, 2, 5, 10, 50, 100, 200]
                ax.set_yticks(np.log(t), t)
                ax.set_xticks(np.log(t), t)
                ax.axis([0, np.log(300), 0, np.log(300)])
                ax.grid(1, alpha = 0.6)
#                ax.scatter(isi[:-1],isi[1:],20,'r')
                # ax.hlines(linedge,min(edges[::3]),max(edges[::3] ),lw=0.3)
                # ax.vlines(linedge,min(edges[:: 3]),max(edges[:: 3]),lw=0.3)
                return (edges, counts)
            return -0.5*np.sum(Pisi[x,y]*np.log2(Pisi[x,y]))

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
