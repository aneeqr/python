#!/usr/bin/env python
'''
Plotting Utilities
This module contains utilities for:
    - spiketrain visualization
    - quick axes manipulation
    - automatic caption generation and printing.
Requires matplotlib.pylab
'''
# fix the ticks direction
from matplotlib import rcParams
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'

import os
import matplotlib.pylab as plt
import numpy as np
# For printFigureWithCaption
import subprocess
import shlex
import shutil

__all__ = ['plotRastergram','printFigureWithCaption','fixAxes']

def plotRastergram(spiketrains,offset=0,value=None,color=[0.1,0.1,0.1]):
    '''
    Plots a rastergram of the spike trains
    Spiketrains should be a list of spiketrains
    '''
    if value is None:
        for i,sp in enumerate(spiketrains):
            X = np.tile(sp,(2,1))
            Y = np.transpose(np.tile(np.transpose([i,i+0.8]),(len(sp),1))+offset)
            plt.plot(X,Y,color=color,lw=0.6)
    else:
        mindist = np.min(np.diff(value))
        for i,sp in enumerate(spiketrains):
            X = np.tile(sp,(2,1))
            Y = np.transpose(np.tile(np.transpose([value[i],value[i]+mindist]),(len(sp),1))+offset)
            plt.plot(X,Y,color=color,lw=0.7)

def printFigureWithCaption(filename,caption,fig=None,savefig_args={}, hide_output=True):
    '''Uses pylab.savefig to save the file and appends a caption using latex.
    Requires a working installation of pdflatex.
    WARNING: Files are overwritten without prompt.
    '''
    if fig is None:
        fig = plt.gcf()
    if not 'format' in savefig_args.keys():
        savefig_args['format']='pdf'
    TMPFIGNAME = 'tmpFig.'+savefig_args['format']
    TMPLATEXNAME = 'tmpLatex'
    latex_code = r'''\documentclass{article}
\usepackage{amsmath}
\usepackage[active,tightpage,textmath,displaymath,floats,graphics,previewborder=0.05cm]{preview}
\usepackage{graphicx}
\begin{document}
    \begin{figure}
        \centering
        \includegraphics[width=\textwidth]{FILENAME}
        \caption{CAPTION}
    \end{figure}
\end{document}'''
    latex_code = latex_code.replace('FILENAME',TMPFIGNAME)
    latex_code = latex_code.replace('CAPTION',caption)
    with open(TMPLATEXNAME+'.tex','w') as f:
        f.write(latex_code)
        f.close()
    fig.savefig(TMPFIGNAME, **savefig_args)
    # plt.show() # for debuging
    # Use latexpdf to append caption
    if os.path.isfile(TMPFIGNAME):
        stdout = None
        if hide_output:
            stdout = open('/dev/null', 'w')
        proc = subprocess.Popen(shlex.split('pdflatex -interaction=nonstopmode '
                                            +TMPLATEXNAME+'.tex'),stdout=stdout)
        proc.wait()
        # Remove unwanted files and move to requested location
        os.unlink(TMPLATEXNAME+'.tex')   
        os.unlink(TMPLATEXNAME+'.aux')   
        os.unlink(TMPLATEXNAME+'.log')
        shutil.move(TMPLATEXNAME+'.pdf',filename)
        os.unlink(TMPFIGNAME)   
    else:
        print('Figure not saved...HELP...')
        plt.show()
def fixAxes(ax=None,xlabel='',ylabel='',xloc='bottom',yloc='left',xposition=('outward',1.5),yposition=('outward',1.5),lw=0.7,fontsize=8,xcolor='black',ycolor='black'):
    notxloc='bottom'
    notyloc='right'
    if xloc=='bottom':
        notxloc='top'
    if yloc=='right':
        notyloc='left'
    if ax==None:
        ax=plt.gca()
    ax.set_axisbelow(True)
    ax.spines[notxloc].set_visible(False)
    ax.spines[notyloc].set_visible(False)
    ax.spines[xloc].set_visible(True)
    ax.spines[yloc].set_visible(True)
    ax.spines[xloc].set_lw(lw)
    ax.spines[xloc].set_color(xcolor)
    ax.spines[xloc].set_position(xposition)
    ax.spines[yloc].set_lw(lw)
    ax.spines[yloc].set_color(ycolor)
    ax.spines[yloc].set_position(yposition)
    ax.yaxis.set_ticks_position(yloc)
    ax.xaxis.set_ticks_position(xloc)
    xax = ax.xaxis
    xax.set_label_text(xlabel)
    xax.set_label_position(xloc)
    xax.set_label_position(xloc)
    xax.label.set_color(xcolor)
    xax.label.set_size(fontsize)
    yax = ax.yaxis
    yax.set_label_text(ylabel)
    yax.set_label_position(yloc)
    yax.set_label_position(yloc)
    yax.label.set_color(ycolor)
    yax.label.set_size(fontsize)
    for label in ax.get_xticklabels():
            label.set_color(xcolor)
            label.set_fontsize(fontsize)
    for label in ax.get_yticklabels():
            label.set_color(ycolor)
            label.set_fontsize(fontsize)
    for label in ax.get_xticklines():
            label.set_color(xcolor)
            label.set_lw(lw)
            label.set_markersize(4)
    for label in ax.get_yticklines():
            label.set_color(ycolor)
            label.set_lw(lw)
            label.set_markersize(4)
