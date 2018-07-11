#!/usr/bin/env python

import csv
import numpy as np
from scipy import stats
import os
import re
from matplotlib import pyplot as plt
from astropy import units as u
from astropy import constants as const
from astropy.table import Table, Column
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# Radial profile for one galaxy
def pltcoprof(gal, ax, rtype='arcsec', ebar=True):
    print('Plotting {0}...'.format(gal))
    b0dat = Table.read('rprof/'+gal+'.b0.rprof.txt', format='ascii.ecsv')
    b1dat = Table.read('rprof/'+gal+'.b1.rprof.txt', format='ascii.ecsv')
    b2dat = Table.read('rprof/'+gal+'.b2.rprof.txt', format='ascii.ecsv')
    if rtype == 'kpc':
        rad = b0dat['r_pc']/1000.
        rmax = 15.
    elif rtype == 'r25':
        rad = b0dat['r_r25']
        rmax = 1.
    else:
        rad = b0dat['radius']
        rmax = 90.
    plt.plot(rad, b0dat['wtmean'], color='b', ls='-', marker='^', ms=7)
    plt.plot(rad, b1dat['wtmean'], color='g', ls='--', marker=None)
    #plt.plot(rad,b2dat['wtmean'],color='r',ls=':',marker=None)
    plt.fill_between(rad,b0dat['wtmean'],
        b1dat['wtmean'], facecolor='green', alpha=0.5)
    if ebar == True:
        plt.errorbar(rad, b0dat['wtmean'], yerr=b0dat['rms'], 
            ecolor='dimgray', capsize=0, 
            zorder=1, marker=None, ls='None', lw=1, label=None)
    rcoeff = np.interp(0.5*b0dat['cummass'][-1], b0dat['cummass'], rad)
    plt.plot(rad, b0dat['detlim'], color='k', lw=2, ls='--', marker=None)
    ax.axvline(x=rcoeff, lw=4, color='r', alpha=0.5)
    ax.set_yscale('log')
    ax.set_ylim(10**-1, 10**3)
    ax.set_xlim(0,rmax)
    ax.set_title(gal, fontsize='large')
#     Do a linear regression fit if requested
#     if linfit is not None:
#         goodidx = (xdata>0) & (ydata>0)
#         xdata = xdata[goodidx]
#         ydata = ydata[goodidx]
#         sorted=np.argsort(xdata)
#         m, b, rval, pval, std_err = stats.linregress(np.log10(xdata[sorted]),
#             np.log10(ydata[sorted]))
#         xmod = np.logspace(-3,6,100)
#         ymod = b+m*np.log10(xmod)
#         axes.plot(xmod, 10**(ymod), linestyle='--', color=linfit)
#         axes.text(0.03,0.95,'slope = $%4.2f$' % m,size=10,transform=axes.transAxes)
#         axes.text(0.03,0.90,'intercept = $%4.2f$' % b,size=10,transform=axes.transAxes)
    return

# Make the plot
dir = '/Users/tonywong/Work/projects/EDGE/edge-sql-base/global_values/external/'
data = Table.read(dir+'edge_leda.csv', format='ascii.ecsv')
data.sort('ledaRA')
list=data['Name'].tolist()
nx=7
ny=5
pages = int(np.ceil(float(len(list)) / (nx*ny)))

for num in range(0,pages):
    aa = nx*ny*num
    bb = nx*ny+aa
    gals = list[aa:bb]
    nrows = len(gals)//nx
    #fig = plt.figure(figsize=(20,14))
    figure = plt.figure(0)
    figure.set_size_inches(nx*4.5, nrows*4.)

    for i, gal in enumerate(gals):
        row,col = divmod(i,nx)
        ax = plt.subplot2grid((nrows,nx),(row,col))
        if (gal not in ['UGC05498NED01','NGC4211NED02','NGC6394']):
            pltcoprof(gal, ax, rtype='r25', ebar=False)
    if num == 3:
        figure.text(0.5, 0.02, 'Radius / '+r'$R_{25}$', ha='center', fontsize=24)
    else:
        figure.text(0.5, 0.06, 'Radius / '+r'$R_{25}$', ha='center', fontsize=24)
    figure.text(0.09, 0.5, 'log I [K km/s]', va='center', 
        rotation='vertical', fontsize=24)

    plt.savefig('coprof-'+str(num)+'.pdf', bbox_inches='tight')
    plt.close()

os.system("pdfunite coprof-*.pdf edge_coprof.pdf")
#os.system("rm -f "+set+'_'+type+"_flux-*.pdf")
