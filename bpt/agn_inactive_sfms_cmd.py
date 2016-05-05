
import numpy as N
import matplotlib.pyplot as P

import os

from astropy.cosmology import FlatLambdaCDM

from astropy import units as un
from astropy.table import Table, vstack

import math

os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'


P.rc('figure', facecolor='none', edgecolor='none', autolayout=True)
P.rc('path', simplify=True)
P.rc('text', usetex=True)
P.rc('font', family='serif')
P.rc('axes', labelsize='large', facecolor='none', linewidth=0.7, color_cycle = ['k', 'r', 'g', 'b', 'c', 'm', 'y'])
P.rc('xtick', labelsize='medium')
P.rc('ytick', labelsize='medium')
P.rc('lines', markersize=4, linewidth=1, markeredgewidth=0.2)
P.rc('legend', numpoints=1, frameon=False, handletextpad=0.3, scatterpoints=1, handlelength=2, handleheight=0.1)
P.rc('savefig', facecolor='none', edgecolor='none', frameon='False')

params =   {'font.size' : 11,
            'xtick.major.size': 8,
            'ytick.major.size': 8,
            'xtick.minor.size': 3,
            'ytick.minor.size': 3,
            }
P.rcParams.update(params) 

agn = Table.read('bpt_identified_type2_agn_seyferts_no_liners_gz2_galex_matched_h_alpha_eqw_sigma_eddington_lum_ratios_stellar_mass_extra_axis_ratio_ab_MASS_SFR.fit', format='fits')
inac = Table.read('inactive_sample_matched_seyferts_mass_pd_ps_pm_5pc_MPA_JHU_MASS_SFR.fit', format='fits')
sdss = Table.read('MPA_JHU_MASS_SFR.fit', format='fits')
alls = Table.read('/Users/becky/Projects/Green-Valley-Project/data/GalaxyZoo1_DR_extra.fits', format='fits')

from plothist2d import plothist2d
import triangle
Hi, Xi, Yi = N.histogram2d(inac['AVG_MASS'], inac['AVG_SFR'], bins=30, range=((8,12.5),(-2,2)), normed=True)
Hic, Xic, Yic = N.histogram2d(alls['MR'], alls['MU']-alls['MR'], bins=30, range=((-24,-18),(0.5,3.5)), normed=True)

Hs, Xs, Ys = N.histogram2d(sdss['AVG_MASS'], sdss['AVG_SFR'], bins=50, range=((8,12.5),(-2,2)), normed=True)

Z = N.genfromtxt('/Users/becky/Projects/Green-Valley-Project/bayesian/find_t_tau/colour_mag_diag/baldry_hist_data.txt')
X = N.genfromtxt('/Users/becky/Projects/Green-Valley-Project/bayesian/find_t_tau/colour_mag_diag/baldry_x_hist.txt')
Y = N.genfromtxt('/Users/becky/Projects/Green-Valley-Project/bayesian/find_t_tau/colour_mag_diag/baldry_y_hist.txt')

def peng_sfr(m,t):
    return (2.5*((m/10**10)**(-0.1))*((t/3.5)**(-2.2)))*(ms/1E9)

ms = N.array((1E8, 1E12, 51))
sfr_138 = N.log10(peng_sfr(ms, 13.8))


P.figure(figsize=(3,3))
ax1 = P.subplot(111)
ax1.pcolor(Xs[:-1], Ys[:-1], Hs.T, cmap=P.cm.binary, alpha=0.8)
ax1.contour(Xs[:-1], Ys[:-1], Hs.T, 10, colors='k', alpha=0.8)
#ax1.contour(Xi[:-1], Yi[:-1], Hi.T, colors='b')
ax1.scatter(agn['MASS'], N.log10(agn['SFR']), marker='o', color='r', alpha=0.4)
ax1.plot(N.log10(ms), sfr_138, linestyle='dashed', color='k')
ax1.set_xlabel(r'$\log_{10}[M_*/M_{\odot}]$')
ax1.set_ylabel(r'$\log_{10} SFR_{H\alpha} [M_{\odot} yr^{-1}]$')
ax1.set_xlim(9.1, 11.9)
ax1.set_ylim(-2,1.5)
ax1.minorticks_on()
#ax1.set_yscale('log')
P.tight_layout()
P.savefig('mass_sfr_agn_sdss_inac.pdf')

Mr = N.linspace(-24.5, -17, 200)
C_dash = 2.06 - 0.244*N.tanh((Mr + 20.07)/1.09)
upper = C_dash + 0.128
lower = C_dash - 0.128

P.figure(figsize=(3,3))
ax1 = P.subplot(111)
ax1.plot(Mr, C_dash, color='k')
ax1.plot(Mr, upper, color='k', linestyle='dashed')
ax1.plot(Mr, lower, color='k', linestyle='dashed')
#ax1.contour(Xic[:-1], Yic[:-1], Hic.T, colors='k', linestyle='dashed', alpha=0.7)
#ax1.pcolor(Xic[:-1], Yic[:-1], Hic.T, cmap=P.cm.binary, linestyle='dashed', alpha=0.3)
plothist2d(Z.T, X, Y, ax=ax1, plot_contours=True, plot_datapoints=False, levels=[80, 160, 320, 640, 1280, 1920, 2560, 5120])
ax1.scatter(agn['MR'], agn['MU_MR'], marker='o', color='r', alpha=0.4)
ax1.set_xlabel(r'$M_r$')
ax1.set_ylabel(r'$u-r$')
ax1.set_xlim(-17.5, -23.5)
ax1.set_ylim(0.5,3.5)
ax1.minorticks_on()
P.tight_layout()
P.savefig('cmd_agn_sdss_inac.pdf')


