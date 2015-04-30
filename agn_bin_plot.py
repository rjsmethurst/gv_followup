import numpy as N
import pylab as P
import scipy as S
import pyfits as F

ags = F.open('/Users/becky/Projects/Green-Valley-Project/data/GZ2_sample_match_with_schawinski_agn_catalogue.fits')
a = ags[1].data

bx = N.linspace(0, 14, 10)
by = N.linspace(0, 4, 10)

agn = a[N.where(a['bpt_class']==3)]
not_agn = a[N.where(a['bpt_class']< 2)]

idx = N.random.randint(0, len(not_agn), len(agn))
non_agn = not_agn[idx]


bdx = N.digitize(non_agn['best_t'], bx)
bdy = N.digitize(non_agn['best_tau'], by)

bdxagn = N.digitize(agn['best_t'], bx)
bdyagn = N.digitize(agn['best_tau'], by)


def bin_calc(array, bdx, bdy):
	grid = N.zeros((len(by), len(bx)))
	array_grid = N.zeros((len(by), len(bx)))
	for n in range(len(array)):
		grid[bdy[n]-1, bdx[n]-1] += 1
		array_grid[bdy[n]-1, bdx[n]-1] += array[n]
	avg_array = array_grid/grid
	return array_grid, avg_array, grid


def bin_plot(a, c, b, array):
 	fig = P.figure(figsize=(12,4))
 	ax1 = P.subplot(131)
 	ap1 = ax1.imshow(N.log10(a), origin='lower', aspect='auto', vmin=0, vmax=N.max(N.log10(a)), extent=(0, 13.8, 0, 4), interpolation='nearest', cmap=P.cm.gist_heat_r)
 	ax1.set_xlabel(r'$t_q$ [Gyr]')
 	ax1.set_ylabel(r'$\tau$ [Gyr]')
 	cbar1 = P.colorbar(ap1)
 	cbar1.set_label(r'log count of galaxies in each bin', fontsize=8)
 	ax2 = P.subplot(132)
 	ap2 = ax2.imshow(N.log10(b), origin='lower', aspect='auto', vmin=0, vmax=N.max(N.log10(a)), extent=(0, 13.8, 0, 4), interpolation='nearest', cmap=P.cm.gist_heat_r)
 	ax2.set_xlabel(r'$t_q$ [Gyr]')
 	ax2.set_ylabel(r'$\tau$ [Gyr]')
 	cbar2 = P.colorbar(ap2)
 	cbar2.set_label(r'log count of '+array+' galaxies in each bin', fontsize=8)
 	ax3 = P.subplot(133)
 	ap3 = ax3.imshow(c, origin='lower', aspect='auto', extent=(0, 13.8, 0, 4), interpolation='nearest', cmap=P.cm.gist_heat_r)
 	ax3.set_xlabel(r'$t_q$ [Gyr]')
 	ax3.set_ylabel(r'$\tau$ [Gyr]')
 	cbar3 = P.colorbar(ap3)
 	cbar3.set_label(r'fractional count of '+array+' galaxies in each bin', fontsize=8)
 	P.tight_layout()
 	P.savefig('fractional_'+array+'_new.png')

alpha = 8.13
beta = 4.02
sigma_0 = 200
sigma_agn = N.nan_to_num(agn['sigma'])
sigma_non_agn = N.nan_to_num(non_agn['sigma'])
log_mbh_agn = alpha + beta*(N.log10(sigma_agn/sigma_0))
log_mbh_non_agn = alpha + beta*(N.log10(sigma_non_agn/sigma_0))

total_mbh_agn, avg_mbh_agn, count_agn = bin_calc(log_mbh_agn, bdxagn, bdyagn)
total_mbh_non_agn, avg_mbh_non_agn, count_non_agn = bin_calc(log_mbh_non_agn, bdx, bdy)

P.figure(figsize=(8,4))
ax3 = P.subplot(121)
ap3 = ax3.imshow(avg_mbh_agn, origin='lower', aspect='auto', extent=(0, 13.8, 0, 4), interpolation='nearest', vmin=5.0, vmax=8.8)
ax3.set_xlabel(r'$t_q$ [Gyr]')
ax3.set_ylabel(r'$\tau$ [Gyr]')
cbar3 = P.colorbar(ap3)
cbar3.set_label(r'average $log [M_{BH}/M_{\odot}]$ of AGN galaxies', fontsize=8)
ax3 = P.subplot(122)
ap3 = ax3.imshow(avg_mbh_non_agn, origin='lower', aspect='auto', extent=(0, 13.8, 0, 4), interpolation='nearest', vmin=5.0, vmax=8.8)
ax3.set_xlabel(r'$t_q$ [Gyr]')
ax3.set_ylabel(r'$\tau$ [Gyr]')
cbar3 = P.colorbar(ap3)
cbar3.set_label(r'average $log [M_{BH}/M_{\odot}]$ of non AGN galaxies', fontsize=8)
P.tight_layout()
P.savefig('average_log_bh_mass_binned_agn_non_agn.png')

P.figure(figsize=(8,4))
ax3 = P.subplot(121)
ap3 = ax3.imshow(N.log10(count_agn), origin='lower', aspect='auto', extent=(0, 13.8, 0, 4), interpolation='nearest')
ax3.set_xlabel(r'$t_q$ [Gyr]')
ax3.set_ylabel(r'$\tau$ [Gyr]')
cbar3 = P.colorbar(ap3)
cbar3.set_label(r'count of AGN galaxies', fontsize=8)
ax3 = P.subplot(122)
ap3 = ax3.imshow(N.log10(count_non_agn), origin='lower', aspect='auto', extent=(0, 13.8, 0, 4), interpolation='nearest')
ax3.set_xlabel(r'$t_q$ [Gyr]')
ax3.set_ylabel(r'$\tau$ [Gyr]')
cbar3 = P.colorbar(ap3)
cbar3.set_label(r'count of non AGN galaxies', fontsize=8)
P.tight_layout()
P.savefig('average_count_binned_agn_non_agn.png')


P.figure(figsize=(8,4))
ax3 = P.subplot(121)
ap3 = ax3.imshow(count_agn/count_non_agn, origin='lower', aspect='auto', extent=(0, 13.8, 0, 4), interpolation='nearest')
ax3.set_xlabel(r'$t_q$ [Gyr]')
ax3.set_ylabel(r'$\tau$ [Gyr]')
cbar3 = P.colorbar(ap3)
cbar3.set_label(r'fraction of AGN galaxies', fontsize=8)
ax3 = P.subplot(122)
ap3 = ax3.imshow(N.log10(count_non_agn), origin='lower', aspect='auto', extent=(0, 13.8, 0, 4), interpolation='nearest')
ax3.set_xlabel(r'$t_q$ [Gyr]')
ax3.set_ylabel(r'$\tau$ [Gyr]')
cbar3 = P.colorbar(ap3)
cbar3.set_label(r'count of non AGN galaxies', fontsize=8)
P.tight_layout()
P.savefig('average_fractional_count_binned_agn_non_agn.png')

total_lum_agn, avg_lum_agn, count_agn = bin_calc(N.log10(agn['L_03']), bdxagn, bdyagn)
total_lum_non_agn, avg_lum_non_agn, count_non_agn = bin_calc(N.log10(non_agn['L_03']), bdx, bdy)

P.figure(figsize=(8,4))
ax3 = P.subplot(121)
ap3 = ax3.imshow(avg_lum_agn, origin='lower', aspect='auto', extent=(0, 13.8, 0, 4), interpolation='nearest')
ax3.set_xlabel(r'$t_q$ [Gyr]')
ax3.set_ylabel(r'$\tau$ [Gyr]')
cbar3 = P.colorbar(ap3)
cbar3.set_label(r'average $log [L_[OIII] erg/s]$ of AGN galaxies', fontsize=8)
ax3 = P.subplot(122)
ap3 = ax3.imshow(avg_lum_non_agn, origin='lower', aspect='auto', extent=(0, 13.8, 0, 4), interpolation='nearest')
ax3.set_xlabel(r'$t_q$ [Gyr]')
ax3.set_ylabel(r'$\tau$ [Gyr]')
cbar3 = P.colorbar(ap3)
cbar3.set_label(r'average $log [L_[OIII] erg/s]$ of non AGN galaxies', fontsize=8)
P.tight_layout()
P.savefig('average_luminosity_binned_agn_non_agn.png')



