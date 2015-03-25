import numpy as N
import pylab as P
import scipy as S
import pyfits as F

d = F.open('/Users/becky/Projects/Green-Valley-Project/data/GZ2_all_GALEX_match_GZ1_k_correct_MPA_JHU_SFR_mass_with_best_fit')
data = d[1].data
#x = da[N.where(da['AVG_MASS'] < 12.0)]
#data = x[N.where(x['AVG_MASS'] > 11.5)]

ags = F.open('/Users/becky/Projects/Green-Valley-Project/data/GZ2_sample_match_with_schawinski_agn_catalogue.fits')
agns = ags[1].data

bx = N.linspace(0, 14, 10)
by = N.linspace(0, 4, 10)

bdx = N.digitize(data['best_t'], bx)
bdy = N.digitize(data['best_tau'], by)

bdxagn = N.digitize(agns['best_t'], bx)
bdyagn = N.digitize(agns['best_tau'], by)


def bin_calc(array, bdx, bdy):
	grid = N.zeros((len(by), len(bx)))
	array_grid = N.zeros((len(by), len(bx)))
	for n in range(len(array)):
		grid[bdy[n]-1, bdx[n]-1] += 1
		array_grid[bdy[n]-1, bdx[n]-1] += array[n]
	#grid[N.where(grid<100)] = 0
	#array_grid[N.where(array_grid<10)] = 0
	grid_array = array_grid/grid
	return array_grid, grid_array, grid

def bin_avg(values, bdx, bdy):
	cum_grid = N.zeros((len(by), len(bx)))
	count_grid = N.zeros((len(by), len(bx)))
	for n in range(len(values)):
		cum_grid[bdy[n]-1, bdx[n]-1] += values[n]
		count_grid[bdy[n]-1, bdx[n]-1] += 1
	return cum_grid / count_grid




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

mass = data['avg_mass']
avg_mass_bin = bin_avg(mass, bdx, bdy)
P.figure()
ax3 = P.subplot(111)
ap3 = ax3.imshow(avg_mass_bin, origin='lower', aspect='auto', extent=(0, 13.8, 0, 4), interpolation='nearest')
ax3.set_xlabel(r'$t_q$ [Gyr]')
ax3.set_ylabel(r'$\tau$ [Gyr]')
cbar3 = P.colorbar(ap3)
cbar3.set_label(r'average stellar mass $[log M_{\odot}]$ in each bin', fontsize=8)
P.tight_layout()
P.savefig('average_mass_binned_all.png')

colour = ['red', 'green', 'blue']
mass = [9.0, 9.5, 10.0, 10.5, 11.0, 11.5]
morph = ['t01_smooth_or_features_a01_smooth_flag', 't01_smooth_or_features_a02_features_or_disk_flag']
morphl = ['smooth', 'disc']
for n in range(len(colour)):
	j = 0
	print j
	P.figure(figsize=(7, 15))
	for m in range(len(mass)):
		if m == (len(mass)-1):
			pass
		else:
			for k in range(len(morph)):
				j+=1
				try:
					bdxagns = N.digitize(agns['best_t'], bx)
					bdyagns = N.digitize(agns['best_tau'], by)
					bdxagni = N.digitize(agns['best_t'][N.where((agns[colour[n]]==1)& (agns[morph[k]]==1)& (agns['log_mstellar']>mass[m])& (agns['log_mstellar']<mass[m+1]))], bx)
					bdyagni = N.digitize(agns['best_tau'][N.where((agns[colour[n]]==1)& (agns[morph[k]]==1)& (agns['log_mstellar']>mass[m])& (agns['log_mstellar']<mass[m+1]))], by)
				except ValueError:
					continue
				agn = agns['agn?']
				agn[N.where(agn < 0)]=0
				agni = agn[N.where((agns[colour[n]]==1)& (agns[morph[k]]==1)& (agns['log_mstellar']>mass[m])& (agns['log_mstellar']<mass[m+1]))]
				a, frac, c = bin_calc(agn, bdxagns, bdyagns)
				b, frac, d = bin_calc(agni, bdxagni, bdyagni)
				ax3 = P.subplot(len(mass)-1,len(morph),j)
				ap3 = ax3.imshow(b/c, origin='lower', aspect='auto', extent=(0, 13.8, 0, 4), interpolation='nearest', cmap=P.cm.gist_heat_r)
				ax3.set_xlabel(r'$t_q$ [Gyr]')
				ax3.set_ylabel(r'$\tau$ [Gyr]')
				ax3.text(5.0, 3.5, str(mass[m])+' < M < '+str(mass[m+1]), fontsize=9)
				cbar3 = P.colorbar(ap3)
				cbar3.set_label(r'fractional count of AGN', fontsize=9)
				print j, k, m, n
	P.tight_layout()
	P.savefig('fractional_count_AGN_binned_for '+str(colour[n])+'_by_morph_and_mass.png')

lum = agns['L_03']
colour = ['red', 'green', 'blue']
mass = [9.0, 9.5, 10.0, 10.5, 11.0, 11.5]
morph = ['t01_smooth_or_features_a01_smooth_flag', 't01_smooth_or_features_a02_features_or_disk_flag']
morphl = ['smooth', 'disc']
for n in range(len(colour)):
	j = 0
	print j
	P.figure(figsize=(7, 15))
	for m in range(len(mass)):
		if m == (len(mass)-1):
			pass
		else:
			for k in range(len(morph)):
				j+=1
				try:
					bdxagn = N.digitize(agns[N.where((agns['L_03']>0) & (agns[colour[n]]==1)& (agns[morph[k]]==1)& (agns['log_mstellar']>mass[m])& (agns['log_mstellar']<mass[m+1]))]['best_t'], bx)
					bdyagn = N.digitize(agns[N.where((agns['L_03']>0) & (agns[colour[n]]==1)& (agns[morph[k]]==1)& (agns['log_mstellar']>mass[m])& (agns['log_mstellar']<mass[m+1]))]['best_tau'], by)
				except ValueError:
					continue
				agn_lum = lum[N.where((agns['L_03']>0) & (agns[colour[n]]==1)& (agns[morph[k]]==1)& (agns['log_mstellar']>mass[m])& (agns['log_mstellar']<mass[m+1]))]
				avg_agn_bin = bin_avg(agn_lum, bdxagn, bdyagn)
				ax3 = P.subplot(len(mass)-1,len(morph),j)
				ap3 = ax3.imshow(N.log10(avg_agn_bin), origin='lower', aspect='auto', vmin=38, vmax=55, extent=(0, 13.8, 0, 4), interpolation='nearest', cmap=P.cm.gist_heat_r)
				ax3.set_xlabel(r'$t_q$ [Gyr]')
				ax3.set_ylabel(r'$\tau$ [Gyr]')
				ax3.text(5.0, 3.5, str(mass[m])+' < M < '+str(mass[m+1]), fontsize=9)
				cbar3 = P.colorbar(ap3)
				cbar3.set_label(r'average log $L[OIII]$ (erg/s)', fontsize=9)
				print j, k, m, n
	P.tight_layout()
	P.savefig('average_agn_lum_binned_for '+str(colour[n])+'_by_morph_and_mass.png')

smooth = data['t01_smooth_or_features_a01_smooth_flag']
smooth[N.where(smooth==-99)]=0
smooth_bin, frac_smooth, frac_all = bin_calc(smooth, bdx, bdy)
bin_plot(frac_all, frac_smooth, smooth_bin, 'smooth')


disc = data['t01_smooth_or_features_a02_features_or_disk_flag']
disc[N.where(disc==-99)]=0
disc_bin, frac_disc, frac_all = bin_calc(disc, bdx, bdy)
bin_plot(frac_all, frac_disc, disc_bin, 'disc')

agn = agns['agn?']
agn[N.where(agn==-32768)]=0
agn_bin, frac_agn, frac_all = bin_calc(agn, bdxagn, bdyagn)
bin_plot(frac_all, frac_agn, agn_bin, 'agn')

disc_agn = agns['t01_smooth_or_features_a02_features_or_disk_flag']
disc_agn_bin, frac_agn_disc, frac_discs = bin_calc(disc_agn, bdxagn, bdyagn)
bin_plot(frac_discs, frac_agn_disc, disc_agn_bin, 'agn in discs')

smooth_agn = agns['t01_smooth_or_features_a01_smooth_flag']
smooth_agn_bin, frac_agn_smooth, frac_smooth = bin_calc(smooth_agn, bdxagn, bdyagn)
bin_plot(frac_smooth, frac_agn_smooth, smooth_agn_bin, 'agn in ellipticals')

red = agn[N.where(agns['RED'] == 1)]
bdxagn = N.digitize(agns[N.where(agns['RED']==1)]['best_t'], bx)
bdyagn = N.digitize(agns[N.where(agns['RED']==1)]['best_tau'], by)
red_agn_bin, frac_agn_red, frac_red = bin_calc(red, bdxagn, bdyagn)
bin_plot(frac_red, frac_agn_red, red_agn_bin, 'agn in red galaxies')

green = agn[N.where(agns['GREEN'] == 1)]
bdxagn = N.digitize(agns[N.where(agns['GREEN']==1)]['best_t'], bx)
bdyagn = N.digitize(agns[N.where(agns['GREEN']==1)]['best_tau'], by)
green_agn_bin, frac_agn_green, frac_green = bin_calc(green, bdxagn, bdyagn)
bin_plot(frac_green, frac_agn_green, green_agn_bin, 'agn in green galaxies')

blue = agn[N.where(agns['BLUE'] == 1)]
bdxagn = N.digitize(agns[N.where(agns['BLUE']==1)]['best_t'], bx)
bdyagn = N.digitize(agns[N.where(agns['BLUE']==1)]['best_tau'], by)
blue_agn_bin, frac_agn_blue, frac_blue = bin_calc(blue, bdxagn, bdyagn)
bin_plot(frac_blue, frac_agn_blue, blue_agn_bin, 'agn in blue galaxies')

lum = agns['L_03']
bdxagn = N.digitize(agns[N.where((agns['L_03']>0) & (agns['RED']==1)& (agns['t01_smooth_or_features_a01_smooth_flag']==1))]['best_t'], bx)
bdyagn = N.digitize(agns[N.where((agns['L_03']>0) & (agns['RED']==1)& (agns['t01_smooth_or_features_a01_smooth_flag']==1))]['best_tau'], by)

agn_lum = lum[N.where((agns['L_03']>0) & (agns['RED']==1)& (agns['t01_smooth_or_features_a01_smooth_flag']==1))]
avg_agn_bin = bin_avg(agn_lum, bdxagn, bdyagn)
bin_plot(frac_all, N.log10(avg_agn_bin), agn_bin, 'average agn luminosity for red smooth')

lum = agns['L_03']
bdxagn = N.digitize(agns[N.where((agns['L_03']>0) & (agns['GREEN']==1)& (agns['t01_smooth_or_features_a01_smooth_flag']==1))]['best_t'], bx)
bdyagn = N.digitize(agns[N.where((agns['L_03']>0) & (agns['GREEN']==1)& (agns['t01_smooth_or_features_a01_smooth_flag']==1))]['best_tau'], by)

agn_lum = lum[N.where((agns['L_03']>0) & (agns['GREEN']==1)& (agns['t01_smooth_or_features_a01_smooth_flag']==1))]
avg_agn_bin = bin_avg(agn_lum, bdxagn, bdyagn)
bin_plot(frac_all, N.log10(avg_agn_bin), agn_bin, 'average agn luminosity for green smooth')

lum = agns['L_03']
bdxagn = N.digitize(agns[N.where((agns['L_03']>0) & (agns['BLUE']==1)& (agns['t01_smooth_or_features_a01_smooth_flag']==1))]['best_t'], bx)
bdyagn = N.digitize(agns[N.where((agns['L_03']>0) & (agns['BLUE']==1)& (agns['t01_smooth_or_features_a01_smooth_flag']==1))]['best_tau'], by)

agn_lum = lum[N.where((agns['L_03']>0) & (agns['BLUE']==1)& (agns['t01_smooth_or_features_a01_smooth_flag']==1))]
avg_agn_bin = bin_avg(agn_lum, bdxagn, bdyagn)
bin_plot(frac_all, N.log10(avg_agn_bin), agn_bin, 'average agn luminosity for blue smooth')


mstellar = agns['log_mstellar']

bdxagn = N.digitize(agns[N.where(agns['log_mstellar']>0)]['best_t'], bx)
bdyagn = N.digitize(agns[N.where(agns['log_mstellar']>0)]['best_tau'], by)

agn_mstellar = mstellar[N.where(agns['log_mstellar']>0)]
avg_mstellar_bin = bin_avg(agn_mstellar, bdxagn, bdyagn)
bin_plot(frac_all, avg_mstellar_bin, agn_bin, 'average stellar mass of agn')

sigma = agns['sigma']

bdxagn = N.digitize(agns[N.where(agns['sigma']>0)]['best_t'], bx)
bdyagn = N.digitize(agns[N.where(agns['sigma']>0)]['best_tau'], by)

agn_sigma = sigma[N.where(agns['sigma']>0)]
avg_sigma_bin = bin_avg(agn_sigma, bdxagn, bdyagn)
bin_plot(frac_all, avg_sigma_bin, agn_bin, 'average sigma of agn')

edge = data['t02_edgeon_a04_yes_flag']
edge[N.where(edge==-99)]=0
edge_bin, frac_edge, frac_all = bin_calc(edge, bdx, bdy)
bin_plot(disc_bin, edge_bin/(disc_bin*len(disc)), edge_bin, 'edge on')


bar = data['t03_bar_a06_bar_flag']
bar[N.where(bar==-99)]=0
bar_bin, frac_bar, frac_all = bin_calc(bar, bdx, bdy)
bin_plot(disc_bin, bar_bin/(disc_bin*len(disc)), bar_bin, 'barblue')

unbar = data['t03_bar_a07_no_bar_flag']
unbar[N.where(unbar==-99)]=0
unbar_bin, frac_unbar, frac_all = bin_calc(unbar, bdx, bdy)
bin_plot(disc_bin, unbar_bin/(disc_bin*len(disc)), unbar_bin, 'unbarred')

spiral = data['t04_spiral_a08_spiral_flag']
spiral[N.where(spiral==-99)]=0
spiral_bin, frac_spiral, frac_all = bin_calc(spiral, bdx, bdy)
bin_plot(disc_bin, spiral_bin/(disc_bin*len(disc)), spiral_bin, 'spiral')

nospiral = data['t04_spiral_a09_no_spiral_flag']
nospiral[N.where(nospiral==-99)]=0
nospiral_bin, frac_nospiral, frac_all = bin_calc(nospiral, bdx, bdy)
bin_plot(disc_bin, nospiral_bin/(disc_bin*len(disc)), nospiral_bin, 'no spiral')


spirals = data[N.where(data['t04_spiral_a08_spiral_flag']==1)]

odd = spirals['t11_arms_number_a31_1_flag'] + spirals['t11_arms_number_a33_3_flag']
odd_bin, frac_odd, frac_all = bin_calc(odd, bdx, bdy)
bin_plot(spiral_bin, odd_bin/(spiral_bin*len(spiral)), odd_bin, 'odd number of arms')

even = spirals['t11_arms_number_a32_2_flag'] +spirals['t11_arms_number_a34_4_flag']
even_bin, frac_even, frac_all = bin_calc(even, bdx, bdy)
bin_plot(spiral_bin, even_bin/(spiral_bin*len(spiral)), even_bin, 'even number of arms')

nobulge = data['t05_bulge_prominence_a10_no_bulge_flag']
nobulge[N.where(nobulge==-99)]=0
nobulge_bin, frac_nobulge, frac_all = bin_calc(nobulge, bdx, bdy)
bin_plot(disc_bin, nobulge_bin/(disc_bin*len(disc)), nobulge_bin, 'bulgeless')

dombulge = data['t05_bulge_prominence_a13_dominant_flag']
dombulge[N.where(dombulge==-99)]=0
dombulge_bin, frac_dombulge, frac_all = bin_calc(dombulge, bdx, bdy)
bin_plot(disc_bin, dombulge_bin/(disc_bin*len(disc)), dombulge_bin, 'dominant bulge')

odd = data['t06_odd_a14_yes_flag']
odd[N.where(odd==-99)]=0
odd_bin, frac_odd, frac_all = bin_calc(odd, bdx, bdy)
bin_plot(frac_all, frac_odd, odd_bin, 'odd')



