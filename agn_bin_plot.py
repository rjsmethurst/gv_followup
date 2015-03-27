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