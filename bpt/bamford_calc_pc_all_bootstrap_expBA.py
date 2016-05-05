import numpy as N
import os
import time
import glob
from astropy.table import Table, vstack

dir1='/usersVol1/smethurst/hyper/samples/'
# l11 = os.listdir(dir1)
dir2 = '/usersVol1/smethurst/hyper/prob/'
# l22 = os.listdir(dir2)
# l22 = l22[2:]
pm = Table.read('bpt_identified_type2_agn_seyferts_no_liners_gz2_galex_matched_h_alpha_eqw_sigma_eddington_lum_ratios_stellar_mass_extra_axis_ratio_ab_MASS_SFR.fit', format='fits')

alow = pm[N.where(pm['expAB_g'] < 0.7)]
N.save('faceon_agn_matched_sample.npy', alow)
alows = []
alowp = []
lowpd = N.zeros((1,1))
lowps = N.zeros((1,1))
print 'probs...'
for n in range(len(alow)):
    alows.append(dir1+'samples_*_'+str(alow['col15'][n])+'_'+str(alow['col16'][n])+'.npy')
    alowp.append(dir2+'log_probability_*_'+str(alow['col15'][n])+'_'+str(alow['col16'][n])+'.npy')
    #lowpd = N.append(lowpd, alow['t01_smooth_or_features_a02_features_or_disk_debiased'][n]*N.ones((40000,1)), axis=0)
    #lowps = N.append(lowps, alow['t01_smooth_or_features_a01_smooth_debiased'][n]*N.ones((40000,1)), axis=0)
print 'globbing...'
alowg = map(glob.glob, alows)
alowgp = map(glob.glob, alowp)
print 'face on agn :', len(alowg)


ahigh = pm[N.where(pm['expAB_g'] > 0.7)]
N.save('edgeon_agn_matched_sample.npy', ahigh)
ahighs = []
ahighp = []
ahighg = []
ahighgp =[]
highpd = N.zeros((1,1))
highps = N.zeros((1,1))
print 'probs...'
for n in range(len(ahigh)):
    ahighs.append(dir1+'samples_*_'+str(ahigh['col15'][n])+'_'+str(ahigh['col16'][n])+'.npy')
    ahighp.append(dir2+'log_probability_*_'+str(ahigh['col15'][n])+'_'+str(ahigh['col16'][n])+'.npy')
    #highpd = N.append(highpd, ahigh['t01_smooth_or_features_a02_features_or_disk_debiased'][n]*N.ones((40000,1)), axis=0)
    #highps = N.append(highps, ahigh['t01_smooth_or_features_a01_smooth_debiased'][n]*N.ones((40000,1)), axis=0)
print 'globbing...'
ahighg = map(glob.glob, ahighs)
ahighgp = map(glob.glob, ahighp)
print 'edge on agn :', len(ahighg)

X = N.linspace(0, 14, 100)
Xs = X[:-1] + N.diff(X)
Y = N.linspace(0, 4, 100)
Ys = Y[:-1] + N.diff(Y)

count = 0

bf = N.zeros((1,6))
vf = N.zeros((1,2))




lowpcd = N.zeros((1000, 3))
medpcd = N.zeros((1000, 3))
highpcd = N.zeros((1000, 3))
lowpcs = N.zeros((1000, 3))
medpcs = N.zeros((1000, 3))
highpcs = N.zeros((1000, 3))

for k in range(1):
    sums=N.zeros((99,99))
    sumd=N.zeros((99,99))
    t = time.time()
    index = N.arange(len(alow))
    idx = N.random.permutation(index)
    alowgi = N.array(alowg)[idx].reshape(-1,)
    alowgpi = N.array(alowgp)[idx].reshape(-1,)
    alowi = alow[idx]
    s = N.array(map(N.load, alowgi))
    p = N.exp(N.array(map(N.load, alowgpi)))
    #s = N.array(s).reshape(-1,2)
    #p = N.array(N.exp(p)).reshape(-1,)
    #s = s[N.where(p>0.2)]
    #p = p[N.where(p>0.2)]
    #lowpds = lowpd[N.where(p>0.2)]
    #lowpss = lowps[N.where(p>0.2)]
    #wd = N.log(p)*lowpds.reshape(-1,)
    #ws = N.log(p)*lowpss.reshape(-1,)
    # del p
    # del lowpds
    # del lowpss
    for n in range(len(s)):
        ss, pp = s[n,:,:], p[n,:]
        ss = ss[N.where(pp > 0.2)]
        pp = pp[N.where(pp > 0.2)]
        if len(ss) == 0:
           pass
        else:
            Hs, Xss, Yss = N.histogram2d(ss[:,0], ss[:,1], bins=(X, Y), normed=True, weights=N.log(pp))
            sumd += (Hs*alowi['t01_smooth_or_features_a02_features_or_disk_debiased'][n].astype(float))
            sums += (Hs*alowi['t01_smooth_or_features_a01_smooth_debiased'][n].astype(float))
    agnla = N.nan_to_num(N.log10(sumd))
    nagnla = N.nan_to_num(N.log10(sums))
    agnla[agnla < 0] = 0
    nagnla[nagnla < 0] = 0
    sumdmi = agnla
    sumsmi = nagnla
    ys = N.sum(sumsmi, axis=0)
    xs = N.sum(sumsmi, axis=1)
    N.save('/users/smethurst/hyper/bootstrap/faceon_agn_smooth_tq_sum_bootstrap_all_'+str(k)+'.npy', xs)
    N.save('/users/smethurst/hyper/bootstrap/faceon_agn_smooth_tau_sum_bootstrap_all_'+str(k)+'.npy', ys)
    # ysmi = (ys-N.min(ys))/(N.max(ys)-N.min(ys))
    # lowpcs[k, 0] = (N.sum(ysmi[N.where(Ys < 1.0)])/N.sum(ysmi))*100
    # lowpcs[k, 1] =  (N.sum(ysmi[N.where(N.logical_and(Ys > 1.0, Ys < 2.0))])/N.sum(ysmi))*100
    # lowpcs[k, 2] =  (N.sum(ysmi[N.where(Ys > 2.0)])/N.sum(ysmi))*100
    yd = N.sum(sumdmi, axis=0)
    xd = N.sum(sumdmi, axis=1)
    N.save('/users/smethurst/hyper/bootstrap/faceon_agn_disc_tq_sum_bootstrap_all_'+str(k)+'.npy', xd)
    N.save('/users/smethurst/hyper/bootstrap/faceon_agn_disc_tau_sum_bootstrap_all_'+str(k)+'.npy', yd)
    # ydmi = (yd-N.min(yd))/(N.max(yd)-N.min(yd))
    # lowpcd[k,0] = (N.sum(ydmi[N.where(Ys < 1.0)])/N.sum(ydmi))*100
    # lowpcd[k,1] =  (N.sum(ydmi[N.where(N.logical_and(Ys > 1.0, Ys < 2.0))])/N.sum(ydmi))*100
    # lowpcd[k,2] =  (N.sum(ydmi[N.where(Ys > 2.0)])/N.sum(ydmi))*100
    # if k == 0:
    #     N.save('sum_weight_first_90pc_edgeon_smooth.npy', sums)
    #     N.save('sum_weight_first_90pc_edgeon_disc.npy', sumd)
    # else:
    #     pass
    print 'face on all agn boostrap: ', (k/1000.0)*100, '% complete'

for k in range(1000):
    sums=N.zeros((99,99))
    sumd=N.zeros((99,99))
    t = time.time()
    index = N.arange(len(alow))
    idx = N.random.permutation(index)[:0.9*len(index)]
    alowgi = N.array(alowg)[idx].reshape(-1,)
    alowgpi = N.array(alowgp)[idx].reshape(-1,)
    alowi = alow[idx]
    s = N.array(map(N.load, alowgi))
    p = N.exp(N.array(map(N.load, alowgpi)))
    #s = N.array(s).reshape(-1,2)
    #p = N.array(N.exp(p)).reshape(-1,)
    #s = s[N.where(p>0.2)]
    #p = p[N.where(p>0.2)]
    #lowpds = lowpd[N.where(p>0.2)]
    #lowpss = lowps[N.where(p>0.2)]
    #wd = N.log(p)*lowpds.reshape(-1,)
    #ws = N.log(p)*lowpss.reshape(-1,)
    # del p
    # del lowpds
    # del lowpss
    for n in range(len(s)):
        ss, pp = s[n,:,:], p[n,:]
        ss = ss[N.where(pp > 0.2)]
        pp = pp[N.where(pp > 0.2)]
        if len(ss) == 0:
           pass
        else:
            Hs, Xss, Yss = N.histogram2d(ss[:,0], ss[:,1], bins=(X, Y), normed=True, weights=N.log(pp))
            sumd += (Hs*alowi['t01_smooth_or_features_a02_features_or_disk_debiased'][n].astype(float))
            sums += (Hs*alowi['t01_smooth_or_features_a01_smooth_debiased'][n].astype(float))
    agnla = N.nan_to_num(N.log10(sumd))
    nagnla = N.nan_to_num(N.log10(sums))
    agnla[agnla < 0] = 0
    nagnla[nagnla < 0] = 0
    sumdmi = agnla
    sumsmi = nagnla
    ys = N.sum(sumsmi, axis=0)
    xs = N.sum(sumsmi, axis=1)
    N.save('/users/smethurst/hyper/bootstrap/faceon_agn_smooth_tq_sum_bootstrap_'+str(k)+'.npy', xs)
    N.save('/users/smethurst/hyper/bootstrap/faceon_agn_smooth_tau_sum_bootstrap_'+str(k)+'.npy', ys)
    # ysmi = (ys-N.min(ys))/(N.max(ys)-N.min(ys))
    # lowpcs[k, 0] = (N.sum(ysmi[N.where(Ys < 1.0)])/N.sum(ysmi))*100
    # lowpcs[k, 1] =  (N.sum(ysmi[N.where(N.logical_and(Ys > 1.0, Ys < 2.0))])/N.sum(ysmi))*100
    # lowpcs[k, 2] =  (N.sum(ysmi[N.where(Ys > 2.0)])/N.sum(ysmi))*100
    yd = N.sum(sumdmi, axis=0)
    xd = N.sum(sumdmi, axis=1)
    N.save('/users/smethurst/hyper/bootstrap/faceon_agn_disc_tq_sum_bootstrap_'+str(k)+'.npy', xd)
    N.save('/users/smethurst/hyper/bootstrap/faceon_agn_disc_tau_sum_bootstrap_'+str(k)+'.npy', yd)
    # ydmi = (yd-N.min(yd))/(N.max(yd)-N.min(yd))
    # lowpcd[k,0] = (N.sum(ydmi[N.where(Ys < 1.0)])/N.sum(ydmi))*100
    # lowpcd[k,1] =  (N.sum(ydmi[N.where(N.logical_and(Ys > 1.0, Ys < 2.0))])/N.sum(ydmi))*100
    # lowpcd[k,2] =  (N.sum(ydmi[N.where(Ys > 2.0)])/N.sum(ydmi))*100
    # if k == 0:
    #     N.save('sum_weight_first_90pc_edgeon_smooth.npy', sums)
    #     N.save('sum_weight_first_90pc_edgeon_disc.npy', sumd)
    # else:
    #     pass
    print 'face on agn boostrap: ', (k/1000.0)*100, '% complete'

# N.save('edgeon_pc_r_i_s_summed_smooth.npy', lowpcs)
# N.save('edgeon_pc_r_i_s_summed_disc.npy', lowpcd)

sums=N.zeros((99,99))
sumd=N.zeros((99,99))


for k in range(1):
    sums=N.zeros((99,99))
    sumd=N.zeros((99,99))
    t = time.time()
    index = N.arange(len(ahigh))
    idx = N.random.permutation(index)
    ahighgi = N.array(ahighg)[idx].reshape(-1,)
    ahighgpi = N.array(ahighgp)[idx].reshape(-1,)
    ahighi = ahigh[idx]
    s = N.array(map(N.load, ahighgi))
    p = N.exp(N.array(map(N.load, ahighgpi)))
    #s = N.array(s).reshape(-1,2)
    #p = N.array(N.exp(p)).reshape(-1,)
    #s = s[N.where(p>0.2)]
    #p = p[N.where(p>0.2)]
    #lowpds = lowpd[N.where(p>0.2)]
    #lowpss = lowps[N.where(p>0.2)]
    #wd = N.log(p)*lowpds.reshape(-1,)
    #ws = N.log(p)*lowpss.reshape(-1,)
    # del p
    # del lowpds
    # del lowpss
    for n in range(len(s)):
        ss, pp = s[n,:,:], p[n,:]
        ss = ss[N.where(pp > 0.2)]
        pp = pp[N.where(pp > 0.2)]
        if len(s) == 0:
           pass
        else:
            Hs, Xss, Yss = N.histogram2d(ss[:,0], ss[:,1], bins=(X, Y), normed=True, weights=N.log(pp))
            sumd += (Hs*ahighi['t01_smooth_or_features_a02_features_or_disk_debiased'][n].astype(float))
            sums += (Hs*ahighi['t01_smooth_or_features_a01_smooth_debiased'][n].astype(float))
    agnmi = N.nan_to_num(N.log10(sumd))
    nagnmi = N.nan_to_num(N.log10(sums))
    agnmi[agnmi < 0] = 0
    nagnmi[nagnmi < 0] = 0
    sumdmi = agnmi
    sumsmi =nagnmi
    ys = N.sum(sumsmi, axis=0)
    xs = N.sum(sumsmi, axis=1)
    N.save('/users/smethurst/hyper/bootstrap/edgeon_agn_smooth_tq_sum_bootstrap_all_'+str(k)+'.npy', xs)
    N.save('/users/smethurst/hyper/bootstrap/edgeon_agn_smooth_tau_sum_bootstrap_all_'+str(k)+'.npy', ys)
    # ysmi = (ys-N.min(ys))/(N.max(ys)-N.min(ys))
    # lowpcs[k, 0] = (N.sum(ysmi[N.where(Ys < 1.0)])/N.sum(ysmi))*100
    # lowpcs[k, 1] =  (N.sum(ysmi[N.where(N.logical_and(Ys > 1.0, Ys < 2.0))])/N.sum(ysmi))*100
    # lowpcs[k, 2] =  (N.sum(ysmi[N.where(Ys > 2.0)])/N.sum(ysmi))*100
    yd = N.sum(sumdmi, axis=0)
    xd = N.sum(sumdmi, axis=1)
    N.save('/users/smethurst/hyper/bootstrap/edgeon_agn_disc_tq_sum_bootstrap_all_'+str(k)+'.npy', xd)
    N.save('/users/smethurst/hyper/bootstrap/edgeon_agn_disc_tau_sum_bootstrap_all_'+str(k)+'.npy', yd)
    # ydmi = (yd-N.min(yd))/(N.max(yd)-N.min(yd))
    # lowpcd[k,0] = (N.sum(ydmi[N.where(Ys < 1.0)])/N.sum(ydmi))*100
    # lowpcd[k,1] =  (N.sum(ydmi[N.where(N.logical_and(Ys > 1.0, Ys < 2.0))])/N.sum(ydmi))*100
    # lowpcd[k,2] =  (N.sum(ydmi[N.where(Ys > 2.0)])/N.sum(ydmi))*100
    # if k == 0:
    #     N.save('sum_weight_first_90pc_edgeon_smooth.npy', sums)
    #     N.save('sum_weight_first_90pc_edgeon_disc.npy', sumd)
    # else:
    #     pass
    print 'edge on all agn boostrap: ', (k/1000.0)*100, '% complete'
# N.save('edgeon_pc_r_i_s_summed_smooth.npy', highpcs)
# N.save('edgeon_pc_r_i_s_summed_disc.npy', highpcd)

for k in range(1000):
    sums=N.zeros((99,99))
    sumd=N.zeros((99,99))
    t = time.time()
    index = N.arange(len(ahigh))
    idx = N.random.permutation(index)[:0.9*len(index)]
    ahighgi = N.array(ahighg)[idx].reshape(-1,)
    ahighgpi = N.array(ahighgp)[idx].reshape(-1,)
    ahighi = ahigh[idx]
    s = N.array(map(N.load, ahighgi))
    p = N.exp(N.array(map(N.load, ahighgpi)))
    #s = N.array(s).reshape(-1,2)
    #p = N.array(N.exp(p)).reshape(-1,)
    #s = s[N.where(p>0.2)]
    #p = p[N.where(p>0.2)]
    #lowpds = lowpd[N.where(p>0.2)]
    #lowpss = lowps[N.where(p>0.2)]
    #wd = N.log(p)*lowpds.reshape(-1,)
    #ws = N.log(p)*lowpss.reshape(-1,)
    # del p
    # del lowpds
    # del lowpss
    for n in range(len(s)):
        ss, pp = s[n,:,:], p[n,:]
        ss = ss[N.where(pp > 0.2)]
        pp = pp[N.where(pp > 0.2)]
        if len(s) == 0:
           pass
        else:
            Hs, Xss, Yss = N.histogram2d(ss[:,0], ss[:,1], bins=(X, Y), normed=True, weights=N.log(pp))
            sumd += (Hs*ahighi['t01_smooth_or_features_a02_features_or_disk_debiased'][n].astype(float))
            sums += (Hs*ahighi['t01_smooth_or_features_a01_smooth_debiased'][n].astype(float))
    agnmi = N.nan_to_num(N.log10(sumd))
    nagnmi = N.nan_to_num(N.log10(sums))
    agnmi[agnmi < 0] = 0
    nagnmi[nagnmi < 0] = 0
    sumdmi = agnmi
    sumsmi =nagnmi
    ys = N.sum(sumsmi, axis=0)
    xs = N.sum(sumsmi, axis=1)
    N.save('/users/smethurst/hyper/bootstrap/edgeon_agn_smooth_tq_sum_bootstrap_'+str(k)+'.npy', xs)
    N.save('/users/smethurst/hyper/bootstrap/edgeon_agn_smooth_tau_sum_bootstrap_'+str(k)+'.npy', ys)
    # ysmi = (ys-N.min(ys))/(N.max(ys)-N.min(ys))
    # lowpcs[k, 0] = (N.sum(ysmi[N.where(Ys < 1.0)])/N.sum(ysmi))*100
    # lowpcs[k, 1] =  (N.sum(ysmi[N.where(N.logical_and(Ys > 1.0, Ys < 2.0))])/N.sum(ysmi))*100
    # lowpcs[k, 2] =  (N.sum(ysmi[N.where(Ys > 2.0)])/N.sum(ysmi))*100
    yd = N.sum(sumdmi, axis=0)
    xd = N.sum(sumdmi, axis=1)
    N.save('/users/smethurst/hyper/bootstrap/edgeon_agn_disc_tq_sum_bootstrap_'+str(k)+'.npy', xd)
    N.save('/users/smethurst/hyper/bootstrap/edgeon_agn_disc_tau_sum_bootstrap_'+str(k)+'.npy', yd)
    # ydmi = (yd-N.min(yd))/(N.max(yd)-N.min(yd))
    # lowpcd[k,0] = (N.sum(ydmi[N.where(Ys < 1.0)])/N.sum(ydmi))*100
    # lowpcd[k,1] =  (N.sum(ydmi[N.where(N.logical_and(Ys > 1.0, Ys < 2.0))])/N.sum(ydmi))*100
    # lowpcd[k,2] =  (N.sum(ydmi[N.where(Ys > 2.0)])/N.sum(ydmi))*100
    # if k == 0:
    #     N.save('sum_weight_first_90pc_edgeon_smooth.npy', sums)
    #     N.save('sum_weight_first_90pc_edgeon_disc.npy', sumd)
    # else:
    #     pass
    print 'edge on agn boostrap: ', (k/1000.0)*100, '% complete'
# N.save('edgeon_pc_r_i_s_summed_smooth.npy', highpcs)
# N.save('edgeon_pc_r_i_s_summed_disc.npy', highpcd)