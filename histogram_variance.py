#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from scipy.spatial import distance
from sklearn.neighbors import BallTree


extents = {"US":[-140,-60,10,75], "EU":[-35,60,10,70], "SAm":[-120,-60,-70,-5], "SAs":[60,140,0,60], "Alaska":[-170,-140, 50, 75], "Global":[-175,175,-85,85]}
region = "EU"

depths = np.arange(0, 2900, 100)
depths = np.append(depths, 2891)

fmin = 0.10
fmax = 0.20
# test_file = f"Plotting_file_{fmin:.2f}_{fmax:.2f}_1000.0km_variance.txt"
# df_test = pd.read_csv(test_file, sep=',', index_col=False)

# rand_variances = []
freq_band_labels = ['0.10 $-$ 0.20', '0.15 $-$ 0.30', '0.20 $-$ 0.40']


depths_hists_r_freqs = []
depths_hists_s_freqs = []

# for f in [[0.1, 0.2], [0.15, 0.3], [0.2, 0.4]]:
#
#     fmin = f[0]
#     fmax = f[1]

test_file = f"Plotting_file_{fmin:.2f}_{fmax:.2f}_1000.0km_variance.txt"
df_test = pd.read_csv(test_file, sep=',', index_col=False)
# df_test = df_test[df_test['mag'] > 0.5]
rand_variances = []
for i in range(1000):
    sub_sample = df_test.sample(20)
    mags_random = sub_sample['mag'].to_numpy()
    vecs_x = sub_sample['del_x_slow'].to_numpy()
    vecs_y = sub_sample['del_y_slow'].to_numpy()

    mean_x = np.mean(vecs_x)
    mean_y = np.mean(vecs_y)

    diff_vecs = np.array(list(zip(vecs_x, vecs_y)))
    means = np.array([[mean_x,mean_y]])

    dists_from_mean = distance.cdist(diff_vecs, means, metric="euclidean")

    vect_variance = np.mean(np.power(dists_from_mean, 2))

    variance_sample = np.var(mags_random)
    rand_variances.append(vect_variance)
    # rand_variances.append(variance_sample)


variance_estimate = np.mean(rand_variances)
variance_std = np.std(rand_variances)
print("mean_variance_sample", variance_estimate)
print("variance_std", np.std(rand_variances))

## get 95% variance value
var_95 = variance_estimate - (2*variance_std)
print("95% variance value", var_95)


df_test = df_test[(df_test['stlo_mean'] > extents[region][0]) & (df_test['stlo_mean'] < extents[region][1]) & (df_test['stla_mean'] > extents[region][2]) & (df_test['stla_mean'] < extents[region][3]) ]
# df_test = df_test[df_test['mag'] > 0.5]

r_reloc_variances_wt_depth = np.empty((len(df_test), len(depths)))
s_reloc_variances_wt_depth = np.empty((len(df_test), len(depths)))
r_gcp_variances_wt_depth = np.empty((len(df_test), len(depths)))
s_gcp_variances_wt_depth = np.empty((len(df_test), len(depths)))


for i, depth in enumerate(depths):

    Plotting_file = f"./Plotting_file_{fmin:.2f}_{fmax:.2f}_{depth}.0km_variance.txt"
    df = pd.read_csv(Plotting_file, sep=',', index_col=False)

    # filter based on station location
    df = df[(df['stlo_mean'] > extents[region][0]) & (df['stlo_mean'] < extents[region][1]) & (df['stla_mean'] > extents[region][2]) & (df['stla_mean'] < extents[region][3]) ]
    # df = df[df['mag'] > 0.5]


    stlas = df['stla_mean'].to_numpy()
    stlos = df['stlo_mean'].to_numpy()
    evlas = df['evla'].to_numpy()
    evlos = df['evlo'].to_numpy()
    evdps = df['evdp'].to_numpy()
    var_s_gcp = df['vec_var_s_gcp'].to_numpy()
    var_r_gcp = df['vec_var_r_gcp'].to_numpy()
    var_s_reloc = df['vec_var_s_reloc'].to_numpy()
    var_r_reloc = df['vec_var_r_reloc'].to_numpy()

    r_reloc_variances_wt_depth[:,i] = var_r_reloc
    s_reloc_variances_wt_depth[:,i] = var_s_reloc
    r_gcp_variances_wt_depth[:,i] = var_r_gcp
    s_gcp_variances_wt_depth[:,i] = var_s_gcp

r_reloc_variances_wt_depth[r_reloc_variances_wt_depth == -1] = 10
s_reloc_variances_wt_depth[s_reloc_variances_wt_depth == -1] = 10

r_gcp_variances_wt_depth[r_gcp_variances_wt_depth == -1] = 10
s_gcp_variances_wt_depth[s_gcp_variances_wt_depth == -1] = 10

mags = df['mag'].to_numpy()


indexes_r = []
indexes_s = []
depth_variances_r = []
depth_variances_s = []

good_cross_sects_r = []
good_cross_sects_s = []

depths_var_r = []
depths_var_s = []

depths_hist_r = []
depths_hist_s = []


los_r = []
las_r = []

los_s = []
las_s = []

for j,cross_sect in enumerate(r_reloc_variances_wt_depth):
     if cross_sect.min() < var_95:
         depths_variance = depths[np.where(cross_sect < var_95)[0]]
         depth_variance = depths[np.where(cross_sect==cross_sect.min())[0]]

         # good_cross_sects_r.append(cross_sect[np.where(cross_sect < var_95)[0]])
         good_cross_sects_r.append(cross_sect)

         depths_var_r.append(depths_variance)
         depths_hist_r.extend(depths_variance)

         indexes_r.append(j)

         # get reloc points
         d_var = depth_variance[0]
         Plotting_file = f"./Plotting_file_{fmin:.2f}_{fmax:.2f}_{d_var}.0km_variance.txt"
         df = pd.read_csv(Plotting_file, sep=',', index_col=False)
         df = df[(df['stlo_mean'] > extents[region][0]) & (df['stlo_mean'] < extents[region][1]) & (df['stla_mean'] > extents[region][2]) & (df['stla_mean'] < extents[region][3]) ]

         r_gcp_pierce_las = df['r_pierce_la'].to_numpy()
         r_gcp_pierce_los = df['r_pierce_lo'].to_numpy()

         if len(depths_variance) <= 10:
             depth_variances_r.append(depth_variance[0])
             los_r.append(r_gcp_pierce_los[j])
             las_r.append(r_gcp_pierce_las[j])



for k,cross_sect in enumerate(s_reloc_variances_wt_depth):
     if cross_sect.min() < var_95:
         depths_variance = depths[np.where(cross_sect < var_95)[0]]
         depth_variance = depths[np.where(cross_sect==cross_sect.min())[0]]

         # good_cross_sects_s.append(cross_sect[np.where(cross_sect < var_95)[0]])
         good_cross_sects_s.append(cross_sect)

         depths_var_s.append(depths_variance)
         depths_hist_s.extend(depths_variance)
         indexes_s.append(k)

         d_var = depth_variance[0]
         Plotting_file = f"./Plotting_file_{fmin:.2f}_{fmax:.2f}_{d_var}.0km_variance.txt"
         df = pd.read_csv(Plotting_file, sep=',', index_col=False)
         df = df[(df['stlo_mean'] > extents[region][0]) & (df['stlo_mean'] < extents[region][1]) & (df['stla_mean'] > extents[region][2]) & (df['stla_mean'] < extents[region][3]) ]

         s_gcp_pierce_las = df['s_pierce_la'].to_numpy()
         s_gcp_pierce_los = df['s_pierce_lo'].to_numpy()

         if len(depths_variance) <= 10:

            los_s.append(s_gcp_pierce_los[k])
            las_s.append(s_gcp_pierce_las[k])
            depth_variances_s.append(depth_variance[0])

depths_hist_s = np.array(depths_hist_s)
depths_hist_r = np.array(depths_hist_r)

depths_hists_r_freqs.append(depths_hist_r)
depths_hists_s_freqs.append(depths_hist_s)

## plot histogram of depths which have variance less than estimate

depths = np.arange(0, 3100, 100)
# depths = np.append(depths, 2891)
depths -= 50

# depths = np.append(depths, 2941)

print(depths_hists_r_freqs)

plt.style.use('ggplot')
colors = ['red', 'blue', 'green']

# fig = plt.figure(figsize=(10,8))
# ax = fig.add_subplot(111)
# ax.hist(depths_hists_r_freqs, bins=depths, histtype='bar', color=colors, label=freq_band_labels, edgecolor='white', align='mid')
# plt.savefig(f"histogram_variance_{fmin}_{fmax}_{region}_r_side.pdf")
# plt.show()
#
# fig = plt.figure(figsize=(10,8))
# ax = fig.add_subplot(111)
# ax.hist(depths_hists_s_freqs, bins=depths, histtype='bar', color=colors, label=freq_band_labels, edgecolor='white', align='mid')
# plt.savefig(f"histogram_variance_{fmin}_{fmax}_{region}_s_side.pdf")
# plt.show()


fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)
ax.hist(depths_hist_r, bins=depths, histtype='bar', color='grey', edgecolor='white', align='mid')
ax.set_xlabel(f"Depths (km)")
ax.set_ylabel(f"Number of values in bin")
ax.set_title(f"{region} Receiver-Side: {fmin:.2f} $-$ {fmax:.2f} Hz")
plt.savefig(f"histogram_variance_{fmin:.2f}_{fmax:.2f}_{region}_r_side_new.pdf")
plt.show()

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)
ax.hist(depths_hist_s, bins=depths, histtype='bar', color='grey', edgecolor='white', align='mid')
ax.set_xlabel(f"Depths (km)")
ax.set_ylabel(f"Number of values in bin")
ax.set_title(f"{region} Source-Side: {fmin:.2f} $-$ {fmax:.2f} Hz")
plt.savefig(f"histogram_variance_{fmin:.2f}_{fmax:.2f}_{region}_s_side_new.pdf")
plt.show()





## Now bin the depths??
#
# lats_lons_r = np.array(list(zip(np.deg2rad(las_r), np.deg2rad(los_r))))
# lats_lons_s = np.array(list(zip(np.deg2rad(las_s), np.deg2rad(los_s))))
#
# tree_r = BallTree(
#     lats_lons_r, leaf_size=lats_lons_r.shape[0] / 2, metric="haversine"
# )
#
# tree_s = BallTree(
#     lats_lons_s, leaf_size=lats_lons_s.shape[0] / 2, metric="haversine"
# )
#
# bin_radius = 5
#
# spacing = 10
# min_lo = -180 + bin_radius
# max_lo = 180 - bin_radius
# n_lo = (max_lo - min_lo) / spacing
#
# min_la = -90 + bin_radius
# max_la = 90 - bin_radius
# n_la = (max_la - min_la) / spacing
#
# los_centre_grid = np.linspace(min_lo, max_lo, int(n_lo+1))
# las_centre_grid = np.linspace(min_la, max_la, int(n_la+1))
# los,las = np.meshgrid(los_centre_grid, las_centre_grid)
# los = los.flatten()
# las = las.flatten()
# la_lo_grid = np.array(list(zip(np.deg2rad(las), np.deg2rad(los))))
#
# indices_r = tree_r.query_radius(la_lo_grid, np.radians(bin_radius), count_only=False, return_distance = False)
# counts_r = tree_r.query_radius(la_lo_grid, np.radians(bin_radius), count_only=True, return_distance = False)
#
# indices_s = tree_s.query_radius(la_lo_grid, np.radians(bin_radius), count_only=False, return_distance = False)
# counts_s = tree_s.query_radius(la_lo_grid, np.radians(bin_radius), count_only=True, return_distance = False)
#
# depths_var_r = np.array(depths_var_r)
# good_cross_sects_r = np.array(good_cross_sects_r)
# good_cross_sects_s = np.array(good_cross_sects_s)
# print(depths_var_r)
#
# for point_bin_r in indices_r:
#     depths_bin = depths_var_r[point_bin_r]
#     cross_s = good_cross_sects_r[point_bin_r]
#     print(depths_bin)
#     if len(depths_bin) > 0:
#         for a, depths_temp in enumerate(depths_bin):
#             no_value = np.where(cross_s[a] == 10)[0]
#             cross_s[a][no_value] = np.nan
#             print(depths)
#             # depths[no_value] = np.nan
#             plt.plot(cross_s[a], depths)
#             plt.xlim(0,5)
#             plt.ylim(0,3000)
#             plt.axvline(x=var_95, ymin=0, ymax=3000, label='significant variance threshold')
#         plt.show()
#
# exit()
#
# fig = plt.figure(figsize=(12,8))
#
# ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
# ax.add_feature(cartopy.feature.LAND, edgecolor='black', color='lightgrey')
#
# a = ax.scatter(
# los_r,
# las_r,
# transform=ccrs.PlateCarree(),
# marker='o',
# s=80,
# linewidth=0.1,
# c=depth_variances_r,
# edgecolors='black',
# cmap='viridis'
# )
#
# plt.colorbar(a, label='Depth (km)')
#
# ax.set_title(f"")
# ax.set_xlabel('Longitude ($^{\circ}$)')
# ax.set_ylabel('Latitude ($^{\circ}$)')
#
# # ax.set_extent(extents_s['SA'])
# plt.legend(loc='best')
#
# gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                   linewidth=1, color='gray', alpha=0.5, linestyle='--')
#
# gl.xlines = True
# gl.ylines = True
# gl.xlabels_top = False
# gl.ylabels_right = True
# gl.xlabels_bottom = True
# gl.ylabels_left = True
#
# # plt.savefig(f"./img{i}.png", dpi=600)
# # plt.savefig(f"./vector_bin_{depth}.pdf")
# # os.rename(f"./img{i}.png", f"{dir}/img{i}.png")
# #     os.rename(f"./vector_bin_{depth}.pdf", f"{dir}/vector_bin_{depth}.pdf")
#
# plt.savefig(f"./receiver_minimum_variance_depth_{fmin:.2f}_{fmax:.2f}_{region}.pdf")
#
#
# fig = plt.figure(figsize=(12,8))
#
# ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
# ax.add_feature(cartopy.feature.LAND, edgecolor='black', color='lightgrey')
#
# a = ax.scatter(
# los_s,
# las_s,
# transform=ccrs.PlateCarree(),
# marker='o',
# s=80,
# linewidth=0.1,
# c=depth_variances_s,
# edgecolors='black',
# cmap='viridis'
# )
#
# plt.colorbar(a)
#
# ax.set_title(f"")
# ax.set_xlabel('Longitude ($^{\circ}$)')
# ax.set_ylabel('Latitude ($^{\circ}$)')
#
# # ax.set_extent(extents_s['SA'])
# plt.legend(loc='best')
#
# gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                   linewidth=1, color='gray', alpha=0.5, linestyle='--')
#
# gl.xlines = True
# gl.ylines = True
# gl.xlabels_top = False
# gl.ylabels_right = True
# gl.xlabels_bottom = True
# gl.ylabels_left = True

# os.rename(f"./img{i}.png", f"{dir}/img{i}.png")
#     os.rename(f"./vector_bin_{depth}.pdf", f"{dir}/vector_bin_{depth}.pdf")
