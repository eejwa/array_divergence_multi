#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from scipy.spatial import distance

extents = {"US":[-140,-60,10,75], "EU":[-35,60,10,70], "SAm":[-120,-60,-70,-5], "SAs":[60,140,0,60], "Alaska":[-170,-140, 50, 75]}
region = "EU"

depths = np.arange(0, 2900, 100)
depths = np.append(depths, 2891)

fmin = 0.20
fmax = 0.40
test_file = f"Plotting_file_{fmin:.2f}_{fmax:.2f}_1000.0km_variance.txt"
df_test = pd.read_csv(test_file, sep=',', index_col=False)

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

    print(vect_variance)
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
df_test = df_test[df_test['mag'] < 0.5]
r_reloc_variances_wt_depth = np.empty((len(df_test), len(depths)))
s_reloc_variances_wt_depth = np.empty((len(df_test), len(depths)))
r_gcp_variances_wt_depth = np.empty((len(df_test), len(depths)))
s_gcp_variances_wt_depth = np.empty((len(df_test), len(depths)))


for i, depth in enumerate(depths):

    Plotting_file = f"./Plotting_file_{fmin:.2f}_{fmax:.2f}_{depth}.0km_variance.txt"
    df = pd.read_csv(Plotting_file, sep=',', index_col=False)

    # filter based on station location
    df = df[(df['stlo_mean'] > extents[region][0]) & (df['stlo_mean'] < extents[region][1]) & (df['stla_mean'] > extents[region][2]) & (df['stla_mean'] < extents[region][3]) ]

    # filter based on mag
    df = df[df['mag'] < 0.5]

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

print(r_gcp_variances_wt_depth.max())

mags = df['mag'].to_numpy()


indexes_r = []
indexes_s = []
depth_variances_r = []
depth_variances_s = []

print(len(df), r_reloc_variances_wt_depth.shape)

los_r = []
las_r = []

los_s = []
las_s = []

for j,cross_sect in enumerate(r_reloc_variances_wt_depth):
     if cross_sect.min() < var_95:
         depths_variance = depths[np.where(cross_sect < var_95)[0]]
         depth_variance = depths[np.where(cross_sect==cross_sect.min())[0]]

         print(j, cross_sect.min(), depths[np.where(cross_sect==cross_sect.min())[0]], var_95)
         print(depths_variance)
         indexes_r.append(j)

         # get reloc points
         d_var = depth_variance[0]
         Plotting_file = f"./Plotting_file_{fmin:.2f}_{fmax:.2f}_{d_var}.0km_variance.txt"
         df = pd.read_csv(Plotting_file, sep=',', index_col=False)
         df = df[(df['stlo_mean'] > extents[region][0]) & (df['stlo_mean'] < extents[region][1]) & (df['stla_mean'] > extents[region][2]) & (df['stla_mean'] < extents[region][3]) ]

         r_reloc_pierce_las = df['r_reloc_pierce_la'].to_numpy()
         r_reloc_pierce_los = df['r_reloc_pierce_lo'].to_numpy()

         if len(depths_variance) <= 10:
             depth_variances_r.append(depth_variance[0])
             los_r.append(r_reloc_pierce_los[j])
             las_r.append(r_reloc_pierce_las[j])



for k,cross_sect in enumerate(s_reloc_variances_wt_depth):
     if cross_sect.min() < var_95:
         depths_variance = depths[np.where(cross_sect < var_95)[0]]
         depth_variance = depths[np.where(cross_sect==cross_sect.min())[0]]
         print(k, cross_sect.min(), depths[np.where(cross_sect==cross_sect.min())[0]], var_95)

         indexes_s.append(k)

         d_var = depth_variance[0]
         Plotting_file = f"./Plotting_file_{fmin:.2f}_{fmax:.2f}_{d_var}.0km_variance.txt"
         df = pd.read_csv(Plotting_file, sep=',', index_col=False)
         df = df[(df['stlo_mean'] > extents[region][0]) & (df['stlo_mean'] < extents[region][1]) & (df['stla_mean'] > extents[region][2]) & (df['stla_mean'] < extents[region][3]) ]

         s_reloc_pierce_las = df['s_reloc_pierce_la'].to_numpy()
         s_reloc_pierce_los = df['s_reloc_pierce_lo'].to_numpy()

         if len(depths_variance) <= 10:

            los_s.append(s_reloc_pierce_los[k])
            las_s.append(s_reloc_pierce_las[k])
            depth_variances_s.append(depth_variance[0])



same_indices = np.intersect1d(indexes_r, indexes_s)
different_indices_r = np.setdiff1d(indexes_r, indexes_s)
different_indices_s = np.setdiff1d(indexes_s, indexes_r)

print(same_indices)
print(different_indices_r)
print(different_indices_s)

## screw it lets plot these up coloured by depth!!


fig = plt.figure(figsize=(12,8))

ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND, edgecolor='black', color='lightgrey')

a = ax.scatter(
los_r,
las_r,
transform=ccrs.PlateCarree(),
marker='o',
s=80,
linewidth=0.1,
c=depth_variances_r,
edgecolors='black',
cmap='viridis'
)

plt.colorbar(a, label='Depth (km)')

ax.set_title(f"")
ax.set_xlabel('Longitude ($^{\circ}$)')
ax.set_ylabel('Latitude ($^{\circ}$)')

# ax.set_extent(extents_s['SA'])
plt.legend(loc='best')

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')

gl.xlines = True
gl.ylines = True
gl.xlabels_top = False
gl.ylabels_right = True
gl.xlabels_bottom = True
gl.ylabels_left = True

# plt.savefig(f"./img{i}.png", dpi=600)
# plt.savefig(f"./vector_bin_{depth}.pdf")
# os.rename(f"./img{i}.png", f"{dir}/img{i}.png")
#     os.rename(f"./vector_bin_{depth}.pdf", f"{dir}/vector_bin_{depth}.pdf")

# plt.savefig(f"./receiver_minimum_variance_depth_{fmin:.2f}_{fmax:.2f}_{region}.pdf")

plt.show()
fig = plt.figure(figsize=(12,8))

ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND, edgecolor='black', color='lightgrey')

a = ax.scatter(
los_s,
las_s,
transform=ccrs.PlateCarree(),
marker='o',
s=80,
linewidth=0.1,
c=depth_variances_s,
edgecolors='black',
cmap='viridis'
)

plt.colorbar(a)

ax.set_title(f"")
ax.set_xlabel('Longitude ($^{\circ}$)')
ax.set_ylabel('Latitude ($^{\circ}$)')

# ax.set_extent(extents_s['SA'])
plt.legend(loc='best')

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')

gl.xlines = True
gl.ylines = True
gl.xlabels_top = False
gl.ylabels_right = True
gl.xlabels_bottom = True
gl.ylabels_left = True
plt.show()
# os.rename(f"./img{i}.png", f"{dir}/img{i}.png")
#     os.rename(f"./vector_bin_{depth}.pdf", f"{dir}/vector_bin_{depth}.pdf")
