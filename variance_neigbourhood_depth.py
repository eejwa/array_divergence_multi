#!/usr/bin/env python

import cartopy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd
import numpy as np
import matplotlib.ticker as mticker
from sklearn.cluster import KMeans
import os
import shutil
from sklearn.neighbors import BallTree
from scipy.spatial import distance

depths = np.arange(0, 2900, 100)
depths = np.append(depths, 2891)

fmin = 0.10
fmax = 0.20
cmap = plt.cm.get_cmap('summer')
gird_lons = np.arange(-360,375, 15)
gird_lats = np.arange(-90,105, 15)

bin_distance = 100 # km

for d,depth in enumerate(depths):

    rad = 6371 - depth

    # calculate bin distance in radians
    bin_radians = bin_distance / rad

    Plotting_file = f"./Plotting_file_{fmin:.2f}_{fmax:.2f}_{depth}.0km.txt"

    df = pd.read_csv(Plotting_file, sep=' ', index_col=False)
    # df = df[(df['stlo_mean'] > extents[region][0]) & (df['stlo_mean'] < extents[region][1]) & (df['stla_mean'] > extents[region][2]) & (df['stla_mean'] < extents[region][3])]
    df = df.sort_values(by='multi')
    df = df.dropna()

    stlas = df['stla_mean'].to_numpy()
    stlos = df['stlo_mean'].to_numpy()

    evlas = df['evla'].to_numpy()
    evlos = df['evlo'].to_numpy()
    evdps = df['evdp'].to_numpy()

    mean_lo = np.mean([np.mean(stlos), np.mean(evlos)])
    mean_la = np.mean([np.mean(stlas), np.mean(evlas)])

    multis = df['multi'].to_numpy()
    baz_diffs = df['baz_diff'].to_numpy()
    baz_stds = df['baz_std_dev'].to_numpy()
    slow_diffs = df['slow_diff'].to_numpy()
    slow_stds = df['slow_std_dev'].to_numpy()

    vecs_x = df['del_x_slow'].to_numpy()
    vecs_y = df['del_y_slow'].to_numpy()
    mags = df['mag'].to_numpy()
    azs = df['az'].to_numpy()

    s_gcp_pierce_las = df['s_pierce_la'].to_numpy()
    s_gcp_pierce_los = df['s_pierce_lo'].to_numpy()

    s_reloc_pierce_las = df['s_reloc_pierce_la'].to_numpy()
    s_reloc_pierce_los = df['s_reloc_pierce_lo'].to_numpy()

    r_reloc_pierce_las = df['r_reloc_pierce_la'].to_numpy()
    r_reloc_pierce_los = df['r_reloc_pierce_lo'].to_numpy()

    r_gcp_pierce_las = df['r_pierce_la'].to_numpy()
    r_gcp_pierce_los = df['r_pierce_lo'].to_numpy()

    lats_lons_all_r_gcp = np.array(list(zip(np.deg2rad(r_gcp_pierce_las), np.deg2rad(r_gcp_pierce_los))))
    lats_lons_all_s_gcp = np.array(list(zip(np.deg2rad(s_gcp_pierce_las), np.deg2rad(s_gcp_pierce_los))))

    lats_lons_all_r_reloc = np.array(list(zip(np.deg2rad(r_reloc_pierce_las), np.deg2rad(r_reloc_pierce_los))))
    lats_lons_all_s_reloc = np.array(list(zip(np.deg2rad(s_reloc_pierce_las), np.deg2rad(s_reloc_pierce_los))))

    # indices where event depth is greater than the current depth
    id_ev = evdps > depth

    mag_variances_points = []
    vec_variances_points = []

    for i,points in enumerate([lats_lons_all_r_gcp, lats_lons_all_s_gcp, lats_lons_all_r_reloc, lats_lons_all_s_reloc]):
        tree = BallTree(
            points, leaf_size=points.shape[0] / 2, metric="haversine"
        )

        indices = tree.query_radius(points, bin_radians, count_only=False, return_distance = False)
        counts = tree.query_radius(points, bin_radians, count_only=True, return_distance = False)
        mags_var = []
        vecs_var = []
        for point_bin in indices:

            mags_neighbourhood = mags[point_bin]
            mags_variance = np.var(mags_neighbourhood)
            mags_var.append(mags_variance)

            slow_xs = vecs_x[point_bin]
            slow_ys = vecs_y[point_bin]

            # find mean x and y in bin
            mean_x = np.mean(slow_xs)
            mean_y = np.mean(slow_ys)

            diff_vecs = np.array(list(zip(slow_xs, slow_ys)))
            means = np.array([[mean_x,mean_y]])

            dists_from_mean = distance.cdist(diff_vecs, means, metric="euclidean")

            vect_variance = np.mean(np.power(dists_from_mean, 2))
            # print(diff_xs, diff_ys, dists_from_mean, vect_variance)
            vecs_var.append(vect_variance)



        mags_var = np.array(mags_var)
        vecs_var = np.array(vecs_var)
        counts_10 = np.where(counts < 10)[0]
        mags_var[counts_10] = -1
        vecs_var[counts_10] = -1

        if i == 1 or i == 3:
            # fig = plt.figure(figsize=(12,8))
            #
            # ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
            # ax.add_feature(cartopy.feature.LAND, edgecolor='black', color='lightgrey')
            #
            # plt.scatter(points[:,1], points[:,0], transform=ccrs.PlateCarree())
            # plt.show()
            mags_var[id_ev] = -1
            vecs_var[id_ev] = -1


        mag_variances_points.append(mags_var)
        vec_variances_points.append(vecs_var)

    df['mag_var_r_gcp'] = mag_variances_points[0]
    df['mag_var_s_gcp'] = mag_variances_points[1]
    df['mag_var_r_reloc'] = mag_variances_points[2]
    df['mag_var_s_reloc'] = mag_variances_points[3]

    df['vec_var_r_gcp'] = vec_variances_points[0]
    df['vec_var_s_gcp'] = vec_variances_points[1]
    df['vec_var_r_reloc'] = vec_variances_points[2]
    df['vec_var_s_reloc'] = vec_variances_points[3]

    out_file = f"./Plotting_file_{fmin:.2f}_{fmax:.2f}_{depth}.0km_variance.txt"

    df.to_csv(out_file)

    print(depth)
