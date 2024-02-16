#!/usr/bin/env python

import cartopy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd
import numpy as np
import matplotlib.ticker as mticker
import os
import shutil
from sklearn.neighbors import BallTree
from scipy.stats import circvar

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

    Plotting_file = f"./Locus_Lines_{fmin}_{fmax}_{depth}.0.txt"

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

    s_gcp_pierce_las = df['s_pierce_la'].to_numpy()
    s_gcp_pierce_los = df['s_pierce_lo'].to_numpy()

    s_reloc_pierce_las = df['s_reloc_pierce_la'].to_numpy()
    s_reloc_pierce_los = df['s_reloc_pierce_lo'].to_numpy()

    r_reloc_pierce_las = df['r_reloc_pierce_la'].to_numpy()
    r_reloc_pierce_los = df['r_reloc_pierce_lo'].to_numpy()

    r_gcp_pierce_las = df['r_pierce_la'].to_numpy()
    r_gcp_pierce_los = df['r_pierce_lo'].to_numpy()

    Phis_1 = df['Phi_1']
    Phis_2 = df['Phi_2']

    lats_lons_all_r_gcp = np.array(list(zip(np.deg2rad(r_gcp_pierce_las), np.deg2rad(r_gcp_pierce_los))))
    lats_lons_all_s_gcp = np.array(list(zip(np.deg2rad(s_gcp_pierce_las), np.deg2rad(s_gcp_pierce_los))))

    lats_lons_all_r_reloc = np.array(list(zip(np.deg2rad(r_reloc_pierce_las), np.deg2rad(r_reloc_pierce_los))))
    lats_lons_all_s_reloc = np.array(list(zip(np.deg2rad(s_reloc_pierce_las), np.deg2rad(s_reloc_pierce_los))))


    # indices where event depth is greater than the current depth
    id_ev = evdps > depth

    locus_variances_points = []

    for i,points in enumerate([lats_lons_all_r_gcp, lats_lons_all_s_gcp, lats_lons_all_r_reloc, lats_lons_all_s_reloc]):
        tree = BallTree(
            points, leaf_size=points.shape[0] / 2, metric="haversine"
        )

        indices = tree.query_radius(points, bin_radians, count_only=False, return_distance = False)
        counts = tree.query_radius(points, bin_radians, count_only=True, return_distance = False)
        loci_var = []
        for point_bin in indices:

            phi_1_neighbourhood = Phis_1[point_bin]
            phi_2_neighbourhood = Phis_2[point_bin]

            min_phi_neighbourhood = np.minimum(phi_1_neighbourhood, phi_2_neighbourhood)
            loc_variance = circvar(np.radians(min_phi_neighbourhood))
            loci_var.append(loc_variance)

            # print("loci_var: ", loc_variance)



        loci_var = np.array(loci_var)
        counts_10 = np.where(counts < 10)[0]
        loci_var[counts_10] = np.nan

        if i == 1 or i == 3:
            # fig = plt.figure(figsize=(12,8))
            #
            # ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
            # ax.add_feature(cartopy.feature.LAND, edgecolor='black', color='lightgrey')
            #
            # plt.scatter(points[:,1], points[:,0], transform=ccrs.PlateCarree())
            # plt.show()
            loci_var[id_ev] = np.nan

        print(np.degrees(loci_var))
        locus_variances_points.append(loci_var)


    df['locus_var_r_gcp'] = locus_variances_points[0]
    df['locus_var_s_gcp'] = locus_variances_points[1]
    df['locus_var_r_reloc'] = locus_variances_points[2]
    df['locus_var_s_reloc'] = locus_variances_points[3]

    out_file = f"./Locus_file_{fmin:.2f}_{fmax:.2f}_{depth}.0_variance.txt"

    df.to_csv(out_file)

    print(depth)
