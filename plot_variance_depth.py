#!/usr/bin/env python

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
from sklearn.neighbors import BallTree
import os
from scipy.stats import circvar

extents = {"US":[-140,-60,10,75], "EU":[-35,60,10,70], "SAm":[-120,-60,-70,-5], "SAs":[60,140,0,60], "Alaska":[-170,-140, 50, 75]}
region = "EU"
dir = f"figures_animations_variance_{region}"
try:
    os.mkdir(dir)
except:
    pass


depths = np.arange(500, 2900, 100)
depths = np.append(depths, 2891)

fmin = 0.20
fmax = 0.40
cmap = plt.cm.get_cmap('summer')


for i, depth in enumerate(depths):

    Plotting_file = f"./Plotting_file_{fmin:.2f}_{fmax:.2f}_{depth}.0km_variance.txt"

    df = pd.read_csv(Plotting_file, sep=',', index_col=False)

    df = df[(df['stlo_mean'] > extents[region][0]) & (df['stlo_mean'] < extents[region][1]) & (df['stla_mean'] > extents[region][2]) & (df['stla_mean'] < extents[region][3]) ]

    df = df.sort_values(by='multi')

    multis = df['multi'].to_numpy()
    multi_values = np.array(multis)
    multi_values[multi_values=='n'] = 0.0
    multi_values[multi_values=='y'] = 1.0
    multi_values[multi_values=='m'] = 1.0

    multi_values = multi_values.astype(float)

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

    var_s_gcp = df['var_s_gcp']
    var_r_gcp = df['var_r_gcp']
    var_s_reloc = df['var_s_reloc']
    var_r_reloc = df['var_r_reloc']


    fig = plt.figure(figsize=(12,8))

    ax = fig.add_subplot(111, projection=ccrs.PlateCarree(central_longitude=10))
    ax.add_feature(cartopy.feature.LAND, edgecolor='black', color='lightgrey')

    a = ax.scatter(
    r_reloc_pierce_los[var_r_reloc<0],
    r_reloc_pierce_las[var_r_reloc<0],
    transform=ccrs.Geodetic(),
    marker='o',
    s=80,
    linewidth=0.1,
    c='grey',
    edgecolors='black'
    )

    a = ax.scatter(
    r_reloc_pierce_los[var_r_reloc>0],
    r_reloc_pierce_las[var_r_reloc>0],
    transform=ccrs.Geodetic(),
    marker='o',
    s=80,
    linewidth=0.1,
    c=var_r_reloc[var_r_reloc>0],
    edgecolors='black',
    cmap='viridis'
    )

    plt.colorbar(a)

    ax.set_title(f"Depth {depth}")

    ax.set_xlabel('Longitude ($^{\circ}$)')
    ax.set_ylabel('Latitude ($^{\circ}$)')

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
