#!/usr/bin/env python

import cartopy
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')
import cartopy.crs as ccrs
import pandas as pd
import numpy as np
import matplotlib.ticker as mticker
from sklearn.cluster import KMeans
import os
import shutil
from sklearn.neighbors import BallTree
from scipy.spatial import distance
import netCDF4 as nc
from scipy.stats import circvar, circmean

import argparse


land_mask = cartopy.feature.NaturalEarthFeature(
    'physical',
    'land',
    scale='50m',
    edgecolor='black',
    facecolor='lightgrey',
    alpha=0.5
)


parser = argparse.ArgumentParser(description='Preprocessing script for data retrieved from obspy DMT')

parser.add_argument("-fmin","--minimuum_frequency", help="Enter the lowest frequency to be plotted.", type=float, required=True, action="store", default = 0.20)

parser.add_argument("-fmax","--maximuum_frequency", help="Enter the highest frequency to be plotted.", type=float, required=True, action="store", default = 0.40)

parser.add_argument("-r","--region", help="Enter the region you want to make the animation of (US/EU).", type=str, required=True, action="store", default='EU')

parser.add_argument("-v", "--verbose", help="Increase verbosity of output, just --verbose is enough to turn on.",action="store_true")

args = parser.parse_args()

fmin = args.minimuum_frequency
fmax = args.maximuum_frequency
region = args.region

plt.style.use('fivethirtyeight')



extents = {"US":[-140,-60,25,50], "EU":[-15,55,15,70], "SAm":[-120,-60,-70,-5], "SAs":[60,140,0,60], "Alaska":[-170,-140, 50, 75]}

extents_s = {"Asia":[100,160,-20,60], "SA":[-80,-20,-40,20], "Fiji":[160,-160,-40,0], 'SPacific':[160,-140,-40,15]}

side='r'

depths = np.arange(0, 2900, 100)
depths = np.append(depths, 2891)

# depths = [1800]
gird_lons = np.arange(-360,375, 15)
gird_lats = np.arange(-90,105, 15)

dir = f"vec_multi_maps{region}_{fmin}_{fmax}"

try:
    os.mkdir(dir)
except:
    pass
cmap = plt.cm.get_cmap('summer')

bin_radius = 200 #km


# tomo = "/Users/jamieward/PhD/Figures/Tomography_Models/GYPSUM_percent.nc"
# tomo = "/Users/jamieward/PhD/Figures/Tomography_Models/S362ANI_percent.nc"
# tomo = "/Users/jamieward/PhD/Figures/Tomography_Models/TX2019slab_percent.nc"
# tomo = "/Users/jamieward/PhD/Figures/Tomography_Models/SGLOBE-rani-voigt_percent.nc"
# tomo = "/Users/jamieward/PhD/Figures/Tomography_Models/csem-europe-2019.12.01.nc"
# tomo = "/Users/jamieward/PhD/Figures/Tomography_Models/csem-eastmed-2019.12.01.nc"
# tomo = "/Users/jamieward/PhD/Figures/Tomography_Models/csem-north-america-2019.12.01.nc"
# tomo = "/Users/jamieward/PhD/Figures/Tomography_Models/SEMum-NA14_kmps.nc"
# tomo = "/Users/jamieward/PhD/Figures/Tomography_Models/NA07_percent.nc"


# sst = nc.Dataset(tomo)
# depths_tomo = sst.variables['depth'][:]
# loc_depth = np.absolute(depths_tomo - depths[0]).argmin()
# print(depths_tomo[loc_depth])
# lon = sst.variables['longitude'][:]
# lat = sst.variables['latitude'][:]
# try:
#     dvs = sst.variables['dvs'][loc_depth,:,:]
# except:
#     try:
#         vsv = sst.variables['vsv'][loc_depth,:,:]
#         vsh = sst.variables['vsh'][loc_depth,:,:]
#
#         dvs = np.sqrt(np.power(vsh, 2) + np.power(vsv, 2))
#
#         dvs = ((dvs - dvs.mean()) / dvs.mean()) * 100
#
#     except:
#         dvs = sst.variables['Vs'][loc_depth,:,:]
#
# print(dvs)



# tomo = "/Users/jamieward/PhD/Figures/UCB_a3d_dist.SEMUCB-WM1.r20151019/SEMUCB_2891.xyz"
#
# tmp = np.loadtxt(tomo)
# lon = tmp[:,1].reshape((91,181))
# lat = tmp[:,2].reshape((91,181))
# dvs = tmp[:,3].reshape((91,181))

# tomo = f"/Users/jamieward/PhD/Figures/SP12RTS_plotting_GMT5/S40_2891km.nc"
# tomo = '/Users/jamieward/PhD/Figures/Tomography_Models/EU60_500km.nc'

# sst = nc.Dataset(tomo)
# lon = sst.variables['x'][:]
# lat = sst.variables['y'][:]
# dvs = sst.variables['z'][:,:]
# #
# lon2d, lat2d = np.meshgrid(lat, lon)


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
variance_limit = variance_estimate - (2*variance_std)



for d,depth in enumerate(depths):

    Plotting_file = f"./Plotting_file_{fmin:.2f}_{fmax:.2f}_{depth}.0km_variance.txt"
    Locus_file = f"Locus_Lines_{fmin}_{fmax}_{depth}.0.txt"
    df_locus = pd.read_csv(Locus_file, sep=' ', index_col=False)
    df_locus = df_locus[(df_locus['stlo_mean'] > extents[region][0]) & (df_locus['stlo_mean'] < extents[region][1]) & (df_locus['stla_mean'] > extents[region][2]) & (df_locus['stla_mean'] < extents[region][3])]
    df_locus = df_locus.dropna()

    df = pd.read_csv(Plotting_file, sep=',', index_col=False)
    df = df[(df['stlo_mean'] > extents[region][0]) & (df['stlo_mean'] < extents[region][1]) & (df['stla_mean'] > extents[region][2]) & (df['stla_mean'] < extents[region][3])]
    # df = df.sort_values(by='multi')
    df = df.dropna()

    # df = df[(df['baz_pred'] >= 180) & (df['baz_pred'] <= 270)]
    # df_locus = df_locus[(df_locus['baz_pred'] >= 180) & (df_locus['baz_pred'] <= 270)]


    df_multi = df.loc[df['multi'] == 'y']
    df_no_multi = df.loc[df['multi'] == 'n']

    if side == 'r':
        df_var = df[df['vec_var_r_reloc'] <= variance_limit]
    elif side == 's':
        df_var = df[df['vec_var_s_reloc'] <= variance_limit]


    radius = 6371 - depth
    bin_radians = bin_radius / radius

    bin_deg = np.degrees(bin_radians)

    spacing = bin_deg/2
    min_lo = -180 + bin_deg
    max_lo = 180 - bin_deg
    n_lo = (max_lo - min_lo) / spacing

    min_la = -90 + bin_deg
    max_la = 90 - bin_deg
    n_la = (max_la - min_la) / spacing

    los_centre_grid = np.linspace(min_lo, max_lo, int(n_lo+1))
    las_centre_grid = np.linspace(min_la, max_la, int(n_la+1))
    los,las = np.meshgrid(los_centre_grid, las_centre_grid)
    los = los.flatten()
    las = las.flatten()

    print(los)
    print(las)


    la_lo_grid = np.array(list(zip(np.deg2rad(las), np.deg2rad(los))))

    stlas = df['stla_mean'].to_numpy()
    stlos = df['stlo_mean'].to_numpy()

    evlas = df['evla'].to_numpy()
    evlos = df['evlo'].to_numpy()

    mean_lo = np.mean([np.mean(stlos), np.mean(evlos)])
    mean_la = np.mean([np.mean(stlas), np.mean(evlas)])

    multis = df['multi'].to_numpy()
    baz_diffs = df['baz_diff'].to_numpy()
    baz_stds = df['baz_std_dev'].to_numpy()
    slow_diffs = df['slow_diff'].to_numpy()
    slow_stds = df['slow_std_dev'].to_numpy()

    vecs_x = df_var['del_x_slow'].to_numpy()
    vecs_y = df_var['del_y_slow'].to_numpy()
    mags = df_var['mag'].to_numpy()
    azs = df_var['az'].to_numpy()

    s_gcp_pierce_las = df['s_pierce_la'].to_numpy()
    s_gcp_pierce_los = df['s_pierce_lo'].to_numpy()

    s_reloc_pierce_las = df_var['s_reloc_pierce_la'].to_numpy()
    s_reloc_pierce_los = df_var['s_reloc_pierce_lo'].to_numpy()

    r_reloc_pierce_las = df_var['r_reloc_pierce_la'].to_numpy()
    r_reloc_pierce_los = df_var['r_reloc_pierce_lo'].to_numpy()

    r_gcp_vec_pierce_las = df_var['r_pierce_la'].to_numpy()
    s_gcp_vec_pierce_las = df_var['s_pierce_la'].to_numpy()

    r_gcp_vec_pierce_los = df_var['r_pierce_lo'].to_numpy()
    s_gcp_vec_pierce_los = df_var['s_pierce_lo'].to_numpy()

    r_gcp_pierce_las = df['r_pierce_la'].to_numpy()
    r_gcp_pierce_los = df['r_pierce_lo'].to_numpy()

    # process multipathing data
    multi_values = np.array(multis)
    multi_values[multi_values=='n'] = 0.0
    multi_values[multi_values=='y'] = 1.0
    multi_values[multi_values=='m'] = 1.0

    multi_values = multi_values.astype(float)

    # cluster magnitudes
    # k=2
    # kmeans = KMeans(n_clusters=k)
    # df['mag_cluster'] = kmeans.fit_predict(df[['mag']])
    #
    # mag_clusters = df['mag_cluster'].values.astype(float)

    # create locus points for plotting
    Phi_1 = df_locus['Phi_1'].to_numpy()
    Phi_2 = df_locus['Phi_2'].to_numpy()

    s_gcp_pierce_las_locus = df_locus['s_pierce_la'].to_numpy()
    s_gcp_pierce_los_locus = df_locus['s_pierce_lo'].to_numpy()

    s_reloc_pierce_las_locus = df_locus['s_reloc_pierce_la'].to_numpy()
    s_reloc_pierce_los_locus = df_locus['s_reloc_pierce_lo'].to_numpy()

    r_reloc_pierce_las_locus = df_locus['r_reloc_pierce_la'].to_numpy()
    r_reloc_pierce_los_locus = df_locus['r_reloc_pierce_lo'].to_numpy()

    r_gcp_pierce_las_locus = df_locus['r_pierce_la'].to_numpy()
    r_gcp_pierce_los_locus = df_locus['r_pierce_lo'].to_numpy()


    dist = 1.5 # degrees
    dys = np.cos(np.deg2rad(Phi_1))*dist
    dxs = np.sin(np.deg2rad(Phi_1))*dist


    ## ----- bin the vectors and proportion multi----- ##

    if side == 'r':
        lats_lons_all = np.array(list(zip(np.deg2rad(r_reloc_pierce_las), np.deg2rad(r_reloc_pierce_los))))
        lats_lons_gcp = np.array(list(zip(np.deg2rad(r_gcp_vec_pierce_las), np.deg2rad(r_gcp_vec_pierce_los))))
        lats_lons_locus = np.array(list(zip(np.deg2rad(r_gcp_pierce_las_locus), np.deg2rad(r_gcp_pierce_los_locus))))

    elif side == 's':
        lats_lons_all = np.array(list(zip(np.deg2rad(s_reloc_pierce_las), np.deg2rad(s_reloc_pierce_los))))
        lats_lons_gcp = np.array(list(zip(np.deg2rad(s_gcp_vec_pierce_las), np.deg2rad(s_gcp_vec_pierce_las))))
        lats_lons_locus = np.array(list(zip(np.deg2rad(s_gcp_pierce_las_locus), np.deg2rad(s_gcp_pierce_los_locus))))

    else:
        print('side needs to be "s" or "r"')
        exit()

    tree_all = BallTree(
        lats_lons_all, leaf_size=lats_lons_all.shape[0] / 2, metric="haversine"
    )

    tree_gcp = BallTree(
        lats_lons_gcp, leaf_size=lats_lons_gcp.shape[0] / 2, metric="haversine"
    )

    tree_locus = BallTree(
        lats_lons_locus, leaf_size=lats_lons_locus.shape[0] / 2, metric="haversine"
    )

    indices = tree_all.query_radius(la_lo_grid, np.radians(bin_deg), count_only=False, return_distance = False)
    counts_all = tree_all.query_radius(la_lo_grid, np.radians(bin_deg), count_only=True, return_distance = False)

    indices_gcp = tree_gcp.query_radius(la_lo_grid, np.radians(bin_deg), count_only=False, return_distance = False)
    counts_gcp = tree_gcp.query_radius(la_lo_grid, np.radians(bin_deg), count_only=True, return_distance = False)

    indices_locus = tree_locus.query_radius(la_lo_grid, np.radians(bin_deg), count_only=False, return_distance = False)
    counts_locus = tree_locus.query_radius(la_lo_grid, np.radians(bin_deg), count_only=True, return_distance = False)



    counts_10 = np.where(counts_all >= 10)
    indices_filt = indices[counts_10]
    counts_10_gcp = np.where(counts_gcp >= 10)
    indices_filt_gcp = indices_gcp[counts_10_gcp]

    counts_10_locus = np.where(counts_locus >= 10)
    indices_filt_locus = indices_locus[counts_10_locus]

    los_filt = los[counts_10]
    las_filt = las[counts_10]

    los_filt_gcp = los[counts_10_gcp]
    las_filt_gcp = las[counts_10_gcp]

    los_filt_locus = los[counts_10_locus]
    las_filt_locus = las[counts_10_locus]

    multi_proportions = []
    n_multi = []
    mean_mags = []
    mag_cluster_means = []
    mag_stds = []
    slow_xs_mean = []
    slow_ys_mean = []
    slow_xs_max = []
    slow_ys_max = []

    slow_xs_mean_gcp = []
    slow_ys_mean_gcp = []


    az_vars = []
    az_means = []
    mean_slows = []

    Phis = []
    dx_bin = []
    dy_bin = []

    for point_bin_locus in indices_filt_locus:
        Phi_mean = circmean(np.radians(Phi_1[point_bin_locus]))
        Phis.append(Phi_mean)
        dy = np.cos(Phi_mean)*dist
        dx = np.sin(Phi_mean)*dist

        dx_bin.append(dx)
        dy_bin.append(dy)


    for point_bin_gcp in indices_filt_gcp:
        prop_multi = np.sum(multi_values[point_bin_gcp])/len(point_bin_gcp)
        multi_proportions.append(prop_multi)
        n_multi.append(np.sum(multi_values[point_bin_gcp]))
        print(multi_values[point_bin_gcp], np.sum(multi_values[point_bin_gcp]))

        slow_xs_mean_gcp.append(np.mean(vecs_x[point_bin_gcp]))
        slow_ys_mean_gcp.append(np.mean(vecs_y[point_bin_gcp]))

    for point_bin in indices_filt:

        mean_mag = np.mean(mags[point_bin])
        mag_std = np.std(mags[point_bin])

        slow_xs_max.append(np.sum(vecs_x[point_bin]))
        slow_ys_max.append(np.sum(vecs_y[point_bin]))

        slow_xs_mean.append(np.mean(vecs_x[point_bin]))
        slow_ys_mean.append(np.mean(vecs_y[point_bin]))
        az_variance = circvar(np.radians(azs[point_bin]))
        az_variance = np.degrees(az_variance)
        az_mean = circmean(np.radians(azs[point_bin]))

        mean_mags.append(mean_mag)
        mag_stds.append(mag_std)
        az_vars.append(az_variance)
        az_means.append(az_mean)
        mean_slow = np.mean(slow_diffs[point_bin])
        mean_slows.append(mean_slow)

    slow_xs_mean = np.array(slow_xs_mean)
    slow_ys_mean = np.array(slow_ys_mean)
    mean_slows = np.array(mean_slows)
    slow_colours = np.copy(mean_slows)
    slow_colours[slow_colours>0] = 1
    slow_colours[slow_colours<0] = 0
    slow_colours = slow_colours.astype(int)
    slow_colours = slow_colours.astype(str)
    slow_colours[slow_colours == '0'] = 'blue'
    slow_colours[slow_colours == '1'] = 'red'
    print(slow_colours)


    levels = [-3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0.5, 1., 1.5, 2., 2.5, 3., 3.5]
    colours = ['darkred', 'darkred', 'darkred', 'darkred', 'darkred', 'darkred',
               'darkred', 'darkblue', 'darkblue', 'darkblue', 'darkblue',
               'darkblue', 'darkblue', 'darkblue']

    # make plot

    fig1 = plt.figure(figsize=(12,8), facecolor='none')

    ax1 = fig1.add_subplot(111, projection=ccrs.PlateCarree())

    # t = ax1.contourf(lon, lat, dvs, 50,
    #              transform=ccrs.PlateCarree(), cmap='seismic_r',
    #              vmin=-4, vmax=4, zorder=1, alpha=0.8, levels=np.linspace(-4,4,25))

    # t = ax1.contourf(lon, lat, dvs, 50,
    #              transform=ccrs.PlateCarree(), cmap='seismic_r',
    #              zorder=1, alpha=0.8)


    # c = ax2.contourf(lon, lat, dvs, 50,
    #              transform=ccrs.PlateCarree(), cmap='seismic_r',
    #              zorder=1, alpha=0.8, vmin=-4, vmax=4,levels=np.linspace(-4,4,25))




    # plt.colorbar(t, label="$\delta V_{s}$ (%)")
    # ax1.add_feature(cartopy.feature.LAND, facecolor='lightgrey')
    # ax1.add_feature(land_mask)
    # ax1.add_feature(cartopy.feature.OCEAN  , facecolor='lightgrey')

    ax1.coastlines(zorder=2, resolution='50m', color='black', linewidth=1, alpha=0.5)



    if side == 's':
        ax1.set_extent(extents_s['SPacific'])
        pass
    else:
        ax1.set_extent([los_filt.min() - 5, los_filt.max() + 5, las_filt.min() - 5, las_filt.max() + 5])


    # for i, dx in enumerate(dxs):
    #     dy = dys[i]
    #     if side == 'r':
    #         pierce_la = r_gcp_pierce_las_locus[i]
    #         pierce_lo = r_gcp_pierce_los_locus[i]
    #
    #     elif side == 's':
    #         pierce_la = s_gcp_pierce_las_locus[i]
    #         pierce_lo = s_gcp_pierce_los_locus[i]
    #
    #     p1_y = pierce_la + dy
    #     p2_y = pierce_la - dy
    #
    #     p1_x = pierce_lo + dx
    #     p2_x = pierce_lo - dx
    #
    #     ax1.plot(
    #         [p1_x, p2_x],
    #         [p1_y, p2_y],
    #         transform=ccrs.Geodetic(),
    #         c='black',
    #         zorder=2,
    #         label="Multi",
    #         linewidth=2,
    #     )
    # if side == 'r':
    #     ax1.scatter(
    #         r_gcp_pierce_los_locus,
    #         r_gcp_pierce_las_locus,
    #         transform=ccrs.PlateCarree(),
    #         marker="o",
    #         c='grey',
    #         zorder=3,
    #         s=80,
    #         label="Multi",
    #         ec='black',
    #         linewidth=0.5,
    #     )
    # elif side == 's':
    #     ax1.scatter(
    #         s_gcp_pierce_los_locus,
    #         s_gcp_pierce_las_locus,
    #         transform=ccrs.PlateCarree(),
    #         marker="o",
    #         c='grey',
    #         zorder=3,
    #         s=80,
    #         label="Multi",
    #         ec='black',
    #         linewidth=0.5,
    #     )


    m = ax1.scatter(
        los_filt_gcp,
        las_filt_gcp,
        transform=ccrs.PlateCarree(),
        marker="o",
        c=multi_proportions,
        zorder=3,
        s=100,
        label="Multi",
        ec='black',
        linewidth=1,
        vmin=0,
        vmax=0.8,
        cmap = 'viridis'
    )

    # m = ax1.tricontourf(los_filt_gcp,las_filt_gcp,multi_proportions, 20)


    for i, dx in enumerate(dx_bin):
        dy = dy_bin[i]

        pierce_la = las_filt_locus[i]
        pierce_lo = los_filt_locus[i]


        p1_y = pierce_la + dy
        p2_y = pierce_la - dy

        p1_x = pierce_lo + dx
        p2_x = pierce_lo - dx

        # ax1.plot(
        #     [p1_x, p2_x],
        #     [p1_y, p2_y],
        #     transform=ccrs.Geodetic(),
        #     c='white',
        #     zorder=3,
        #     label="Multi",
        #     linewidth=2,
        #     alpha=0.8
        # )

        # ax12.scatter(
        #     pierce_lo,
        #     pierce_la,
        #     transform=ccrs.Geodetic(),
        #     c='grey',
        #     zorder=3,
        #     label="Multi",
        #     s=80,
        #     ec='black'
        # )


    plt.colorbar(m,  label = 'Multipathing proportion', orientation = 'horizontal')



    ax1.set_title(f"Multipathing Proportions")
    ax1.set_xlabel('Longitude ($^{\circ}$)')
    ax1.set_ylabel('Latitude ($^{\circ}$)')


    # ax.set_extent(extents[region])
    gl1 = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

    gl1.xlines = True
    gl1.ylines = True
    gl1.xlabels_top = False
    gl1.ylabels_right = False
    gl1.xlabels_bottom = True
    gl1.ylabels_left = True
    #
    gl1.xlocator = mticker.FixedLocator(gird_lons)
    gl1.ylocator = mticker.FixedLocator(gird_lats)

    plt.savefig(f"./multi_prop_{depth}_{fmin}_{fmax}.pdf")
    os.rename(f"./multi_prop_{depth}_{fmin}_{fmax}.pdf", f"{dir}/multi_prop_{depth}_{fmin}_{fmax}.pdf")


    fig2 = plt.figure(figsize=(12,8), facecolor='none')
    ax2 = fig2.add_subplot(111, projection=ccrs.PlateCarree())

    # c = ax2.contourf(lon, lat, dvs, 50,
    #              transform=ccrs.PlateCarree(), cmap='seismic_r',
    #              vmin=-4, vmax=4, zorder=1, alpha=0.8, levels=np.linspace(-4,4,25))

    # c = ax2.contourf(lon, lat, dvs, 50,
    #              transform=ccrs.PlateCarree(), cmap='seismic_r',
    #              zorder=1, alpha=0.8, vmin=-4, vmax=4,levels=np.linspace(-4,4,25))

    # plt.colorbar(c, label="$\delta V_{s}$ (%)")
    # ax2.add_feature(land_mask)
    # ax1.add_feature(cartopy.feature.OCEAN  , facecolor='lightgrey')

    ax2.coastlines(zorder=2, resolution='50m', color='black', linewidth=1, alpha=0.5)

    if side == 's':
        ax2.set_extent(extents_s['SPacific'])
        pass
    else:
        ax2.set_extent([los_filt.min() - 5, los_filt.max() + 5, las_filt.min() - 5, las_filt.max() + 5])

    q = ax2.quiver(
    los_filt,
    las_filt,
    np.array(slow_xs_mean),
    np.array(slow_ys_mean),
    # mean_mags,
    zorder=3,
    label='Binned slowness vector residuals',
    transform=ccrs.PlateCarree(),
    fc = 'orange',
    ec = 'black',
    lw=0.75,
    scale=4,
    width=0.0005,
    headwidth=20,
    headlength = 20,
    headaxislength=15,
    minshaft=1,
    minlength=0.5,
    scale_units='inches'
    # cmap='hot_r'
    )

    plt.quiverkey(q,
    X=0.04,
    Y=1.02,
    U=1,
    labelpos='E',
    label='1 s/$^{\circ}$')

    ax2.set_title(f"Binned Vectors")
    ax2.set_xlabel('Longitude ($^{\circ}$)')
    ax2.set_ylabel('Latitude ($^{\circ}$)')

    # ax.set_extent(extents[region])
    gl2 = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

    gl2.xlines = True
    gl2.ylines = True
    gl2.xlabels_top = False
    gl2.ylabels_right = False
    gl2.xlabels_bottom = True
    gl2.ylabels_left = True
    #
    gl2.xlocator = mticker.FixedLocator(gird_lons)
    gl2.ylocator = mticker.FixedLocator(gird_lats)

    plt.savefig(f"./vector_bin_{depth}_{fmin}_{fmax}.pdf")
    os.rename(f"./vector_bin_{depth}_{fmin}_{fmax}.pdf", f"{dir}/vector_bin_{depth}_{fmin}_{fmax}.pdf")



    # n_numbers = len(Phi_1)
    # bins_number = 36  # the [0, 360) interval will be subdivided into this
    # # number of equal bins
    # bins = np.linspace(0.0, 2 * np.pi, bins_number + 1)
    # n, _, _ = plt.hist(np.radians(Phi_2), bins)

    # width = 2 * np.pi / bins_number

    # fig3 = plt.figure(figsize=(10,10))
    # ax3 = fig3.add_subplot(111, projection='polar')
    # bars = ax3.bar(bins[:bins_number], n, width=width, bottom=0.0)
    # for bar in bars:
    #     bar.set_alpha(0.5)
    # plt.show()


    fig4 = plt.figure(figsize=(12,8), facecolor='none')

    ax4 = fig4.add_subplot(111, projection=ccrs.PlateCarree())

    # t = ax4.contourf(lon, lat, dvs, 50,
    #              transform=ccrs.PlateCarree(), cmap='seismic_r',
    #              vmin=-4, vmax=4, zorder=1, alpha=0.8, levels=np.linspace(-4,4,25))

    # t = ax4.contourf(lon, lat, dvs, 50,
    #              transform=ccrs.PlateCarree(), cmap='seismic_r',
    #              zorder=1, alpha=0.8)


    # plt.colorbar(t, label="$\delta V_{s}$ (%)")
    # ax4.add_feature(cartopy.feature.OCEAN,facecolor=("grey"))
    # ax4.add_feature(cartopy.feature.LAND,facecolor=("grey"))
    # ax4.add_feature(land_mask)
    # ax4.add_feature(cartopy.feature.OCEAN  , facecolor='lightgrey')

    ax4.coastlines(zorder=2, resolution='50m', color='black', linewidth=1, alpha=0.5)

    if side == 's':
        ax4.set_extent(extents_s['SPacific'])
        pass

    else:
        ax4.set_extent([los_filt.min() - 5, los_filt.max() + 5, las_filt.min() - 5, las_filt.max() + 5])


    # m = ax4.scatter(
    #     los_filt_gcp,
    #     las_filt_gcp,
    #     transform=ccrs.PlateCarree(),
    #     marker="o",
    #     c=multi_proportions,
    #     zorder=3,
    #     s=80,
    #     label="Multi",
    #     ec='black',
    #     linewidth=0.5,
    #     vmin=0,
    #     vmax=1
    # )


    for i, dx in enumerate(dx_bin):
        dy = dy_bin[i]

        pierce_la = las_filt_locus[i]
        pierce_lo = los_filt_locus[i]


        p1_y = pierce_la + dy
        p2_y = pierce_la - dy

        p1_x = pierce_lo + dx
        p2_x = pierce_lo - dx

        ax4.plot(
            [p1_x, p2_x],
            [p1_y, p2_y],
            transform=ccrs.Geodetic(),
            c='black',
            zorder=3,
            label="Multi",
            linewidth=2,
        )

        # ax4.scatter(
        #     pierce_lo,
        #     pierce_la,
        #     transform=ccrs.Geodetic(),
        #     c='none',
        #     zorder=3,
        #     label="Multi",
        #     s=80,
        #     ec='black'
        # )


    # plt.colorbar(m,  label = 'Multipathing proportion', orientation = 'horizontal')



    ax4.set_title(f"Mean Loci")
    ax4.set_xlabel('Longitude ($^{\circ}$)')
    ax4.set_ylabel('Latitude ($^{\circ}$)')


    # ax.set_extent(extents[region])
    gl4 = ax4.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

    gl4.xlines = True
    gl4.ylines = True
    gl4.xlabels_top = False
    gl4.ylabels_right = False
    gl4.xlabels_bottom = True
    gl4.ylabels_left = True
    #
    gl4.xlocator = mticker.FixedLocator(gird_lons)
    gl4.ylocator = mticker.FixedLocator(gird_lats)

    plt.savefig(f"./multi_locus_{depth}_{fmin}_{fmax}.pdf")
    os.rename(f"./multi_locus_{depth}_{fmin}_{fmax}.pdf", f"{dir}/multi_locus_{depth}_{fmin}_{fmax}.pdf")

    plt.show()
    plt.close()
