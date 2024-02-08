#!/usr/bin/env python

"""
Python code to bin the vector field and calculate the divergence and convergence
spatially. Then this is plotted using cartopy and matplotlib.
"""

import cartopy
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import pandas as pd
import numpy as np
import matplotlib.ticker as mticker
import matplotlib as mpl
from sklearn.cluster import KMeans
import os
import shutil
from sklearn.neighbors import BallTree
from scipy.spatial import distance
import netCDF4 as nc
from scipy.stats import circvar, circmean

# for interpolating tomography
from scipy.interpolate import interp1d, griddata

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

parser.add_argument("-d", "--depths", help="Enter depths for plotting", nargs='+', type=float, required=True)

parser.add_argument("-v", "--verbose", help="Increase verbosity of output, just --verbose is enough to turn on.",action="store_true")

args = parser.parse_args()

fmin = args.minimuum_frequency
fmax = args.maximuum_frequency
region = args.region
depths = args.depths

plt.style.use('ggplot')

# to find consecutive zeros
def zero_runs(a):
    # Create an array that is 1 where a is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(a, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

extents = {"US":[-140,-60,25,50], "EU":[-15,55,15,70], "SAm":[-120,-60,-70,-5], "SAs":[60,140,0,60], "Alaska":[-170,-140, 50, 75]}

extents_s = {"Asia":[100,160,-20,60], "SA":[-80,-20,-40,20], "Fiji":[160,-160,-40,0], 'SPacific':[160,-140,-40,15]}

side='r'

hotspot_file_courtillot="/Users/jamieward/PhD/Figures/Courtillot_2003_plumes_table.txt"
hotspot_file_montelli="/Users/jamieward/PhD/Figures/Montelli_2006_plumes_table.txt"

hotspot_locations_c = pd.read_csv(hotspot_file_courtillot, sep=',', index_col=None)
hotspot_lons_c = hotspot_locations_c['Lon'].to_numpy().astype(float)
hotspot_lats_c = hotspot_locations_c['Lat'].to_numpy().astype(float)

hotspot_locations_m = pd.read_csv(hotspot_file_montelli, sep=',', index_col=None)
hotspot_lons_m = hotspot_locations_m['Longitude'].to_numpy().astype(float)
hotspot_lats_m = hotspot_locations_m['Latitude'].to_numpy().astype(float)

gird_lons = np.arange(-360,375, 15)
gird_lats = np.arange(-90,105, 15)

dir = f"vec_multi_maps{region}_{fmin}_{fmax}"

try:
    os.mkdir(dir)
except:
    pass
cmap = plt.cm.get_cmap('summer')

bin_radius = 200 #km
min_count = 10

# number of bins to calculate divergence over
divergence_steps = 1


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

variance_limit = variance_estimate - (2*variance_std)

tomo='/Users/jamieward/PhD/Figures/Tomography_Models/NA13.r0.0.nc'
sst = nc.Dataset(tomo)
lon = sst.variables['longitude'][:].data
lat = sst.variables['latitude'][:].data
dvs = sst.variables['dvs'][:].data
dvs_new = np.where(dvs > 100, 0, dvs)
dvs_new = np.where(dvs_new < -100, 0, dvs_new)
depths_tomo = sst.variables['depth'][:].data

tomog_interp_vs = interp1d(depths_tomo, dvs_new, axis=0, bounds_error=False, fill_value='extrapolate')

for d,depth in enumerate(depths):

    if depth <= 1200:
        tomo_model='NA13'
        dvs_depth = tomog_interp_vs(depth)
        

    else:
        tomo_model='SEMUCB'
        tomo = f"/Users/jamieward/PhD/Figures/UCB_a3d_dist.SEMUCB-WM1.r20151019/SEMUCB_{depth}.xyz"
        tmp = np.loadtxt(tomo)
        lon = tmp[:,1].reshape((91,181))
        lat = tmp[:,2].reshape((91,181))
        dvs_depth = tmp[:,3].reshape((91,181))

    Plotting_file = f"./Plotting_file_{fmin:.2f}_{fmax:.2f}_{depth}km_variance.txt"

    df = pd.read_csv(Plotting_file, sep=',', index_col=False)
    df = df[(df['stlo_mean'] > extents[region][0]) & (df['stlo_mean'] < extents[region][1]) & (df['stla_mean'] > extents[region][2]) & (df['stla_mean'] < extents[region][3])]
    # df = df.sort_values(by='multi')
    df = df.dropna()

    # df = df[(df['baz_pred'] >= 180) & (df['baz_pred'] <= 270)]
    # df_locus = df_locus[(df_locus['baz_pred'] >= 180) & (df_locus['baz_pred'] <= 270)]

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
    n_lo = int((max_lo - min_lo) / spacing)

    min_la = -90 + bin_deg
    max_la = 90 - bin_deg
    n_la = int((max_la - min_la) / spacing)

    los_centre_grid = np.linspace(min_lo, max_lo, int(n_lo+1))
    las_centre_grid = np.linspace(min_la, max_la, int(n_la+1))
    los_grid,las_grid = np.meshgrid(los_centre_grid, las_centre_grid)
    los = los_grid.flatten()
    las = las_grid.flatten()


    names = df['Name'].to_numpy()
    la_lo_grid = np.array(list(zip(np.deg2rad(las), np.deg2rad(los))))

    stlas = df['stla_mean'].to_numpy()
    stlos = df['stlo_mean'].to_numpy()

    evlas = df['evla'].to_numpy()
    evlos = df['evlo'].to_numpy()

    mean_lo = np.mean([np.mean(stlos), np.mean(evlos)])
    mean_la = np.mean([np.mean(stlas), np.mean(evlas)])

    vecs_x = df_var['del_x_slow'].to_numpy()
    vecs_y = df_var['del_y_slow'].to_numpy()
    mags = df_var['mag'].to_numpy()
    azs = df_var['az'].to_numpy()
    multis = df['multi'].to_numpy()

    s_gcp_pierce_las = df['s_pierce_la'].to_numpy()
    s_gcp_pierce_los = df['s_pierce_lo'].to_numpy()

    s_reloc_pierce_las = df_var['s_reloc_pierce_la'].to_numpy()
    s_reloc_pierce_los = df_var['s_reloc_pierce_lo'].to_numpy()

    r_reloc_pierce_las = df_var['r_reloc_pierce_la'].to_numpy()
    r_reloc_pierce_los = df_var['r_reloc_pierce_lo'].to_numpy()

    r_gcp_vec_pierce_las = df['r_pierce_la'].to_numpy()
    s_gcp_vec_pierce_las = df['s_pierce_la'].to_numpy()

    r_gcp_vec_pierce_los = df['r_pierce_lo'].to_numpy()
    s_gcp_vec_pierce_los = df['s_pierce_lo'].to_numpy()

    r_gcp_pierce_las = df['r_pierce_la'].to_numpy()
    r_gcp_pierce_los = df['r_pierce_lo'].to_numpy()

    # process multipathing data

    name_numbers = []
    for name in names:
        name_numbers.append(name.split('_')[2])

    name_numbers = np.array(name_numbers)
    multi_count = np.count_nonzero(name_numbers.astype(float) == 2)
    non_multi_count = np.count_nonzero(name_numbers.astype(float) == 1)

    multi_values = np.array(multis)
    multi_values[multi_values=='n'] = 0.0
    multi_values[multi_values=='y'] = 1.0
    multi_values[multi_values=='m'] = 1.0

    multi_values = multi_values.astype(float)


    ## ----- bin the vectors and proportion multi----- ##

    if side == 'r':
        lats_lons_all = np.array(list(zip(np.deg2rad(r_reloc_pierce_las), np.deg2rad(r_reloc_pierce_los))))
        lats_lons_gcp = np.array(list(zip(np.deg2rad(r_gcp_vec_pierce_las), np.deg2rad(r_gcp_vec_pierce_los))))

    elif side == 's':
        lats_lons_all = np.array(list(zip(np.deg2rad(s_reloc_pierce_las), np.deg2rad(s_reloc_pierce_los))))
        lats_lons_gcp = np.array(list(zip(np.deg2rad(s_gcp_vec_pierce_las), np.deg2rad(s_gcp_vec_pierce_las))))

    else:
        print('side needs to be "s" or "r"')
        exit()

    tree_all = BallTree(
        lats_lons_all, leaf_size=lats_lons_all.shape[0] / 2, metric="haversine"
    )

    tree_gcp = BallTree(
        lats_lons_gcp, leaf_size=lats_lons_gcp.shape[0] / 2, metric="haversine"
    )

    indices = tree_all.query_radius(la_lo_grid, np.radians(bin_deg), count_only=False, return_distance = False)
    counts_all = tree_all.query_radius(la_lo_grid, np.radians(bin_deg), count_only=True, return_distance = False)

    indices_gcp = tree_gcp.query_radius(la_lo_grid, np.radians(bin_deg), count_only=False, return_distance = False)
    counts_gcp = tree_gcp.query_radius(la_lo_grid, np.radians(bin_deg), count_only=True, return_distance = False)

    counts_10 = np.where(counts_all >= min_count)
    counts_less_10 = np.where(counts_all < min_count)
    indices_filt = indices[counts_10]

    counts_10_gcp = np.where(counts_gcp >= min_count)
    counts_less_10_gcp = np.where(counts_gcp < min_count)
    indices_filt_gcp = indices_gcp[counts_10_gcp]

    los_filt = los[counts_10]
    las_filt = las[counts_10]

    los_nans = np.copy(los)
    las_nans = np.copy(las)

    los_islands = np.zeros(los.shape)
    las_islands = np.zeros(las.shape)

    los_islands[counts_less_10] = 1
    las_islands[counts_less_10] = 1

    los_islands[counts_10] = 0
    las_islands[counts_10] = 0

    los_filt_gcp = los[counts_10_gcp]
    las_filt_gcp = las[counts_10_gcp]

    slow_xs_max = []
    slow_ys_max = []

    slow_xs_mean_vec = []
    slow_ys_mean_vec = []

    slow_xs_mean_gcp = []
    slow_ys_mean_gcp = []


    n_multi = []


    slow_xs_mean = np.zeros(los.shape)
    slow_ys_mean = np.zeros(los.shape)
    multi_proportions = np.zeros(los.shape)

    for ind, count in enumerate(counts_all):
        if count < min_count:
            slow_xs_mean[ind] = np.nan
            slow_ys_mean[ind] = np.nan

        elif count >= min_count:
            point_bin = indices[ind]
            slow_xs_mean[ind] = np.mean(vecs_x[point_bin])
            slow_ys_mean[ind] = np.mean(vecs_y[point_bin])
            slow_xs_mean_vec.append(np.mean(vecs_x[point_bin]))
            slow_ys_mean_vec.append(np.mean(vecs_y[point_bin]))


    for ind_gcp, count_gcp in enumerate(counts_gcp):
        if count_gcp < min_count:
            multi_proportions[ind_gcp] = np.nan

        elif count_gcp >= min_count:
            point_bin_gcp = indices_gcp[ind_gcp]
            prop_multi = np.sum(multi_values[point_bin_gcp])/len(point_bin_gcp)
            multi_proportions[ind_gcp] = prop_multi



    # for point_bin_gcp in indices_filt_gcp:
    #     prop_multi = np.sum(multi_values[point_bin_gcp])/len(point_bin_gcp)
    #     multi_proportions.append(prop_multi)
    #     n_multi.append(np.sum(multi_values[point_bin_gcp]))
    #     print(multi_values[point_bin_gcp], np.sum(multi_values[point_bin_gcp]))
    #

    # for point_bin in indices_filt:
    #
    #     slow_xs_max.append(np.sum(vecs_x[point_bin]))
    #     slow_ys_max.append(np.sum(vecs_y[point_bin]))
    #
    #     slow_xs_mean.append(np.mean(vecs_x[point_bin]))
    #     slow_ys_mean.append(np.mean(vecs_y[point_bin]))
    #
    #     mean_slow = np.mean(slow_diffs[point_bin])
    #     mean_slows.append(mean_slow)



    # need to sort longitudes into a new order so they are consecutive
    # other arrays also need to be sorted.... argh!
    # latitudes are already sorted
    los_sorted_idx = np.argsort(los)
    los_sorted = los[los_sorted_idx]
    los_islands_sorted = los_islands[los_sorted_idx]
    slow_xs_mean_sorted =  slow_xs_mean[los_sorted_idx]
    slow_ys_mean_sorted = slow_ys_mean[los_sorted_idx]
    # ok replace nans with 1s and numbers with 0

    consecutive_los_bins = zero_runs(los_islands_sorted)
    consecutive_las_bins = zero_runs(las_islands)

    gradients_x = np.zeros(los_sorted.shape)
    gradients_y = np.zeros(las.shape)

    curls_x = np.zeros(los_sorted.shape)
    curls_y = np.zeros(las.shape)

    gradients_x[:] = np.nan
    gradients_y[:] = np.nan
    curls_x[:] = np.nan
    curls_y[:] = np.nan

    test = []
    test_la = []
    for lo_bin in consecutive_los_bins:
        start_bin = lo_bin[0]
        end_bin = lo_bin[1]
        test.append(los_sorted[start_bin:end_bin].size)

        vecs_y_bin_c = slow_ys_mean_sorted[start_bin:end_bin]
        vecs_x_bin = slow_xs_mean_sorted[start_bin:end_bin]
        # vecs_y_bin_c = slow_ys_mean[start_bin:end_bin]
        # vecs_x_bin = slow_xs_mean[start_bin:end_bin]

        if len(vecs_x_bin) > 1:
            # lat_start = las[start_bin]
            # lat_end = las[end_bin]
            # print('lats',lat_start, lat_end)
            grad_x = np.gradient(vecs_x_bin, spacing)
            curl_x = np.gradient(vecs_y_bin_c, spacing)

            gradients_x[start_bin:end_bin] = grad_x
            curls_x[start_bin:end_bin] = curl_x

    for la_bin in consecutive_las_bins:
        start_bin_la = la_bin[0]
        end_bin_la = la_bin[1]
        test_la.append(las[start_bin_la:end_bin_la].size)
        vecs_y_bin = slow_ys_mean[start_bin_la:end_bin_la]
        vecs_x_bin_c = slow_xs_mean[start_bin_la:end_bin_la]

        if len(vecs_y_bin) > 1:
            # lon_start = los[start_bin_la]
            # lon_end = los[end_bin_la]
            # print('lons',lon_start, lon_end)

            grad_y = np.gradient(vecs_y_bin, spacing)
            curl_y = np.gradient(vecs_x_bin_c, spacing)

            gradients_y[start_bin_la:end_bin_la] = grad_y
            curls_y[start_bin_la:end_bin_la] = curl_y


    # now need to resort the gradients to be in the original
    # order... how on earth do i do this!

    # first combing the gradients and sorted los list
    los_and_gradients_x = np.vstack((los_sorted,gradients_x)).T
    los_and_curls_x = np.vstack((los_sorted,curls_x)).T

    idx1 = np.argsort(los)
    idx2_g = np.argsort(los_and_gradients_x[:, 0])
    idx2_c = np.argsort(los_and_curls_x[:, 0])
    idx1_inv = np.argsort(idx1)
    result_g = los_and_gradients_x[idx2_g][idx1_inv]
    result_c = los_and_curls_x[idx2_c][idx1_inv]

    gradients_x = result_g[:,1]
    curls_x = result_c[:,1]


    divergence = []
    curl = []

    for gx,gy in zip(gradients_x, gradients_y):

        if gx is not np.nan and gy is not np.nan:

            divergence.append(gx+gy)
        else:
            divergence.append(np.nan)

    for cx,cy in zip(curls_x, curls_y):
        if cx!=0 and cy!=0:

            curl.append(cx - cy)
        else:
            curl.append(np.nan)

    divergence = np.array(divergence)
    curl = np.array(curl)
    multi_proportions = np.array(multi_proportions)

    slow_x_mean_grid = slow_xs_mean.reshape((las_centre_grid.shape[0],los_centre_grid.shape[0]))
    slow_y_mean_grid = slow_ys_mean.reshape((las_centre_grid.shape[0],los_centre_grid.shape[0]))

    Fx_dx = np.gradient(slow_x_mean_grid, spacing, axis=1)
    Fx_dy = np.gradient(slow_x_mean_grid, spacing, axis=0)

    Fy_dx = np.gradient(slow_y_mean_grid, spacing, axis=1)
    Fy_dy = np.gradient(slow_y_mean_grid, spacing, axis=0)

    divergence = np.add(Fx_dx, Fy_dy)
    curl = np.subtract(Fy_dx, Fx_dy)
    # calculate vector field divergence
    # take divergence of lon and lat in degrees
    # I guess or km as they are spaced roughly
    # every 200 km.

    # Ok! The way to do this is complex but hear me out
    # First, I will replace all non used bins with nans
    # then i will get all bins which are adjacent to
    # another. Then for each of these mini sequences
    # I will find the gradient and if there is both
    # a dx and dy I will calculate the divergence.
    # Great!


    fig = plt.figure(figsize=(12,16), facecolor='white')
    ax = fig.add_subplot(311, projection=ccrs.PlateCarree())

    cmap = mpl.cm.get_cmap("seismic_r").copy()
    cmap.set_bad("white")

    # divergence = divergence.reshape((las_centre_grid.shape[0],los_centre_grid.shape[0]))
    # curl = curl.reshape((las_centre_grid.shape[0],los_centre_grid.shape[0]))

    max_divergence = np.absolute(divergence).max()
    max_curl = np.absolute(curl).max()

    vmin = max_divergence * -1

    # c = ax2.contourf(los_centre_grid, las_centre_grid, divergence, 20,
    #          transform=ccrs.PlateCarree(), cmap='seismic', vmin=vmin, vmax=max_divergence)


    # c = ax.pcolor(los_grid, las_grid, divergence,
    #          transform=ccrs.PlateCarree(), cmap='seismic',
    #          vmin=-1, vmax=1)

    ax.coastlines(zorder=2, resolution='50m', color='black', linewidth=1, alpha=0.5)

    q = ax.quiver(
    los_filt,
    las_filt,
    np.array(slow_xs_mean_vec),
    np.array(slow_ys_mean_vec),
    # mean_mags,
    zorder=3,
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
    Y=1.03,
    U=1,
    labelpos='E',
    label='1 s/$^{\circ}$')

    ## plot plume locations


    plt.scatter(hotspot_lons_c, hotspot_lats_c,
                transform=ccrs.PlateCarree(), color='yellow',
                marker='*', s=300, edgecolor='black', zorder=3,
                label='Plume locations (Courtillot 2003)', 
                linewidth=1.5)
    plt.legend(loc='best')
    # plt.scatter(hotspot_lons_m, hotspot_lats_m,
    #             transform=ccrs.PlateCarree(), color='Purple',
    #             marker='*', s=300, edgecolor='black', zorder=3,
    #             label='Montelli 2006')


    # cax = fig.add_axes([ax.get_position().x1 + 0.025, ax.get_position().y0,0.02, ax.get_position().height])

    # plt.colorbar(c, ax=ax, label='Divergence ($s/^{\circ^2}$)', location='right',
    #              shrink=1)



    # ax.set_xlabel('Longitude ($^{\circ}$)')
    # ax.set_ylabel('Latitude ($^{\circ}$)')



    # ax.set_extent([-131.9792592737045, -62.44319631753537, 23.44829029617145, 50.40109937197289])

    gl = ax.gridlines(draw_labels=True, linewidth=0, color='gray', alpha=0.5, linestyle='--')

    gl.xlines = True
    gl.ylines = True
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabels_bottom = True
    gl.ylabels_left = True

    gl.xlocator = mticker.FixedLocator(gird_lons)
    gl.ylocator = mticker.FixedLocator(gird_lats)

    if side == 's':
        ax.set_extent(extents_s['SPacific'])
        pass
    else:
        ax.set_extent([los_filt.min() - 15, los_filt.max() + 2, las_filt.min() - 3, las_filt.max() + 0.5])
        print([los_filt.min() - 15, los_filt.max() + 2, las_filt.min() - 3, las_filt.max() + 0.5])
    # plt.savefig(f"divergence_plot_r{region}_fl{fmin}_fh{fmax}_d{depth}_side_{side}.pdf",
    #             facecolor='white', edgecolor='black')
    # plt.tight_layout()
    #plt.show()

    # ## plot curl
    # fig = plt.figure(figsize=(14,6), facecolor='white')
    # ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

    # cmap = mpl.cm.get_cmap("seismic_r").copy()
    # cmap.set_bad("white")

    # if side == 's':
    #     ax.set_extent(extents_s['SPacific'])
    #     pass
    # else:
    #     ax.set_extent([los_filt.min() - 5, los_filt.max() + 5, las_filt.min() - 5, las_filt.max() + 0.5])

    # # divergence = divergence.reshape((las_centre_grid.shape[0],los_centre_grid.shape[0]))
    # # curl = curl.reshape((las_centre_grid.shape[0],los_centre_grid.shape[0]))

    # # c = ax2.contourf(los_centre_grid, las_centre_grid, divergence, 20,
    # #          transform=ccrs.PlateCarree(), cmap='seismic', vmin=vmin, vmax=max_divergence)

    # print(los_grid.shape, las_grid.shape, divergence.shape)


    # c = ax.pcolor(los_grid, las_grid, curl,
    #          transform=ccrs.PlateCarree(), cmap='seismic',
    #          vmin=-1, vmax=1)

    # ax.coastlines(zorder=2, resolution='50m', color='black', linewidth=1, alpha=0.5)

    # cax = fig.add_axes([ax.get_position().x1 + 0.025, ax.get_position().y0,0.02, ax.get_position().height])

    # plt.colorbar(c, cax=cax, label='Curl ($s/^{\circ^2}$)')


    # q = ax.quiver(
    # los_filt,
    # las_filt,
    # np.array(slow_xs_mean_vec),
    # np.array(slow_ys_mean_vec),
    # # mean_mags,
    # zorder=3,
    # label='Binned slowness vector residuals',
    # transform=ccrs.PlateCarree(),
    # fc = 'orange',
    # ec = 'black',
    # lw=0.75,
    # scale=4,
    # width=0.0005,
    # headwidth=20,
    # headlength = 20,
    # headaxislength=15,
    # minshaft=1,
    # minlength=0.5,
    # scale_units='inches'
    # # cmap='hot_r'
    # )

    # plt.quiverkey(q,
    # X=0.04,
    # Y=1.03,
    # U=1,
    # labelpos='E',
    # label='1 s/$^{\circ}$')

    # ax.set_xlabel('Longitude ($^{\circ}$)')
    # ax.set_ylabel('Latitude ($^{\circ}$)')

    # # ax.set_extent([-131.9792592737045, -62.44319631753537, 23.44829029617145, 50.40109937197289])

    # gl = ax.gridlines(draw_labels=True, linewidth=0, color='gray', alpha=0.5, linestyle='--')

    # gl.xlines = True
    # gl.ylines = True
    # gl.xlabels_top = False
    # gl.ylabels_right = False
    # gl.xlabels_bottom = True
    # gl.ylabels_left = True

    # gl.xlocator = mticker.FixedLocator(gird_lons)
    # gl.ylocator = mticker.FixedLocator(gird_lats)
    # plt.savefig(f"curl_plot_r{region}_fl{fmin}_fh{fmax}_d{depth}_side_{side}.pdf",
    #             facecolor='white', edgecolor='black')

    # #plt.show()

    multi_proportions = multi_proportions.reshape((las_centre_grid.shape[0],los_centre_grid.shape[0]))
    cmap = mpl.cm.get_cmap("viridis").copy()
    cmap.set_bad("white")


    # fig2 = plt.figure(figsize=(14,6), facecolor='white')
    ax2 = fig.add_subplot(313, projection=ccrs.PlateCarree(), sharex=ax)

    m = ax2.pcolor(los_grid, las_grid, multi_proportions,
                   transform=ccrs.PlateCarree(), cmap=cmap, vmin=0, vmax=1)

    # cax2 = fig.add_axes([ax2.get_position().x1 + 0.05, ax2.get_position().y0,0.02, ax2.get_position().height])

    ax2.coastlines(zorder=1, resolution='50m', color='black', linewidth=1, alpha=0.5)


    plt.colorbar(m, ax = ax2, label='Multipathing Proportion', location='right', 
                 shrink=1)

    # plt.savefig(f"./divergence_{depth}_{region}.pdf")
    # os.rename(f"./vector_bin_{depth}.pdf", f"{dir}/vector_bin_{depth}.pdf")


    gl2 = ax2.gridlines(draw_labels=True, linewidth=0.0, color='gray', alpha=0.5, linestyle='--')

    gl2.xlines = True
    gl2.ylines = True
    gl2.xlabels_top = False
    gl2.ylabels_right = False
    gl2.xlabels_bottom = True
    gl2.ylabels_left = True

    gl2.xlocator = mticker.FixedLocator(gird_lons)
    gl2.ylocator = mticker.FixedLocator(gird_lats)

    if side == 's':
        ax2.set_extent(extents_s['SPacific'])
        pass
    else:
        ax2.set_extent([los_filt.min() - 3, los_filt.max() + 2, las_filt.min() - 3, las_filt.max() + 0.5])



    ax3 = fig.add_subplot(312, projection=ccrs.PlateCarree(), sharex=ax)

    max_dvs = max(abs(np.min(dvs_depth)), np.max(dvs_depth))

    if depth <= 1200:
        m = ax3.contourf(lon, lat, dvs_depth, levels=20,
                        transform=ccrs.PlateCarree(), cmap='seismic_r', 
                        vmin=-1*max_dvs, vmax=max_dvs)
    else:
        m = ax3.contourf(lon, lat, dvs_depth, levels=30,
                        transform=ccrs.PlateCarree(), cmap='seismic_r', 
                        vmin=-1*max_dvs, vmax=max_dvs)

    # cax3 = fig.add_axes([ax3.get_position().x1 + 0.05, ax3.get_position().y0,0.02, ax3.get_position().height])

    ax3.coastlines(zorder=1, resolution='50m', color='black', linewidth=1, alpha=1)


    plt.colorbar(m, ax=ax3, label='$\delta V_{S} (\%)$', location='right', shrink=1)

    # plt.savefig(f"./divergence_{depth}_{region}.pdf")
    # os.rename(f"./vector_bin_{depth}.pdf", f"{dir}/vector_bin_{depth}.pdf")


    gl3 = ax3.gridlines(draw_labels=True, linewidth=0.0, color='black', alpha=0.5, linestyle='--')

    gl3.xlines = True
    gl3.ylines = True
    gl3.xlabels_top = False
    gl3.ylabels_right = False
    gl3.xlabels_bottom = True
    gl3.ylabels_left = True

    gl3.xlocator = mticker.FixedLocator(gird_lons)
    gl3.ylocator = mticker.FixedLocator(gird_lats)

    if side == 's':
        ax3.set_extent(extents_s['SPacific'])
        pass
    else:
        ax3.set_extent([los_filt.min() - 3, los_filt.max() + 2, las_filt.min() - 3, las_filt.max() + 0.5])



    ax.set_title("A: Slowness Vector Deviation", fontsize=16)
    ax2.set_title('C: Multipathing Proportion', fontsize=16)
    ax3.set_title(f'B: Tomography model {tomo_model}', fontsize=16)
    ax.set_extent([los_filt.min() - 5, los_filt.max() + 2, las_filt.min() - 3, las_filt.max() + 0.5])
    ax2.set_extent([los_filt.min() - 5, los_filt.max() + 2, las_filt.min() - 3, las_filt.max() + 0.5])
    ax3.set_extent([los_filt.min() - 5, los_filt.max() + 2, las_filt.min() - 3, las_filt.max() + 0.5])

    plt.tight_layout()
    # ax2.set_extent([-131.9792592737045, -62.44319631753537, 23.44829029617145, 50.40109937197289])
    # plt.savefig(f"multipathing_plot_r{region}_fl{fmin}_fh{fmax}_d{depth}_side_{side}.pdf")
    
    
    # UNCOMMENT TO SAVE SUMMARY
    # plt.savefig(f"summary_plot_r{region}_fl{fmin}_fh{fmax}_d{depth}_side_{side}.pdf",
    #             transparent=True)

    # plt.tight_layout()
    plt.show()
