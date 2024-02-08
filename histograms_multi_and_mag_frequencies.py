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
import matplotlib.colors as colors
import matplotlib.patches as mpl_patches

extents = {"US":[-170,-60,10,75], "EU":[-35,60,10,70], "SAm":[-120,-60,-70,-5], "SAs":[60,140,0,60], "Global":[-175,175,-85,85]}

region='EU'

dir_multi = f"histograms_multi_{region}"
dir_mags = f"histograms_mags_{region}"
dir_bazs = f"histograms_bazs_{region}"
dir_slows = f"histograms_slows_{region}"


try:
    os.mkdir(dir)
except:
    pass

depths = np.arange(0, 2900, 100)
depths = np.append(depths, 2891)

cmap = plt.cm.get_cmap('summer')
gird_lons = np.arange(-360,375, 15)
gird_lats = np.arange(-90,105, 15)

freq_bands = ['0.10_0.20', '0.15_0.30', '0.20_0.40']

freq_band_labels = ['0.10 $-$ 0.20', '0.15 $-$ 0.30', '0.20 $-$ 0.40']
multi_counts = []
multi_props = []
mags_freq = []
bazs_freq = []
slows_freq = []

for j,f in enumerate(freq_bands):

    fmin = float(f.split('_')[0])
    fmax = float(f.split('_')[1])

    Plotting_file = f"./Plotting_file_{fmin:.2f}_{fmax:.2f}_2000.0km.txt"

    df = pd.read_csv(Plotting_file, sep=' ', index_col=False)
    df = df[(df['stlo_mean'] > extents[region][0]) & (df['stlo_mean'] < extents[region][1]) & (df['stla_mean'] > extents[region][2]) & (df['stla_mean'] < extents[region][3])]
    df = df.sort_values(by='multi')
    df = df.dropna()

    names = df['Name'].to_numpy()
    multis = df['multi'].to_numpy()
    mags = df['mag'].to_numpy()

    del_bazs = df['baz_diff'].to_numpy()
    del_bazs = (del_bazs + 180) % 360 - 180

    del_slows = df['slow_diff'].to_numpy()

    name_numbers = []
    mags_freq.append(mags)
    bazs_freq.append(del_bazs)
    slows_freq.append(del_slows)

    for name in names:
        name_numbers.append(name.split('_')[2])

    name_numbers = np.array(name_numbers)
    total_obs_multi = np.count_nonzero(name_numbers.astype(float) == 1)
    multi_count = np.count_nonzero(name_numbers.astype(float) == 2)
    multi_prop = multi_count / total_obs_multi

    multi_counts.append(multi_count)
    multi_props.append(multi_prop)


plt.style.use('ggplot')
colors = ['red', 'blue', 'green']
X_axis = np.arange(len(freq_bands))

## Multipathing frequency dependence
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
ax.bar(X_axis, multi_props)
ax.set_xticks(X_axis)
ax.set_xticklabels(freq_band_labels)

ax.set_xlabel("Frequency band (Hz)")
ax.set_ylabel("Multipathing proportion")
ax.set_title(f"Proportion of multipathing observations in each frequency band")
plt.savefig(f"./bar_multi_{region}.pdf")
plt.show()

# Slowness vector magnitudes frequency dependence
# bins_mag = np.linspace(0,2.5,26)
bins_mag = np.linspace(0,2.4,13)



fig2 = plt.figure(figsize=(12,8))
ax = fig2.add_subplot(111)
ax.hist(mags_freq, bins_mag, histtype='bar', color=colors, label=freq_band_labels)
ax.set_xlabel("Slowness vector magnitude (s/$^{\circ}$)")
ax.set_ylabel("Number of values in bin")
ax.set_title(f"binned slowness vector magnitudes in different frequency bands")
ax.legend(loc='best')
plt.savefig(f"./histogram_mag_{region}.pdf")
plt.show()

# Backazimuth with frequency
fig3 = plt.figure(figsize=(12,8))

bins_baz = np.linspace(-34,34,35)

# This is  the colormap I'd like to use.
cm = plt.cm.get_cmap('seismic_r')

col = np.linspace(0,1,35)


ax1 = fig3.add_subplot(311)

# Plot histogram.
n, bins, patches = ax1.hist(bazs_freq[0], color='green', ec='black', align='mid', bins=bins_baz)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
for c, p in zip(col, patches):
    plt.setp(p, 'facecolor', cm(c))

ax2 = fig3.add_subplot(312)

# Plot histogram.
n, bins, patches = ax2.hist(bazs_freq[1], color='green', ec='black', bins=bins_baz)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
for c, p in zip(col, patches):
    plt.setp(p, 'facecolor', cm(c))

ax3 = fig3.add_subplot(313)

# Plot histogram.
n, bins, patches = ax3.hist(bazs_freq[2], color='green', ec='black', bins=bins_baz)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
for c, p in zip(col, patches):
    plt.setp(p, 'facecolor', cm(c))

# # scale values to interval [0,1]
# col = bin_centers - min(bin_centers)
# print(col)
# col /= (max(bazs_freq[0]) - min(bin_centers))
# print(col)
handle = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white",
                                 lw=0, alpha=0)]


ax1.legend(handle,["0.10 $-$ 0.20 Hz"], loc='best', fancybox=True,
           framealpha=1, handlelength=0, handletextpad=0, facecolor='white')
ax2.legend(handle,["0.15 $-$ 0.30 Hz"], loc='best', fancybox=True,
           framealpha=1, handlelength=0, handletextpad=0, facecolor='white')
ax3.legend(handle,["0.20 $-$ 0.40 Hz"], loc='best', fancybox=True,
           framealpha=1, handlelength=0, handletextpad=0, facecolor='white')

ax1.sharex(ax3)
ax2.sharex(ax3)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)

plt.subplots_adjust(wspace=0, hspace=0)

# n, bins, patches = plt.hist(bazs_freq[0], 50, color='green')

fig3.suptitle(f"Backazimuth Deviations", fontsize=22)
fig3.supxlabel("$\delta\\theta$ ($^{\circ}$)", fontsize=18)
fig3.supylabel("Number of values in bin", fontsize=18)


plt.savefig(f"./histogram_bazs_{region}.pdf")
plt.show()

# Horizontal slowness with frequency

fig4 = plt.figure(figsize=(12,8))
cm = plt.cm.get_cmap('PRGn')

col = np.linspace(0,1,21)
bins_slow = np.linspace(-2,2,21)

ax1 = fig4.add_subplot(311)

# Plot histogram.
n, bins, patches = ax1.hist(slows_freq[0], color='green', ec='black', bins=bins_slow)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
for c, p in zip(col, patches):
    plt.setp(p, 'facecolor', cm(c))


ax2 = fig4.add_subplot(312)

# Plot histogram.
n, bins, patches = ax2.hist(slows_freq[1], color='green', ec='black', bins=bins_slow)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
for c, p in zip(col, patches):
    plt.setp(p, 'facecolor', cm(c))

ax3 = fig4.add_subplot(313)

# Plot histogram.
n, bins, patches = ax3.hist(slows_freq[2], color='green', ec='black', bins=bins_slow)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
for c, p in zip(col, patches):
    plt.setp(p, 'facecolor', cm(c))

ax1.legend(handle,["0.10 $-$ 0.20 Hz"], loc='best', fancybox=True,
           framealpha=1, handlelength=0, handletextpad=0, facecolor='white')
ax2.legend(handle,["0.15 $-$ 0.30 Hz"], loc='best', fancybox=True,
           framealpha=1, handlelength=0, handletextpad=0, facecolor='white')
ax3.legend(handle,["0.20 $-$ 0.40 Hz"], loc='best', fancybox=True,
           framealpha=1, handlelength=0, handletextpad=0, facecolor='white')

ax3.xaxis.set_ticks(np.arange(-2,2,0.4))

ax1.sharex(ax3)
ax2.sharex(ax3)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)

plt.subplots_adjust(wspace=0, hspace=0)

fig4.suptitle(f"Horizontal Slowness Deviations", fontsize=22)
fig4.supxlabel("$\delta$$p$ (s/$^{\circ}$)", fontsize=18)
fig4.supylabel("Number of values in bin", fontsize=18)

plt.savefig(f"./histogram_slow_{region}.pdf")

plt.show()
