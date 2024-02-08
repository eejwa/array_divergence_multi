#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
from circ_beam import calculate_locus
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--depth", help="depth to calculate pierce points at", action="store", type=float, required=True)
parser.add_argument("-fmin", "--min_frequency", help="minimum frequency", action="store", type=float, required=True)
parser.add_argument("-fmax", "--max_frequency", help="maximum frequency", action="store", type=float, required=True)
args = parser.parse_args()

depth = args.depth
fmin = args.min_frequency
fmax = args.max_frequency

original_df = pd.read_csv(f"./Plotting_file_{fmin:.2f}_{fmax:.2f}_{depth}km.txt", sep=' ', index_col=False)
multi_df = original_df[(original_df['multi'] == 'y') | (original_df['multi'] == 'm') ]
locus_file = f"Locus_Lines_{fmin}_{fmax}_{depth}.txt"

header_list = list(original_df.columns.values)
header_list.append('Phi_1')
header_list.append('Phi_2\n')
header = ' '.join(header_list)
with open(locus_file, 'w') as header_file:
    header_file.write(header)

slow_pred = original_df['slow_pred'].to_numpy()
baz_pred = original_df['baz_pred'].to_numpy()

pred_x = slow_pred * np.sin(np.radians(baz_pred))
pred_y = slow_pred * np.cos(np.radians(baz_pred))

del_x = original_df['del_x_slow'].to_numpy()
del_y = original_df['del_y_slow'].to_numpy()


found_stlas = []
found_names = []
with open(locus_file, 'a') as infile:

    for i, row in multi_df.iterrows():
        stla = row['stla_mean']
        name = row['Name']
        evt = name.split("_")[0] + "_" + name.split("_")[1]


        # if stla not in found_stlas and evt not in found_names:

            # get rows in dataframe which have the name and stla

        temp_df = multi_df[(multi_df['stla_mean'] == stla) & (multi_df['Name'].str.contains(evt)) ]
        slow_pred = temp_df.iloc[0]['slow_pred']
        baz_pred = temp_df.iloc[0]['baz_pred']

        del_x = temp_df.iloc[0]['del_x_slow']
        del_y = temp_df.iloc[0]['del_y_slow']

        pred_x = slow_pred * np.sin(np.radians(baz_pred))
        pred_y = slow_pred * np.cos(np.radians(baz_pred))

        P1_x = pred_x + del_x
        P1_y = pred_y + del_y

        P1 = [P1_x, P1_y]


        for j, row in temp_df.iterrows():

            slow_pred = row['slow_pred']
            baz_pred = row['baz_pred']

            del_x = row['del_x_slow']
            del_y = row['del_y_slow']

            pred_x = slow_pred * np.sin(np.radians(baz_pred))
            pred_y = slow_pred * np.cos(np.radians(baz_pred))

            P2_x = pred_x + del_x
            P2_y = pred_y + del_y

            if P2_x != P1_x and P2_y != P1_y:
                P2 = [P2_x, P2_y]
                Theta, Midpoint, Phi_1, Phi_2 = calculate_locus(P1, P2)
                line_row = list(row)
                line_row.append(Phi_1)
                line_row.append(Phi_2)
                line_row.append('\n')
                infile.write(' '.join(map(str, line_row)))

        found_stlas.append(stla)
        found_names.append(evt)
