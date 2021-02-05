#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from granatum_sdk import Granatum
from palettable.cmocean.sequential import Amp_3
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colors
from itertools import cycle

import os
import traceback
import sys

from gbox_py_helpers import bug_report

COLORS = ["#3891ea", "#29ad19", "#ac2d58", "#db7580", "#ed2310", "#ca2dc2", "#5f7575", "#7cc1b5", "#c3bd78", "#4ffa24"]

cdict = {'red':   [(0.0,   0.86, 0.86), (0.192, 0.95, 0.95), (0.385, 0.95, 0.95), (0.577, 0.88, 0.88), (1.0,   0.7,  0.7)], 'green': [(0.0,   0.86, 0.86), (0.192, 0.77, 0.77), (0.385, 0.64, 0.64), (0.577, 0.45, 0.45), (1.0,   0.04,  0.04)], 'blue':  [(0.0,   0.86, 0.86), (0.192, 0.63, 0.63), (0.385, 0.43, 0.43), (0.577, 0.32, 0.32), (1.0,   0.1,  0.1)]}

def produce_cdict(color_name, grey=0.9, min_alpha=0.2, max_alpha=1.0):
    clist = colors.to_rgba(color_name)
    g = grey
    return {
            'red': [(0.0, g, g), (0.2, g, g), (1.0, clist[0], clist[0])],
            'green': [(0.0, g, g), (0.2, g, g), (1.0, clist[1], clist[1])],
            'blue': [(0.0, g, g), (0.2, g, g), (1.0, clist[2], clist[2])],
            'alpha': [(0.0, min_alpha, min_alpha), (0.2, min_alpha, min_alpha), (1.0, max_alpha, max_alpha)]
            }

def main():
    gn = Granatum()

    sample_coords = gn.get_import("viz_data")
    df = gn.pandas_from_assay(gn.get_import("assay"))
    gene_ids = gn.get_arg("gene_ids")
    overlay_genes = gn.get_arg("overlay_genes")
    max_colors = gn.get_arg("max_colors")
    min_level = gn.get_arg("min_level")
    max_level = gn.get_arg("max_level")
    convert_to_zscore = gn.get_arg("convert_to_zscore")
    min_marker_area = gn.get_arg("min_marker_area")
    max_marker_area = gn.get_arg("max_marker_area")
    min_alpha = gn.get_arg("min_alpha")
    max_alpha = gn.get_arg("max_alpha")
    grey_level = gn.get_arg("grey_level")

    coords = sample_coords.get("coords")
    dim_names = sample_coords.get("dimNames")

    cmaps = []
    if overlay_genes:
        if max_colors == "":
            numcolors = len(gene_ids.split(','))
            cycol = cycle('bgrcmk')
            for i in range(numcolors):
                cmaps = cmaps + [LinearSegmentedColormap("fire", produce_cdict(next(cycol), grey=grey_level, min_alpha=min_alpha, max_alpha=max_alpha), N=256)]
        else:
            for col in max_colors.split(','):
                col = col.strip()
                cmaps = cmaps + [LinearSegmentedColormap("fire", produce_cdict(col, grey=grey_level, min_alpha=min_alpha, max_alpha=max_alpha), N=256)]

    else:
        cmaps = cmaps + [LinearSegmentedColormap("fire", cdict, N=256)]

    plt.clf()
    gene_index = -1
    for gene_id in gene_ids.split(','):
        gene_id = gene_id.strip()
        gene_index = gene_index + 1
        if gene_id in df.index:
            if not overlay_genes:
                plt.clf()

            transposed_df = df.T

            mean = transposed_df[gene_id].mean()
            stdev = transposed_df[gene_id].std(ddof=0)

            if convert_to_zscore:
                scatter_df = pd.DataFrame(
                    {"x": [a[0] for a in coords.values()], "y": [a[1] for a in coords.values()], "value": (df.loc[gene_id, :]-mean)/stdev},
                    index=coords.keys())
            else:
                scatter_df = pd.DataFrame(
                    {"x": [a[0] for a in coords.values()], "y": [a[1] for a in coords.values()], "value": df.loc[gene_id, :]},
                    index=coords.keys())

            values_df = np.clip(scatter_df["value"], min_level, max_level, out=None)
            min_value = np.nanmin(values_df)
            max_value = np.nanmax(values_df)
            scaled_marker_size = (max_marker_area-min_marker_area)*(scatter_df["value"]-min_value)/(max_value-min_value) + min_marker_area
            scaled_marker_size = scaled_marker_size*scaled_marker_size
            # s = 5000 / scatter_df.shape[0]
            plt.scatter(x=scatter_df["x"], y=scatter_df["y"], s=scaled_marker_size, c=values_df, cmap=cmaps[gene_index % len(cmaps)]) #Amp_3.mpl_colormap)
            cbar = plt.colorbar(orientation='horizontal', pad=0.1, aspect=40)
            cbar.set_label(gene_id, rotation=0)

            plt.xlabel(dim_names[0])
            plt.ylabel(dim_names[1])

            if not overlay_genes:
                plt.tight_layout()

                gn.add_current_figure_to_results("Scatter-plot of {} expression".format(gene_id), dpi=75)


        else:

        # if the gene ID entered is not present in the assay
        # Communicate it to the user and output a table of available gene ID's
        
            description = 'The selected gene is not present in the assay. See the step that generated the assay'
            genes_in_assay = pd.DataFrame(df.index.tolist(), columns=['Gene unavailable in assay: choose from below'])
            gn.add_pandas_df(genes_in_assay, description)
    if overlay_genes:
        plt.tight_layout()
        gn.add_current_figure_to_results("Scatter-plot of {} expression".format(gene_ids), height=650+100*len(gene_ids.split(',')), dpi=75)

    gn.commit()

if __name__ == "__main__":
    # Try except block to send an email about error #
    try:
        main()
    except:
        error_message = traceback.format_exc()
        sys.stderr.write(error_message) # Write the error to stderr anyway so the user can see what went wrong
        bug_report("Color Scatter-Plot", "amantrav@umich.edu", error_message)
