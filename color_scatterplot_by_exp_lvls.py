#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

def invert_dict(my_map):
    inv_map = {}
    for k, v in my_map.items():
        inv_map[v] = inv_map.get(v, []) + [k]
    return inv_map

def show_percentages(gn, invdict, values, threshold):
    percent_df = pd.concat([((values.loc[v, :] >= threshold).sum()+0.0)/len(v) for k, v in invdict.items()], axis=0)
    gn.add_result(" * Shape = {} * ".format(percent_df.shape), "markdown")
    percent_df.index = list(invdict.keys())
    gn.add_pandas_df(percent_df.reset_index())

def main():
    gn = Granatum()

    sample_coords = gn.get_import("viz_data")
    df = gn.pandas_from_assay(gn.get_import("assay"))
    labels = gn.get_import("labels")
    gene_ids = gn.get_arg("gene_ids")
    overlay_genes = gn.get_arg("overlay_genes")
    merge_genes = gn.get_arg("merge_genes")
    max_colors = gn.get_arg("max_colors")
    min_level = gn.get_arg("min_level")
    max_level = gn.get_arg("max_level")
    convert_to_zscore = gn.get_arg("convert_to_zscore")
    min_marker_area = gn.get_arg("min_marker_area")
    max_marker_area = gn.get_arg("max_marker_area")
    min_alpha = gn.get_arg("min_alpha")
    max_alpha = gn.get_arg("max_alpha")
    grey_level = gn.get_arg("grey_level")
    threshold = gn.get_arg("threshold")

    coords = sample_coords.get("coords")
    dim_names = sample_coords.get("dimNames")
    if merge_genes:
        overlay_genes = False

    if labels is not None:
        label_inv = invert_dict(labels)
        label_inv = {k:list(set(v).intersection(set(df.index))) for k, v in label_inv.items()}

    cmaps = []
    if overlay_genes and not merge_genes:         # Multiple colors required
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
        if max_colors == "":
            cmaps = cmaps + [LinearSegmentedColormap("fire", cdict, N=256)]
        else:
            for col in max_colors.split(','):
                col = col.strip()
                cmaps = cmaps + [LinearSegmentedColormap("fire", produce_cdict(col, grey=grey_level, min_alpha=min_alpha, max_alpha=max_alpha), N=256)]

    colorbar_height = 10
    plot_height = 650
    num_cbars = 1
    if overlay_genes:
        num_cbars = len(gene_ids.split(','))
    cbar_height_ratio = plot_height/(num_cbars*colorbar_height)
    fig, ax = plt.subplots() # plt.subplots(1+num_cbars, 1, gridspec_kw={'height_ratios': [cbar_height_ratio] + [1]*num_cbars})
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    divider = make_axes_locatable(ax)

    gene_index = -1
    numgenes = len(gene_ids.split(','))
    transposed_df = df.T

    merge_scatters = []  # will be the list of scatter_dfs to merge

    # Only do merge_genes plot at the end
    for gene_id in gene_ids.split(','):
        gene_id = gene_id.strip()
        gene_ids_equal = gene_id.split("=")
        if len(gene_ids_equal) > 1:
            scale = float(gene_ids_equal[1])
        else:
            scale = 1.0/numgenes 
        gene_id = gene_ids_equal[0]
        gene_index = gene_index + 1
        if gene_id in df.index:
            # Data preparation
            if merge_genes and gene_index < numgenes - 1:
                merge_scatters.append(df.loc[gene_id, :])
                continue
            else:
                if merge_genes:
                    values = pd.concat(merge_scatters, axis=1).sum(axis=1) * scale
                else:
                    values = transposed_df[gene_id] * scale
                values = (values - values.mean()) / values.std(ddof=0) if convert_to_zscore else values
                scatter_df = pd.DataFrame({"x": [a[0] for a in coords.values()], "y": [a[1] for a in coords.values()], "value": values}, index=coords.keys())

            if not overlay_genes:     # New gene appears in separate plot
                plt.clf()
                fig, ax = plt.subplots()  # plt.subplots(1+num_cbars, 1, gridspec_kw={'height_ratios': [cbar_height_ratio] + [1]*num_cbars})
                ax.xaxis.set_label_position("top")
                ax.xaxis.tick_top()
                divider = make_axes_locatable(ax)

            values_df = np.clip(scatter_df["value"], min_level, max_level, out=None)
            min_value = np.nanmin(values_df)
            max_value = np.nanmax(values_df)
            scaled_marker_size = (max_marker_area-min_marker_area)*(values_df-min_value)/(max_value-min_value) + min_marker_area
            scaled_marker_size = scaled_marker_size*scaled_marker_size
            # s = 5000 / scatter_df.shape[0]
            if not (merge_genes and gene_index < numgenes - 1):
                scatter = ax.scatter(x=scatter_df["x"], y=scatter_df["y"], s=scaled_marker_size, c=values_df, cmap=cmaps[gene_index % len(cmaps)]) #Amp_3.mpl_colormap)

                if labels is not None and overlay_genes:
                    msg = "* Percentages of each group expressing {} at threshold {} *".format(gene_id, threshold)
                    gn.add_result(msg, "markdown")
                    show_percentages(gn, label_inv, values_df.to_frame(), threshold)
                cax = divider.append_axes('bottom', size=0.15, pad=0.01)
                cbar = fig.colorbar(scatter, cax=cax, orientation='horizontal', aspect=300)
                #cbar = fig.colorbar(scatter, cax=ax[1+(gene_index%num_cbars)], orientation='horizontal', aspect=40)
                #cbar.set_label(gene_id, rotation=0)

                ax.set_xlabel(dim_names[0])
                ax.set_ylabel(dim_names[1])

                if not merge_genes:
                    cbar.ax.set_ylabel(gene_id, rotation=0)
                    cax.yaxis.set_label_coords(0.08, 0.0)

                ax.xaxis.set_tick_params(labeltop=False)
                ax.xaxis.set_tick_params(labelbottom=False)
                ax.yaxis.set_tick_params(labelleft=False)
                ax.grid(False)

            if merge_genes:
                if gene_index == numgenes - 1:
                    cax.tick_params(axis="x",direction="inout", pad=-1)
                    gn.add_current_figure_to_results("Scatter-plot of {} expression".format(gene_ids), dpi=75)
                    if labels is not None:
                        msg = "* Percentages of each group expressing merge at threshold {} *".format(threshold)
                        gn.add_result(msg, "markdown")
                        show_percentages(gn, label_inv, values_df.to_frame(), threshold)
            elif not overlay_genes:
                cax.tick_params(axis="x",direction="inout", pad=-1)
                gn.add_current_figure_to_results("Scatter-plot of {} expression".format(gene_id), dpi=75)
                if labels is not None:
                    msg = "* Percentages of each group expressing {} at threshold {} *".format(gene_id, threshold)
                    gn.add_result(msg, "markdown")
                    show_percentages(gn, label_inv, values_df.to_frame(), threshold)
            else:
                if gene_index < numgenes-1:
                    cax.tick_params(axis="x",direction="in", pad=-1)
                    cax.xaxis.set_ticklabels([])
                else:
                    cax.tick_params(axis="x",direction="inout", pad=-1)
        else:
            # if the gene ID entered is not present in the assay
            # Communicate it to the user and output a table of available gene ID's
            description = 'The selected gene is not present in the assay. See the step that generated the assay'
            genes_in_assay = pd.DataFrame(df.index.tolist(), columns=['Gene unavailable in assay: choose from below'])
            gn.add_pandas_df(genes_in_assay, description)
    if overlay_genes:
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
