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


def start_plot(dim_names):
    plt.clf()
    fig, ax = plt.subplots()  
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    divider = make_axes_locatable(ax)
    ax.set_xlabel(dim_names[0])
    ax.set_ylabel(dim_names[1])
    ax.xaxis.set_tick_params(labeltop=False)
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.yaxis.set_tick_params(labelleft=False)
    ax.grid(False)
    return fig, ax, divider

def add_annotations(ax, label_centers):
    for grp, center in label_centers.items():
        ax.annotate(grp, xy=center, xytext=center)

def main():
    gn = Granatum()

    sample_coords = gn.get_import("viz_data")
    df = gn.pandas_from_assay(gn.get_import("assay")).T
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
    label_groups = gn.get_arg("label_groups")

    coords = sample_coords.get("coords")
    dim_names = sample_coords.get("dimNames")
    if merge_genes:
        overlay_genes = False

    if labels is not None:
        label_inv = invert_dict(labels)
        label_inv = {k:list(set(v).intersection(set(df.index))) for k, v in label_inv.items()}

    # Set up colors
    cmaps = []
    if overlay_genes:         # Multiple colors required
        if max_colors == "":
            numcolors = len(gene_ids.split(','))
            cycol = cycle('bgrcmk')
            for i in range(numcolors):
                cmaps = cmaps + [LinearSegmentedColormap("fire", produce_cdict(next(cycol), grey=grey_level, min_alpha=min_alpha, max_alpha=max_alpha), N=256)]
        else:
            for col in max_colors.split(','):
                cmaps = cmaps + [LinearSegmentedColormap("fire", produce_cdict(col.strip(), grey=grey_level, min_alpha=min_alpha, max_alpha=max_alpha), N=256)]

    else:
        if max_colors == "":
            cmaps = cmaps + [LinearSegmentedColormap("fire", cdict, N=256)]
        else:
            for col in max_colors.split(','):
                cmaps = cmaps + [LinearSegmentedColormap("fire", produce_cdict(col.strip(), grey=grey_level, min_alpha=min_alpha, max_alpha=max_alpha), N=256)]

    gene_index = -1
    numgenes = len(gene_ids.split(','))    # Probably can drop this

    # Split gene list first into gene names with corresponding scale
    gene_scale_tuples = [(gene_eqs.strip(), 1.0) if len(gene_eqs.strip().split("="))<2 else tuple(gene_eqs.strip().split("=")) for gene_eqs in gene_ids.split(',')]
    gene_list = [gene for gene, scale in gene_scale_tuples]
    cmap_genes = {gene_list[i]:cmaps[i] for i in range(len(cmaps))}

    # if the gene ID entered is not present in the assay
    # Communicate it to the user and output a table of available gene ID's
    if len(set(gene_list).difference(set(df.columns))) > 0:
        description = "The following gene(s) {} is/are not present in the assay. See the step that generated the assay".format(set(gene_list).difference(set(df.columns)))
        genes_in_assay = pd.DataFrame(df.columns.tolist(), columns=['Gene unavailable in assay: choose from below'])
        gn.add_pandas_df(genes_in_assay, description)
        gn.commit()
        return
    
    # Create scatter_df, extra column if merging, and show counts in a column if merging
    merge_scatters = [df[gene]*scale for gene, scale in gene_scale_tuples]  # will be the list of scatter_dfs to merge
    values = pd.concat(merge_scatters, axis=1)
    if merge_genes:
        values[gene_ids] = values.sum(axis=1)
    if convert_to_zscore:
        values = (values - values.mean()) / values.std(ddof=0)   # Will do all columns at once
    if labels is not None:
        percentages = pd.concat([((values.loc[v, :] >= threshold).sum()+0.0).to_frame().T*100.0/len(v) for k, v in label_inv.items()], axis=0).round(2)
        percentages["# Cells"] = [len(v) for k, v in label_inv.items()]
        percentages.index = list(label_inv.keys())
    else:
        percentages = ((values >= threshold).sum()+0.0).to_frame().T*100.0 / values.shape[0]
        percentages["# Cells"] = values.shape[0]
        percentages.index = ["All (no labels provided)"]
    gn.add_result("* Percentages of each group expressing at threshold {} *".format(threshold), "markdown")
    gn.add_pandas_df(percentages.reset_index())

    # Export files of the data for other processing
    gn.export(percentages.to_csv(), "percentage_summaries.csv", kind='raw', meta=None, raw=True)
    if labels is not None:
        for k, v in label_inv.items():
            gn.export(values.loc[v, :].to_csv(), "cell_values_{}.csv".format(k), kind='raw', meta=None, raw=True)
    else:
        gn.export(values.to_csv(), "cell_values.csv", kind='raw', meta=None, raw=True)

    x_df = pd.DataFrame({"x": [a[0] for a in coords.values()]}, index=coords.keys())
    y_df = pd.DataFrame({"y": [a[1] for a in coords.values()]}, index=coords.keys())

    # Get the centers of each cluster on the coordinates
    if labels is not None:
        label_centers = {k:(x_df.loc[v, "x"].mean(), y_df.loc[v, "y"].mean()) for k, v in label_inv.items()}

    values_df = np.clip(values, min_level, max_level, out=None)
    values_df = values_df.loc[x_df.index, :]       # Ensure alignment
    min_value = np.nanmin(values_df)
    max_value = np.nanmax(values_df)
    scaled_marker_size = (max_marker_area-min_marker_area)*(values_df-min_value)/(max_value-min_value) + min_marker_area
    scaled_marker_size = scaled_marker_size*scaled_marker_size
    print(values_df.head(), flush=True)
    print(x_df.head(), flush=True)
    print(y_df.head(), flush=True)

    # Now we can easily generate the scatter plots
    if merge_genes:
        # Generate scatter with just merged genes
        fig, ax, divider = start_plot(dim_names)
        scatter = ax.scatter(x=x_df["x"], y=y_df["y"], s=scaled_marker_size[gene_ids], c=values_df[gene_ids], cmap=cmaps[len(cmaps)-1]) #Amp_3.mpl_colormap)
        cax = divider.append_axes('bottom', size=0.15, pad=0.01)
        cax.tick_params(axis="x",direction="inout", pad=-1)
        cbar = fig.colorbar(scatter, cax=cax, orientation='horizontal', aspect=300)
        if labels is not None and label_groups:
            add_annotations(ax, label_centers)
        gn.add_current_figure_to_results("Scatter-plot of summed {} expression".format(gene_ids), dpi=75)
    elif overlay_genes:
        fig, ax, divider = start_plot(dim_names)
        cax = None # Will be the last color axis added
        for gene in gene_list:
            scatter = ax.scatter(x=x_df["x"], y=y_df["y"], s=scaled_marker_size[gene], c=values_df[gene], cmap=cmap_genes[gene])
            cax = divider.append_axes('bottom', size=0.15, pad=0.01)
            cbar = fig.colorbar(scatter, cax=cax, orientation='horizontal', aspect=300)
            cbar.ax.set_ylabel(gene, rotation=0)
            cax.yaxis.set_label_coords(0.08, 0.0)
            cax.tick_params(axis="x",direction="in", pad=-1)
            cax.xaxis.set_ticklabels([])
        cax.tick_params(axis="x",direction="inout", pad=-1)
        if labels is not None and label_groups:
            add_annotations(ax, label_centers)
        gn.add_current_figure_to_results("Scatter-plot of {} expression".format(gene_ids), dpi=75)
    else:
        for gene in gene_list:
            fig, ax, divider = start_plot(dim_names)
            scatter = ax.scatter(x=x_df["x"], y=y_df["y"], s=scaled_marker_size[gene], c=values_df[gene], cmap=cmap_genes[gene])
            cax = divider.append_axes('bottom', size=0.15, pad=0.01)
            cbar = fig.colorbar(scatter, cax=cax, orientation='horizontal', aspect=300)
            cbar.ax.set_ylabel(gene, rotation=0)
            cax.yaxis.set_label_coords(0.08, 0.0)
            cax.tick_params(axis="x",direction="inout", pad=-1)
            cax.xaxis.set_ticklabels([])
            if labels is not None and label_groups:
                add_annotations(ax, label_centers)
            gn.add_current_figure_to_results("Scatter-plot of {} expression".format(gene), dpi=75)

    gn.commit()

if __name__ == "__main__":
    # Try except block to send an email about error #
    try:
        main()
    except:
        error_message = traceback.format_exc()
        sys.stderr.write(error_message) # Write the error to stderr anyway so the user can see what went wrong
        bug_report("Color Scatter-Plot", "amantrav@umich.edu", error_message)
