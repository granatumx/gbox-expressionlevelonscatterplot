#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from granatum_sdk import Granatum
from palettable.cmocean.sequential import Amp_3
from matplotlib.colors import LinearSegmentedColormap

import os
import traceback
import sys

from gbox_py_helpers import bug_report

COLORS = ["#3891ea", "#29ad19", "#ac2d58", "#db7580", "#ed2310", "#ca2dc2", "#5f7575", "#7cc1b5", "#c3bd78", "#4ffa24"]

cdict = {'red':   [(0.0,   0.86, 0.86), (0.192, 0.95, 0.95), (0.385, 0.95, 0.95), (0.577, 0.88, 0.88), (1.0,   0.7,  0.7)], 'green': [(0.0,   0.86, 0.86), (0.192, 0.77, 0.77), (0.385, 0.64, 0.64), (0.577, 0.45, 0.45), (1.0,   0.04,  0.04)], 'blue':  [(0.0,   0.86, 0.86), (0.192, 0.63, 0.63), (0.385, 0.43, 0.43), (0.577, 0.32, 0.32), (1.0,   0.1,  0.1)]}

def main():
    gn = Granatum()

    sample_coords = gn.get_import("viz_data")
    df = gn.pandas_from_assay(gn.get_import("assay"))
    gene_ids = gn.get_arg("gene_ids")
    min_level = gn.get_arg("min_level")
    max_level = gn.get_arg("max_level")
    convert_to_zscore = gn.get_arg("convert_to_zscore")

    coords = sample_coords.get("coords")
    dim_names = sample_coords.get("dimNames")

    for gene_id in gene_ids.split(',')
        gene_id = gene_id.strip()
        if gene_id in df.index:

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

            plt.scatter(x=scatter_df["x"], y=scatter_df["y"], s=5000 / scatter_df.shape[0], c=np.clip(scatter_df["value"], min_level, max_level, out=None), cmap=LinearSegmentedColormap("fire", cdict, N=256)) #Amp_3.mpl_colormap)
            plt.colorbar()

            plt.xlabel(dim_names[0])
            plt.ylabel(dim_names[1])
            plt.tight_layout()

            gn.add_current_figure_to_results("Scatter-plot", dpi=75)

            gn.commit()

        else:

        # if the gene ID entered is not present in the assay
        # Communicate it to the user and output a table of available gene ID's
        
            description = 'The selected gene is not present in the assay. See the step that generated the assay'
            genes_in_assay = pd.DataFrame(df.index.tolist(), columns=['Gene unavailable in assay: choose from below'])
            gn.add_pandas_df(genes_in_assay, description)

            gn.commit()

if __name__ == "__main__":
    # Try except block to send an email about error #
    try:
        main()
    except:
        error_message = traceback.format_exc()
        sys.stderr.write(error_message) # Write the error to stderr anyway so the user can see what went wrong
        bug_report("Color Scatter-Plot", "amantrav@umich.edu", error_message)
