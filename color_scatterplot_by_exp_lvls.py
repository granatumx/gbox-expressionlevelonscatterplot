#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from granatum_sdk import Granatum

import os
import traceback
import sys

from gbox_py_helpers import bug_report

COLORS = ["#3891ea", "#29ad19", "#ac2d58", "#db7580", "#ed2310", "#ca2dc2", "#5f7575", "#7cc1b5", "#c3bd78", "#4ffa24"]


def main():
    gn = Granatum()

    sample_coords = gn.get_import("viz_data")
    df = gn.pandas_from_assay(gn.get_import("assay"))
    gene_id = gn.get_arg("gene_id")

    coords = sample_coords.get("coords")
    dim_names = sample_coords.get("dimNames")

    if gene_id in df.index:

        scatter_df = pd.DataFrame(
            {"x": [a[0] for a in coords.values()], "y": [a[1] for a in coords.values()], "value": df.loc[gene_id, :]},
            index=coords.keys(),
        )

        plt.scatter(x=scatter_df["x"], y=scatter_df["y"], s=5000 / scatter_df.shape[0], c=scatter_df["value"], cmap="Reds")
        plt.colorbar()

        plt.xlabel(dim_names[0])
        plt.ylabel(dim_names[1])
        plt.tight_layout()

        gn.add_current_figure_to_results("Scatter-plot", dpi=75)

        gn.commit()

    else:

        # if the gene ID entered is not present in the assay
        # Communicate it to the user and output a table of available gene ID's

        """genes_in_assay = np.transpose(df.index.values)
        genes_in_assay = np.reshape(genes_in_assay, (genes_in_assay.shape[0], 1)) # reshape to make it a 2D array
        print(genes_in_assay.shape)

        fig, axs = plt.subplots(2)

        message = 'The selected gene is not present in the assay\nSee the step that generated the assay'
        axs[0].text(0.1, 0.5, message, fontsize=16)
        axs[0].axis('off')

        # create table
        axs[1].table(cellText=genes_in_assay, cellLoc='center', colLabels=['Gene'])
        axs[1].axis('off')

        gn.add_current_figure_to_results('error message', dpi=75)"""

        
        description = 'The selected gene is not present in the assay\nSee the step that generated the assay'
        genes_in_assay = pd.DataFrame(df.index.tolist(), columns=['Gene'])
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
