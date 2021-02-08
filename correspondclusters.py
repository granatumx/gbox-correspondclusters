#!/usr/bin/env python

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from colour import Color
from matplotlib.patches import Polygon
import statistics as st
import time
from random import randrange
import re
import sys

from granatum_sdk import Granatum

def main():
    tic = time.perf_counter()

    gn = Granatum()
    assay = gn.pandas_from_assay(gn.get_import('assay'))
    groups = gn.get_import('groups')
    reflabels = gn.get_import('reflabels')
    remove_cells = gn.get_arg('remove_cells')

    inv_map = {}
    for k, v in groups.items():
        inv_map[v] = inv_map.get(v, []) + [k]

    inv_map_ref = {}
    for k, v in reflabels.items():
        inv_map_ref[v] = inv_map_ref.get(v, []) + [k]
    
    group_relabel = {}
    mislabelled_cells = []
    for k, v in inv_map.items():
        vset = set(v)
        label_scores = {}
        for kref, vref in inv_map_ref.items():
            label_scores[kref] = len(set(vref).intersection(vset))
        group_relabel[k] = max(label_scores, key=label_scores.get)
        mislabelled_cells.append(list(vset.difference(set(inv_map_ref[group_relabel[k]]))))

    if remove_cells:
        gn.add_result("Dropping {} mislabelled cells".format(len(mislabelled)), "markdown")
        assay = assay.drop(mislabelled_cells, axis=1)
        groups = {key:val for key, val in groups.items() if val != ds}

    for cell in groups:
        groups[cell] = group_relabel[groups[cell]]
    
    toc = time.perf_counter()
    time_passed = round(toc - tic, 2)

    gn.export_statically(gn.assay_from_pandas(assay), "Corresponded assay")
    gn.export_statically(groups, "Corresponded labels")

    timing = "* Finished sample coloring step in {} seconds*".format(time_passed)
    gn.add_result(timing, "markdown")

    gn.commit()


if __name__ == "__main__":
    main()
