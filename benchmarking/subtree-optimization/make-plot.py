"""Build a comparison plot showing the average and range of tree accumulations during sample-optimize-merge pipeline.

Put replicate csv logfiles in directories named in `directories` list, and name the field to be plotted agains `Iterations` in `plot_field`."""

import pandas
import csv
import statistics
from pathlib import Path
import matplotlib.pyplot as plt
import sys
from matplotlib.patches import Patch
sys.set_int_max_str_digits(10000)

prefix='20D_intermediate/'
directories = [
    "option_--sample-best-tree",
    "option_--sample-best-tree-s10",
    "option_--sample-best-tree-s0",
    "option_--sample-best-tree-s2",
]

logtransform = (lambda n: statistics.log(n, 10), "Log10")
notransform = (lambda n: n, "")

plot_fields = [("MaxParsimony", notransform),
               ("WorstParsimony", notransform),
               ("NTreesMaxParsimony", logtransform),
               ("NTrees", logtransform),
               ("NEdges", notransform),
               ("NNodes", notransform),
               ("Iteration", notransform),
               ]


colors = [
    '#1b9e77',
    '#d95f02',
    '#7570b3',
    '#e7298a',
    '#66a61e',
    '#e6ab02',
]


fields = ['Iteration', 'NTrees', 'NNodes', 'NEdges', 'MaxParsimony', 'NTreesMaxParsimony', 'WorstParsimony', 'SecondsElapsed']
iterations = 501

def floor_mean(numbers):
    s = sum(numbers)
    l = len(numbers)
    try:
        return s / l
    except OverflowError:
        return s // l

mean = floor_mean

fig, axarr = plt.subplots(len(plot_fields), sharex=False, figsize=(7, 4 * len(plot_fields)))

handles = []
for color, directory in zip(colors, directories):
    handles.append(Patch(color=color, label=directory))
    p = Path(prefix + directory)
    for logfile in p.glob('log*/logfile.csv'):
        data = [
            [] for _ in range(len(fields))
        ]
        with open(logfile, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            try:
                next(reader)
            except StopIteration:
                continue
            for iteration_idx, row in enumerate(reader):
                for field_idx, dat in enumerate(row):
                    data[field_idx].append(0 if dat == '' else int(dat))
        data_dict = dict(zip(fields, data))
        for ax, (plot_field, (vertical_transform, vertical_transform_name)) in zip(axarr, plot_fields):
            ax.plot(data_dict["SecondsElapsed"],
                    [vertical_transform(el) for el in data_dict[plot_field]],
                    color=color,
                    alpha=0.6)
            ax.set_ylabel(vertical_transform_name + ' ' + plot_field)

ax.set_xlabel('SecondsElapsed')
fig.legend(handles=handles, loc="lower center", ncol=2)
fig.suptitle(prefix)
fig.tight_layout(rect=[0, 0.02, 1, 0.98])
fig.savefig(prefix + "accum_comparison.pdf", bbox_inches="tight")
