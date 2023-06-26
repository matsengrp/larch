"""Build a comparison plot showing the average and range of tree accumulations during sample-optimize-merge pipeline.

Put replicate csv logfiles in directories named in `directories` list, and name the field to be plotted agains `Iterations` in `plot_field`."""

import pandas
import csv
import statistics
from pathlib import Path
import matplotlib.pyplot as plt
import sys
from itertools import product

all_data = {}
kwnames = ["input_dag", "Iterations", "replicates", "log_folder_name", "optionslist"]

dir_name = sys.argv[1]
with open(dir_name + "/python_params.txt", "r") as f:
    data = [x.strip() for x in f.readlines()]
    all_data = dict(zip(kwnames[:4], data[:4]))

if len(data) > 4 and data[4].strip() == "+":
    opts = ("...".join(data[5:])).split("+")
    if len(opts) > 1:
        all_data["options"] = ["-".join([z for z in y if z != '']) for y in product(*[x.split('...') for x in opts])]
    else:
        all_data["options"] = opts[0].split("...")
else:
    all_data["options"] = ['']

prefix=dir_name + '/' + all_data["log_folder_name"] + "/"

directories = [
    "option_" + x
    for x in all_data["options"]
]

logtransform = (lambda n: statistics.log(n, 10), "Log10")
notransform = (lambda n: n, "")

plot_fields = [("MaxParsimony", notransform),
               ("WorstParsimony", notransform),
               ("NTreesMaxParsimony", logtransform),
               ("NTrees", logtransform),
               ("NEdges", notransform),
               ("NNodes", notransform),
               ("SecondsElapsed", notransform),
               ]

markers = ["o", "v", ">", "<", "*", "+", "x", "1", "2", "3", "4",]

colors = [
    '#1b9e77',
    '#d95f02',
    '#7570b3',
    '#e7298a',
    '#66a61e',
    '#e6ab02',
    '#33cc33',
    '#6600ff',
    '#00ace6',
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

for ax, (plot_field, (vertical_transform, vertical_transform_name)) in zip(axarr, plot_fields):
    cgen = iter(colors)
    mgen = iter(markers)
    for directory in directories:
        try:
            color = next(cgen)
        except:
            cgen = iter(colors)
            color = next(cgen)
        try:
            marker = next(mgen)
        except:
            mgen = iter(markers)
            marker = next(mgen)

        p = Path(prefix + directory)
        data = [
            [[] for _ in range(iterations)]
            for _ in range(len(fields))
        ]
        for logfile in p.glob('log*/logfile.csv'):
            with open(logfile, newline='') as csvfile:
                reader = csv.reader(csvfile, delimiter='\t')
                try:
                    next(reader)
                except StopIteration:
                    continue
                for iteration_idx, row in enumerate(reader):
                    for field_idx, dat in enumerate(row):
                        data[field_idx][iteration_idx].append(0 if dat == '' else int(dat))
        data_dict = dict(zip(fields, data))
        itlist = list(range(len([row[0] for row in data_dict['Iteration'] if len(row) > 0])))
        meanlist = [vertical_transform(mean(row)) for row in data_dict[plot_field] if len(row) > 0]
        minlist = [vertical_transform(min(row)) for row in data_dict[plot_field] if len(row) > 0]
        maxlist = [vertical_transform(max(row)) for row in data_dict[plot_field] if len(row) > 0]
        assert len(itlist) == len(maxlist) and len(maxlist) == len(minlist) and len(minlist) == len(meanlist)
        ax.plot(itlist,
                meanlist,
                label=directory,
                marker=marker,
                color=color)
        ax.fill_between(itlist,
                        minlist,
                        maxlist,
                        color=color,
                        alpha=0.2,)
        ax.set_ylabel(vertical_transform_name + ' ' + plot_field)


ax.legend(loc="lower center", bbox_to_anchor=(0.5,-.2-.08*len(directories)))
ax.set_xlabel('Iteration')
fig.suptitle("folder: " + prefix+"\n"+"start DAG: "+all_data["input_dag"]+"\n# replicates: %s"%all_data["replicates"])
fig.savefig(prefix + "accum_comparison.pdf")
