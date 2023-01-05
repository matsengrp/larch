"""Build a comparison plot showing the average and range of tree accumulations during sample-optimize-merge pipeline.

Put replicate csv logfiles in directories named in `directories` list, and name the field to be plotted agains `Iterations` in `plot_field`."""

import pandas
import csv
import statistics
from pathlib import Path
import matplotlib.pyplot as plt
import sys
sys.set_int_max_str_digits(10000)

prefix='20Dresults8/'
directories = [
    'option_',
    'option_--sample-best-tree',
    'option_--sample-uniform',
    'option_--sample-best-tree--sample-uniform',
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

# plot_field = "NTrees"
# directories = [
#     'accept-all',
#     'source-target-results',
#     'reject-all',
#     'default-callback',
#     'score-is-new',
#     'score-change',
# ]

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

for ax, (plot_field, (vertical_transform, vertical_transform_name)) in zip(axarr, plot_fields):
    cgen = iter(colors)
    for directory in directories:
        color = next(cgen)
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
                color=color)
        ax.fill_between(itlist,
                        minlist,
                        maxlist,
                        color=color,
                        alpha=0.2,)
        ax.set_ylabel(vertical_transform_name + ' ' + plot_field)

##pandas version
#for directory in directories:
#    color = next(cgen)
#    p = Path(prefix + directory)
#    dflist = []
#    for logfile in p.glob('*.csv'):
#        dflist.append(pandas.read_csv(logfile, sep='\t', header=0, index_col=False))
#    df_concat = pandas.concat(dflist)
#    by_row_df = df_concat.groupby(df_concat.index)
#    mean_df = by_row_df.mean()
#    max_df = by_row_df.max()
#    min_df = by_row_df.min()
#    ax.plot(mean_df['Iteration'], mean_df[plot_field], label=directory, color=color)
#    ax.fill_between(min_df['Iteration'], min_df[plot_field], max_df[plot_field], alpha=0.2, color=color)


ax.legend(loc="upper left")
ax.set_xlabel('Iteration')
fig.suptitle(prefix)
fig.savefig(prefix + "accum_comparison.pdf")
