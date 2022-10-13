"""Build a comparison plot showing the average and range of tree accumulations during sample-optimize-merge pipeline.

Put replicate csv logfiles in directories named in `directories` list, and name the field to be plotted agains `Iterations` in `plot_field`."""

import pandas
from pathlib import Path
import matplotlib.pyplot as plt

plot_field = "NTrees"
directories = [
    'no-callback',
    'source-target-results',
    'target-node-results',
    'reject-all'
]

colors = [
    '#1b9e77',
    '#d95f02',
    '#7570b3',
    '#e7298a',
    '#66a61e',
    '#e6ab02',
]

cgen = iter(colors)
fig, ax = plt.subplots()
for directory in directories:
    color = next(cgen)
    p = Path(directory)
    dflist = []
    for logfile in p.glob('*.csv'):
        dflist.append(pandas.read_csv(logfile, sep='\t', header=0, index_col=False))
    df_concat = pandas.concat(dflist)
    by_row_df = df_concat.groupby(df_concat.index)
    mean_df = by_row_df.mean()
    max_df = by_row_df.max()
    min_df = by_row_df.min()
    ax.plot(mean_df['Iteration'], mean_df[plot_field], label=directory, color=color)
    ax.fill_between(min_df['Iteration'], min_df[plot_field], max_df[plot_field], alpha=0.2, color=color)

ax.legend(loc="upper left")
ax.set_ylabel(plot_field)
ax.set_xlabel('Iteration')
fig.savefig("accum_comparison.pdf")
