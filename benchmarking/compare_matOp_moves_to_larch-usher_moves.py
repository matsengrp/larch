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

fig, ax = plt.subplots(2, sharex=False, figsize=(7, 9))

for directory in directories:
    p = Path(prefix + directory)

    matOp_data = []
    larchUsher_data = []
    for logfile in p.glob('matOp_*.log'):
        with open(logfile, "r") as logfile:
            matOpMoves = []
            larchUsherMoves = []
            for line in logfile.readlines():
                if "Applying " in line:
                    numMoves = int(line.split("Applying ")[-1].split()[0])
                    if numMoves != 0:
                        if ("Larch-Usher" in line):
                             larchUsherMoves.append(numMoves)
                        else:
                             matOpMoves.append(numMoves)

            larchUsher_data.append(larchUsherMoves)
            matOp_data.append(matOpMoves)

    larchUsher_data = pandas.DataFrame(zip(*larchUsher_data))
    matOp_data = pandas.DataFrame(zip(*matOp_data))

    ax[0].plot(list(matOp_data.index), matOp_data.mean(axis=1))
    ax[0].fill_between(list(matOp_data.index), matOp_data.min(axis=1), matOp_data.max(axis=1), alpha=.2)

    ax[1].plot(list(larchUsher_data.index), larchUsher_data.mean(axis=1), label=directory)
    ax[1].fill_between(list(larchUsher_data.index), larchUsher_data.min(axis=1), larchUsher_data.max(axis=1), alpha=.2)
ax[0].set_title("matOptimize")
ax[1].set_title("larch-usher")
ax[0].set_ylabel("#moves applied")
ax[1].set_ylabel("#moves applied")
ax[1].set_xlabel('Iteration')

ax[1].legend(loc="lower center", bbox_to_anchor=(0.5,-.2-.08*len(directories)))
fig.suptitle("comparison of moves applied")
fig.savefig(prefix+"moves_applied.pdf", bbox_inches='tight')
