#!/usr/bin/env python

import argparse
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import cycle

fontsize = 15
plt.style.use('./style.mplstyle')
plt.rcParams.update({'figure.figsize': (13, 10)})
plt.rcParams.update({'font.size': fontsize})

samplesizes = [10, 30, 50]
systems = [3, 5, 10]

styles = cycle(
    ['-', '-', '-', '-', '-', (0, (3, 1, 1, 1, 1, 1)), (0, (5, 10)), (0, (5, 10)), (0, (5, 10)), (0, (5, 10)),
     (0, (5, 10)), (0, (3, 1, 1, 1, 1, 1))])
markers = cycle(['o', 's', '^', 'x', '+', '*', 'o', 's', '^', 'x', '+', 'o'])
marklocations = cycle([7, 8, 9, 10])

algs = {
    "UnadjustedWilcoxon": "unadjusted wilcoxon",
    "BonferroniWilcoxon": "wilcoxon + bonferroni",
    "HolmWilcoxon": "wilcoxon + holm",
    "BHWilcoxon": "wilcoxon + bh",
    "BYWilcoxon": "wilcoxon + by",
    "Tukey": "randomised TukeyHSD",
    "UnadjustedTTest": "unadjusted t-test",
    "BonferroniTTest": "t-test + bonferroni",
    "HolmTTest": "t-test + holm",
    "BHTTest": "t-test + bh",
    "BYTTest": "t-test + by",
    "ANOVA": "ANOVA + TukeyHSD"
}


def plot(args):
    fig, axs = plt.subplots(nrows=3, ncols=3, sharex='all', sharey='all')
    plt.subplots_adjust(wspace=0.07, hspace=0.13)

    for ax, size in zip(axs[0], samplesizes):
        ax.set_title(rf"\texttt{{{size}}} topics", fontsize=fontsize)

    for ax in axs[:, 0]:
        ax.set_ylabel(r"Complete Power", fontsize=fontsize)

    for ax in axs[2]:
        ax.set_xlabel("Effect size", fontsize=fontsize)

    for ax, system in zip(axs[:, 2], systems):
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(rf"$m$ = {system} systems ($k$ = {(system * (system - 1)) / 2:1.0f})", fontsize=fontsize)

    legg = False
    changes = np.arange(0.005, 0.155, 0.005)
    for i, system_number in enumerate(systems):
        for j, samplesize in enumerate(samplesizes):
            colors = cycle(["#003A7D", "#008DFF", "#FF73B6", "#882255", "#D83034", "#FF9D3A"])
            ax = axs[i][j]
            ax.set_xticks(np.arange(0.05, 0.16, 0.05))
            ax.set_yticks(np.arange(0.2, 1.01, 0.2))
            ax.set_xlim([0, 0.15])
            ax.set_ylim([0, 1])

            curves = {}
            for alg, label in algs.items():
                curves[alg] = np.zeros(len(changes))

            for k, change in enumerate(changes):
                files = glob.glob(
                    f"{args.folder}/logreg_{args.collection}_alpha-0.05_{args.metric}_systems-{system_number}_topics-{samplesize}_sims_1000/input*-es-{change:1.3f}.tsv")
                if len(files) == 0:
                    continue
                for file in files:
                    df = pd.read_csv(file, sep="\t")
                    for alg, label in algs.items():
                        svlue = df.iloc[[0]][alg].values[0]
                        value = float(svlue)
                        curves[alg][k] += value

                for alg, label in algs.items():
                    curves[alg][k] = curves[alg][k] / len(files)

            for alg, label in algs.items():
                ax.plot(changes, curves[alg], marker=next(markers), linestyle=next(styles), label=f"{label}",
                        markevery=next(marklocations),
                        markersize=5, color=next(colors), linewidth=0.7)

            if legg == False:
                legg = True
                leg = plt.figlegend(loc='center', ncol=4, bbox_to_anchor=(0.50, 0.97), fontsize=fontsize)
                for line in leg.get_lines():
                    line.set_linewidth(0.6)

            plt.savefig(f"{args.output}/power-{args.collection}-{args.metric}.pdf")

    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', required=True, help="Folder to save the figure")
    parser.add_argument('-m', '--metric', required=True, help="MAP or NDCG")
    parser.add_argument('-f', '--folder', required=True, help="Folder with results of the simulation")
    parser.add_argument('-c', '--collection', required=True, help="name of the collection employed in the simulation")
    args = parser.parse_args()
    plot(args)
