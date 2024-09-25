#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import cycle

fontsize = 15
plt.style.use('./style.mplstyle')
plt.rcParams.update({'figure.figsize': (13, 5)})
plt.rcParams.update({'font.size': fontsize})

styles = cycle(
    ['-', '-', '-', '-', '-', (0, (3, 1, 1, 1, 1, 1)), (0, (5, 10)), (0, (5, 10)), (0, (5, 10)), (0, (5, 10)),
     (0, (5, 10)), (0, (3, 1, 1, 1, 1, 1))])
markers = cycle(['o', 's', '^', 'x', '+', '*', 'o', 's', '^', 'x', '+', 'o'])
marklocations = cycle([7, 8, 9, 10])

algs = {
    "unadjusted wilcoxon": 5,
    "wilcoxon + bonferroni": 6,
    "wilcoxon + holm": 7,
    "wilcoxon + bh": 8,
    "wilcoxon + by": 9,
    "randomised TukeyHSD": 10,
    "unadjusted t-test": 0,
    "t-test + bonferroni": 1,
    "t-test + holm": 2,
    "t-test + bh": 3,
    "t-test + by": 4,
    "ANOVA + TukeyHSD": 11
}

alpha = 0.05


def plot(args):
    sample_sizes = np.arange(10, 101, 1)
    samples = np.arange(1, 2001, 1)

    gold_df = pd.read_csv(f"{args.folder}/gold.tsv", header=None, sep="\t")
    gold_pvalues = gold_df.iloc[0].to_numpy()[1:-1]
    h1_indexes_gold = set([x for x, e in enumerate(gold_pvalues) if e <= alpha])

    print("Plotting, this may take a while...")

    algs_plotting = {}
    for i, sample_size in enumerate(sample_sizes):
        for j, sample in enumerate(samples):
            file = f"num-queries-{sample_size}-sample-{sample}.tsv"
            df = pd.read_csv(f"{args.folder}/{file}", header=None, sep="\t")

            for label, index in algs.items():
                pvalues = df.iloc[index].to_numpy()[1:-1]

                if label not in algs_plotting:
                    algs_plotting[label] = np.zeros(len(sample_sizes))

                for index_gold in h1_indexes_gold:
                    if pvalues[index_gold] <= alpha:
                        algs_plotting[label][i] = algs_plotting[label][i] + 1

        colors = cycle(["#003A7D", "#008DFF", "#FF73B6", "#882255", "#D83034", "#FF9D3A"])
        for label, plotting in algs_plotting.items():
            algs_plotting[label][i] = algs_plotting[label][i] / len(samples)
            plt.plot(np.arange(10, sample_size + 1, 1), plotting[:i + 1] / len(h1_indexes_gold),
                     marker=next(markers), linestyle=next(styles),
                     label=label,
                     markevery=next(marklocations),
                     markersize=5, color=next(colors), linewidth=0.7)

        plt.savefig(f"{args.output}/millionquery-power.pdf")
        print(f"Plotted {i + 1} out of 100")

    plt.clf()

    colors = cycle(["#003A7D", "#008DFF", "#FF73B6", "#882255", "#D83034", "#FF9D3A"])
    for label, plotting in algs_plotting.items():
        plt.plot(np.arange(10, 101, 1), plotting / len(h1_indexes_gold), marker=next(markers),
                 linestyle=next(styles),
                 label=label,
                 markevery=next(marklocations),
                 markersize=5,
                 color=next(colors), linewidth=0.7)

    leg = plt.figlegend(loc='center', ncol=4, bbox_to_anchor=(0.50, 0.97), fontsize=fontsize)
    for line in leg.get_lines():
        line.set_linewidth(0.6)

    plt.ylabel("Average Power")
    plt.xlabel("Sample size")
    plt.xlim(0, 100)
    plt.xticks(np.arange(10, 101, 10))

    plt.savefig(f"{args.output}/millionquery-power.pdf")
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', required=True, help="Folder to save the figure")
    parser.add_argument('-f', '--folder', required=True, help="Folder with results of the experiment")
    args = parser.parse_args()
    plot(args)
