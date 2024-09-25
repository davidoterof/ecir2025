#!/usr/bin/env python
import argparse
import glob
import matplotlib.pyplot as plt
import pandas as pd

fontsize = 15
plt.style.use('./style.mplstyle')
plt.rcParams.update({'figure.figsize': (13, 9)})
plt.rcParams.update({'font.size': fontsize})
plt.rcParams.update({'hatch.linewidth': 5})

samplesizes = [10, 30, 50]
systems = [3, 5, 10]

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
    plt.subplots_adjust(wspace=0.04, hspace=0.06)

    for ax, size in zip(axs[0], samplesizes):
        ax.set_title(rf"\texttt{{{size}}} topics", fontsize=fontsize)

    for ax in axs[:, 0]:
        ax.set_ylabel(r"FWER", fontsize=fontsize)

    for ax, system in zip(axs[:, 2], systems):
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(rf"$m$ = {system} systems ($k$ = {(system * (system - 1)) / 2:1.0f})", fontsize=fontsize)

    legg = False
    for i, system_number in enumerate(systems):
        for j, samplesize in enumerate(samplesizes):
            ax = axs[i][j]
            ax.set_ylim([0, 0.10])
            samplesize = samplesizes[j]

            plotting = pd.DataFrame(columns=["test", "error"])
            for alg, label in algs.items():

                sum = 0
                files = glob.glob(
                    f"{args.folder}/logreg_{args.collection}_alpha-0.05_{args.metric}_systems-{system_number}_topics-{samplesize}_sims_1/*-es-0.000.tsv")
                if len(files) == 0:
                    continue
                for file in files:
                    try:
                        df = pd.read_csv(file, sep="\t")
                    except pd.errors.EmptyDataError:
                        continue
                    svlue = df.iloc[[0]][alg].values[0]
                    value = float(svlue)
                    sum += value

                plotting.loc[len(plotting.index)] = [label, sum / len(files)]

            if plotting.empty:
                continue

            colors = ["#003A7D", "#008DFF", "#FF73B6", "#882255", "#D83034", "#FF9D3A"]
            ax.bar(plotting["test"].to_list(), plotting["error"].to_list(), label=plotting["test"].to_list(),
                   hatch=["", "", "", "", "", "//", "\\", "\\", "\\", "\\", "\\", "X"], color=colors,
                   edgecolor="white")
            ax.axhline(0.05, color='tab:gray', linestyle=(0, (5, 10)), linewidth=1, label=r"confidence level $\alpha$")

            ax.set_yticks([0.02, 0.04, 0.05, 0.06, 0.08, 0.10])
            ax.set_xticklabels([])
            if legg == False:
                legg = True
                plt.figlegend(loc='center', ncol=4, bbox_to_anchor=(0.50, 0.97), fontsize=fontsize)
            plt.savefig(f"{args.output}/fwer-{args.collection}-{args.metric}.pdf")
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', required=True, help="Folder to save the figure")
    parser.add_argument('-m', '--metric', required=True, help="MAP or NDCG")
    parser.add_argument('-f', '--folder', required=True, help="Folder with results of the simulation")
    parser.add_argument('-c', '--collection', required=True, help="name of the collection employed in the simulation")
    args = parser.parse_args()
    plot(args)
