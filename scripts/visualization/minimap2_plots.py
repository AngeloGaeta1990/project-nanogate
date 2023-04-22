import seaborn as sns
from matplotlib import pyplot as plt
import argh
import numpy as np
import pandas as pd
import stracolors.calls




def precision_plot(table, plot_name):
    """
    precision boxplot
    """
    jaccard_plot = sns.catplot(x="Parts number", y="precision", data=table, hue='Coverage',
                             kind='box', color=stracolors.calls.call_palette_paired()[2])
    jaccard_plot.set(ylabel="Precision")
    jaccard_plot.savefig(plot_name, orientation="portrait")
    sns.set_context("paper")
    plt.savefig(plot_name, orientation="portrait")

    # fig, axes = plt.subplots(ncols=3, sharex=True, sharey=True)
    # for ax, (n, grp) in zip(axes, table.groupby("Parts number")):
    #     sns.boxplot(x="Coverage", y="precision", data=grp, whis=np.inf, ax=ax)
    #     sns.swarmplot(x="Coverage", y="precision", data=grp,
    #                   palette=["crimson", "indigo"], ax=ax)
    #     ax.set_title(n)
    #     ax.set(ylabel="Recall")
    # # axes[-1].get_legend().remove()
    # sns.set_context("paper")
    # plt.savefig(plot_name, orientation="portrait")

def recall_plot(table, plot_name):
    """
    boxplot of the recall
    """
    jaccard_plot = sns.catplot(x="Parts number", y="recall", data=table, hue='Coverage',
                             kind='box', color=stracolors.calls.call_palette_paired()[2])
    jaccard_plot.set(ylabel="Recall")
    jaccard_plot.savefig(plot_name, orientation="portrait")
    sns.set_context("paper")
    plt.savefig(plot_name, orientation="portrait")
    # fig, axes = plt.subplots(ncols=3, sharex=True, sharey=True)
    #
    # for ax, (n, grp) in zip(axes, table.groupby("Parts number")):
    #     sns.boxplot(x="Coverage", y="recall", data=grp, whis=np.inf, ax=ax)
    #     sns.swarmplot(x="Coverage", y="recall",  data=grp,
    #                   palette=["crimson", "indigo"], ax=ax)
    #     ax.set_title(n)
    #     ax.set(ylabel="Recall")
    # # axes[-1].get_legend().remove()




def alignment_score_plot(table, plot_name):
    """
    barplot of the true positives
    """
    jaccard_plot = sns.catplot(x="Parts number", y="Alignment score", data=table, hue='Coverage',
                             kind='bar', color=stracolors.calls.call_palette_paired()[3])
    jaccard_plot.set(ylabel="Alignment score")
    jaccard_plot.savefig(plot_name, orientation="portrait")

def mapping_quality_plot(table, plot_name):
    """
    barplot of the true positives
    """
    jaccard_plot = sns.catplot(x="Parts number", y="Mapping quality", data=table, hue='Coverage',
                             kind='bar', color=stracolors.calls.call_palette_paired()[2])
    jaccard_plot.set(ylabel="Mapping quality")
    jaccard_plot.savefig(plot_name, orientation="portrait")

def real_time_plot(table, plot_name):
    """
    barplot of real time
    """
    time_graph = sns.catplot(x="Parts number", y="Real Time", col="Constructs library", data=table,
                             kind="bar", height=5, color=stracolors.calls.call_palette_paired()[1],
                             facet_kws={"gridspec_kws": {"wspace": 0.1}},
                             aspect=1, ci="sd")
    time_graph.set(xlabel='Parts per construct', ylabel='Time(sec)')
    time_graph.savefig(plot_name, orientation="portrait")

def cpu_time_plot(table, plot_name):
    """
    barplot of CPU time
    """
    time_graph = sns.catplot(x="Parts number", y="CPU time", col="Constructs library", data=table,
                             kind="bar", height=5, color=stracolors.calls.call_palette_paired()[1],
                             facet_kws={"gridspec_kws": {"wspace": 0.1}},
                             aspect=1, ci="sd")
    time_graph.set(xlabel='Parts per construct', ylabel='Time(sec)')
    time_graph.savefig(plot_name, orientation="portrait")

def minimap2_plots(minimap2_stats_table: 'minimap2 stats',
                   minimap2_precision_filename: 'minimap2 precision filename',
                   minimap2_recall_filename :'minimap2 recall filename',
                   minimap2_mapping_quality_filename: 'minimap2 mapping quality filename',
                   minimap2_alignment_score_filename: 'minimap2 alignment score filename',
                   minimap2_realtime_filename: 'minimap2 real time filename',
                   minimap2_cputime_filename:  'minimap2 cpu time filename'
                   ):

    table = pd.read_csv(minimap2_stats_table)
    # preicision plot
    precision_plot(table, minimap2_precision_filename)
    # recall plot
    recall_plot(table,minimap2_recall_filename)
    #Mapping quality plot
    mapping_quality_plot(table,minimap2_mapping_quality_filename)
    #Alignment score time
    alignment_score_plot(table,minimap2_alignment_score_filename)
    # Real time plot
    real_time_plot(table, minimap2_realtime_filename)
    # CPU time plot
    cpu_time_plot(table, minimap2_cputime_filename)



def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([minimap2_plots
                            ])


if __name__ == '__main__':
    main()
