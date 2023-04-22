import seaborn as sns
from matplotlib import pyplot as plt
import argh
import pandas as pd
import itertools
import numpy as np
import stracolors.calls


def read_table(table_filename):
    """
    read table from .csv file
    """
    table = pd.read_csv(table_filename)
    return table


def precision_plot(table, plot_name):
    """
    boxplot of the true precision
    """
    if len(table)>1:
        jaccard_plot = sns.catplot(x="depth", y="precision", data=table,
                                   col='parts number', kind='box',
                                   color=stracolors.calls.call_palette_paired()[4])
        jaccard_plot.set(ylabel="Precision")
        sns.set_context("paper")
        jaccard_plot.savefig(plot_name, orientation="portrait")

        #Box plot with dots
        # fig, axes = plt.subplots(ncols=3, sharex=True, sharey=True)
        # for ax, (n, grp) in zip(axes, table.groupby("parts number")):
        #     sns.boxplot(x="depth", y="precision", data=grp, whis=np.inf, ax=ax)
        #     sns.swarmplot(x="depth", y="precision", data=grp,
        #                   palette=["crimson", "indigo"], ax=ax)
        #     ax.set_title(n)
        #     ax.set(ylabel="Precision")
        #     # axes[-1].get_legend().remove()
        # plt.show()
        #
        # sns.set_context("paper")
        # plt.savefig(plot_name, orientation="portrait")
    else:
        plot = open(plot_name, "w")





def recall_plot(table, plot_name):
    """
    boxplot of the true precision
    """
    if len(table) > 1:
        jaccard_plot = sns.catplot(x="depth", y="recall", data=table,
                                   col='parts number', kind='box',
                                   color=stracolors.calls.call_palette_paired()[4])
        jaccard_plot.set(ylabel="Recall")
        sns.set_context("paper")
        jaccard_plot.savefig(plot_name, orientation="portrait")

        #Plot with dots code
        # fig, axes = plt.subplots(ncols=3, sharex=True, sharey=True)
        #
        # for ax, (n, grp) in zip(axes, table.groupby("parts number")):
        #     sns.boxplot(x="depth", y="recall", data=grp, whis=np.inf, ax=ax)
        #     sns.swarmplot(x="depth", y="recall",  data=grp,
        #                   palette=["crimson", "indigo"], ax=ax)
        #     ax.set_title(n)
        #     ax.set(ylabel="Recall")
        # # axes[-1].get_legend().remove()
        #
        # sns.set_context("paper")
        # plt.savefig(plot_name, orientation="portrait")
    else:
        plot = open(plot_name, "w")


def error_within_error_plot(table, plot_name):
    """
    barplot of the true positives
    """

    if len(table) > 1:
        jaccard_plot = sns.catplot(x="jaccard threshold", y="parts mismatch per read mismatch", data=table, col='parts number',
                                   row="constructs library", kind='bar', color=stracolors.calls.call_palette_paired()[4])
        jaccard_plot.set(ylabel="Miss calling per error")
        sns.set_context("paper")
        jaccard_plot.savefig(plot_name, orientation="portrait")
    else:
        plot = open(plot_name, "w")


def reads_amount(table,plot_name):
   if len(table)>1:
        jaccard_plot = sns.catplot(x="depth", y="reads number", data=table, col='parts number',
                                   kind='bar', color=stracolors.calls.call_palette_paired()[3])
        jaccard_plot.set(ylabel="Reads analysed")
        sns.set_context("paper")
        jaccard_plot.savefig(plot_name, orientation="portrait")
   else:
       plot = open(plot_name, "w")


def constructs_10_reads_plot(table,plot_name):
    table = table[table['10 reads constructs'] == True]
    if len(table) > 1:
        jaccard_plot = sns.catplot(x="depth", y=None, data=table, col='parts number',
                                   kind='count', hue='constructs library', color=stracolors.calls.call_palette_paired()[3])
        jaccard_plot.set(ylabel="10 reads Constructs")
        sns.set_context("paper")
        jaccard_plot.savefig(plot_name, orientation="portrait")
    else:
        plot = open(plot_name, "w")


# def bad_reads_ratio(table,plot_name):
#
#     table.plot(x="parts",  y=["random reads ratio", "chimeras ratio", "junk reads ratio", "bad reads ratio"], kind="bar")
#     plt.savefig(plot_name)






def constructs_plot(table_filename: 'input table',
                  precision_plot_name: 'precision plot name',
                  recall_plot_name: 'recall plot name',
                  error_within_error_plot_name: 'error within error plot name',
                  good_reads_ratio_plot_name:'good reads ratio plot name',
                  construcst_10_reads: '10 reads constructs plot name'
                  ):
    """
    main function for positive plot, includes true positive, parts positives and error within error plot
    :param table_filename:
    :param true_positive_plot_name:
    :param true_positive_parts_plot_name:
    :param error_within_error_plot_name:
    :return:
    """

    table = read_table(table_filename)
    precision_plot(table, precision_plot_name)
    recall_plot(table, recall_plot_name)
    error_within_error_plot(table,error_within_error_plot_name)
    constructs_10_reads_plot(table,construcst_10_reads)
    reads_amount(table,good_reads_ratio_plot_name)

def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([constructs_plot
                            ])


if __name__ == '__main__':
    main()
