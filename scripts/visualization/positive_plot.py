import seaborn as sns
from matplotlib import pyplot as plt
import argh
import pandas as pd
import itertools
import stracolors.calls
import numpy as np


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
    jaccard_plot = sns.catplot(x="depth", y="precision", data=table, col='parts',
                                kind='box', color=stracolors.calls.call_palette_paired()[4])
    jaccard_plot.set(ylabel="Precision")
    sns.set_context("paper")
    jaccard_plot.savefig(plot_name, orientation="portrait")


    # fig, axes = plt.subplots(ncols=3, sharex=True, sharey=True)
    # for ax, (n, grp) in zip(axes, table.groupby("parts")):
    #         sns.boxplot(x="depth", y="precision", data=grp, whis=np.inf, ax=ax)
    #         sns.swarmplot(x="depth", y="precision",  data=grp,
    #                       palette=["crimson", "indigo"], ax=ax)
    #         ax.set_title(n)
    #         ax.set(ylabel="Recall")
    #     # axes[-1].get_legend().remove()
    # sns.set_context("paper")
    # plt.savefig(plot_name, orientation="portrait")



def recall_plot(table, plot_name):
    """
    boxplot of the true precision
    """
    """
       barplot of the true positives
       """
    jaccard_plot = sns.catplot(x="depth", y="recall", data=table, col='parts',
                               kind='box', color=stracolors.calls.call_palette_paired()[4])
    jaccard_plot.set(ylabel="Recall")
    sns.set_context("paper")
    jaccard_plot.savefig(plot_name, orientation="portrait")

    # fig, axes = plt.subplots(ncols=3, sharex=True, sharey=True)
    # for ax, (n, grp) in zip(axes, table.groupby("parts")):
    #     sns.boxplot(x="depth", y="recall", data=grp, whis=np.inf, ax=ax)
    #     sns.swarmplot(x="depth", y="recall",  data=grp,
    #                   palette=["crimson", "indigo"], ax=ax)
    #     ax.set_title(n)
    #     ax.set(ylabel="Recall")
    # # axes[-1].get_legend().remove()
    #
    # sns.set_context("paper")
    # plt.savefig(plot_name, orientation="portrait")


def error_within_error_plot(table, plot_name):
    """
    barplot of the true positives
    """
    """
       barplot of the true positives
       """
    jaccard_plot = sns.catplot(x="jaccard threshold", y="parts mismatch per read mismatch", data=table, col='parts',
                               row="constructs number", kind='bar', color=stracolors.calls.call_palette_paired()[4])
    jaccard_plot.set(ylabel="Miss calling per error")
    sns.set_context("paper")
    jaccard_plot.savefig(plot_name, orientation="portrait")


def good_reads_ratio(table,plot_name):

    jaccard_plot = sns.catplot(x="depth", y="good reads ratio", data=table, col='parts',
                               kind='bar', color=stracolors.calls.call_palette_paired()[3])
    jaccard_plot.set(ylabel="Reads analysed")
    sns.set_context("paper")
    jaccard_plot.savefig(plot_name, orientation="portrait")



def reads_10_constructs(table,plot_name):
    
    jaccard_plot = sns.catplot(x="depth", y="10 reads constructs", data=table, col='parts',
                               kind='bar', hue='constructs number', color=stracolors.calls.call_palette_paired()[3])
    jaccard_plot.set(ylabel="Reads analysed")
    sns.set_context("paper")
    jaccard_plot.savefig(plot_name, orientation="portrait")


def parts_order_ratio(table, plot_name):
    jaccard_plot = sns.catplot(x="depth", y="parts_order_ratio", data=table, col='parts',
                               kind='box', color=stracolors.calls.call_palette_paired()[4])
    jaccard_plot.set(ylabel="Parts order ratio")
    sns.set_context("paper")
    jaccard_plot.savefig(plot_name, orientation="portrait")

# def bad_reads_ratio(table,plot_name):
#
#     table.plot(x="parts",  y=["random reads ratio", "chimeras ratio", "junk reads ratio", "bad reads ratio"], kind="bar")
#     plt.savefig(plot_name)

def part_order_precision_plot(table, plot_name):
    """
    boxplot of the parts order precision
    """
    jaccard_plot = sns.catplot(x="depth", y="parts order precision", data=table, col='parts',
                                kind='box', color=stracolors.calls.call_palette_paired()[4])
    jaccard_plot.set(ylabel="Parts order precision")
    sns.set_context("paper")
    jaccard_plot.savefig(plot_name, orientation="portrait")

def part_order_recall_plot(table, plot_name):
    """
    boxplot of the parts order recall
    """
    jaccard_plot = sns.catplot(x="depth", y="parts order recall", data=table, col='parts',
                                kind='box', color=stracolors.calls.call_palette_paired()[4])
    jaccard_plot.set(ylabel="Parts order recall")
    sns.set_context("paper")
    jaccard_plot.savefig(plot_name, orientation="portrait")





def positive_plot(table_filename: 'input table',
                  precision_plot_name: 'precision plot name',
                  recall_plot_name: 'recall plot name',
                  error_within_error_plot_name: 'error within error plot name',
                  good_reads_ratio_plot_name:'good reads ratio plot name',
                  reads_10_construct_plot_name:'reads 10 construct plot name',
                  parts_order_ratio_plot_name:'parts order ratio plot name',
                  parts_order_precision_plot_name:'parts order precision plot name',
                  parts_order_recall_plot_name:'parts order recall plot name'
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
    good_reads_ratio(table,good_reads_ratio_plot_name)
    reads_10_constructs(table, reads_10_construct_plot_name)
    parts_order_ratio(table, parts_order_ratio_plot_name)
    part_order_precision_plot(table,parts_order_precision_plot_name)
    part_order_recall_plot(table,parts_order_recall_plot_name)
    

def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([positive_plot
                            ])


if __name__ == '__main__':
    main()
