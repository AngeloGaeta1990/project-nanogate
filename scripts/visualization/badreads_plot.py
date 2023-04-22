import seaborn as sns
from matplotlib import pyplot as plt
import argh
import pandas as pd
import itertools
import stracolors.calls



def read_table(table_filename):
    """
    read table from .csv file
    """
    table = pd.read_csv(table_filename)
    return table




def false_positives(table,plot_name):
    if len(table)>1:
        jaccard_plot = sns.catplot(x="jaccard threshold", y="false positive ratio", data=table, hue="depth",
                                 kind='bar', color=stracolors.calls.call_palette_paired()[2])
        jaccard_plot.set(ylabel="Contaminants discovery")
        plt.ylim(0, 1)
        jaccard_plot.savefig(plot_name, orientation="portrait")
    else:
        jaccard_plot = open(plot_name, "w")





def badreads_plots(table_filename: 'input table',
                  false_positive_plot_name: 'true positive plot name',
                  ):
    """
    main function for positive plot, includes true positive, parts positives and error within error plot
    :param table_filename:

    """

    table = read_table(table_filename)
    false_positives(table, false_positive_plot_name)





def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([badreads_plots
                            ])


if __name__ == '__main__':
    main()
