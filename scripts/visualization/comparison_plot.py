import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd
import argh

# positives_table = "true_positives_comparison.csv"
# time_table = "time_plot.cdv"




def precision_plot(table, plot_filename):
    jaccard_plot = sns.catplot(x="Parts number", y="precision", data=table, hue="Method",
                               kind='box', col="Coverage")
    jaccard_plot.set(ylabel="True positives")
    plt.ylim(0, 1)
    jaccard_plot.savefig(plot_filename, orientation="portrait")

def recall_plot(table, plot_filename):
    jaccard_plot = sns.catplot(x="Parts number", y="recall", data=table, hue="Method",
                               kind='box', col="Coverage")
    jaccard_plot.set(ylabel="True positives")
    plt.ylim(0, 1)
    jaccard_plot.savefig(plot_filename, orientation="portrait")


def time_plot(table, plot_filename):
    jaccard_plot = sns.catplot(x="Parts number", y="Real Time", data=table, hue="Method",
                               kind='bar', col="Coverage")
    jaccard_plot.set(ylabel="Time complexity")
    jaccard_plot.savefig(plot_filename, orientation="portrait")

def comparisons_plot(
                     #input
                     true_positive_comparison_table: 'true positive comparison table',
                     constructs_comparison_table: 'comparison table',
                     time_comparison_table: 'time comparison table',
                     # output
                     precision_5_parts_plot_filename:' 5 parts precision comparison plot filename',
                     recall_5_parts_plot_filename: ' 5 parts recall comparison plot filename',
                     precision_5_parts_constructs_plot_filename:'5 parts precsion plot filename',
                     recall_5_parts_constructs_plot_filename: '5 parts recall plot filename',
                     time_plot_5_parts_filename:'5 parts tme comparison plot filename',
                     precision_mn_plot_filename: 'minimap precision comparison plot filename',
                     recall_mn_plot_filename: 'minimap recall comparison plot filename',
                     time_mn_plot_filename: 'minimap nanogate time comparison plot filename'
                     ):
    """
    comparison plots
    :param true_positive_comparison_table:
    :param time_comparison_table:
    :param true_positives_plot:
    :param time_plot:
    :return:
    """
    # reading input table
    table_positives = pd.read_csv(true_positive_comparison_table)
    # Create subtable of 5 parts
    parts_number = table_positives["Parts number"].values
    min_parts_number = min(parts_number)
    parts_5_positives =table_positives.loc[(table_positives['Parts number'] == min_parts_number) ]
    # Create subtable Nanogate and Minimap2
    minimap2_nanogate_positives = table_positives.loc[(table_positives["Method"]!="minimap2 combinatorial")]


    #reading construct input tablr
    construct_table = pd.read_csv(constructs_comparison_table)
    parts_number = construct_table["Parts number"].values
    min_parts_number = min(parts_number)
    parts_5_constructs = construct_table.loc[(construct_table['Parts number'] == min_parts_number)]
    minimap2_nanogate_constructs = parts_5_constructs.loc[(parts_5_constructs["Method"]!="minimap2 combinatorial")]



    #Reading time table
    table_time = pd.read_csv(time_comparison_table)
    times_5_parts =table_time.loc[(table_time['Parts number'] == min_parts_number) ]
    time_nanogate_minimap2 =table_time.loc[(table_time["Method"]!="minimap2 combinatorial")]




    # 5 parts subtable plotting
    precision_plot(parts_5_positives, precision_5_parts_plot_filename)
    recall_plot(parts_5_positives, recall_5_parts_plot_filename)
    time_plot(times_5_parts, time_plot_5_parts_filename)
    #construct plotting
    precision_plot( minimap2_nanogate_constructs ,precision_5_parts_constructs_plot_filename)
    recall_plot(minimap2_nanogate_constructs, recall_5_parts_constructs_plot_filename)
    # Minimap2 Nanogate comparison
    precision_plot(minimap2_nanogate_positives, precision_mn_plot_filename)
    recall_plot(minimap2_nanogate_positives, recall_mn_plot_filename)
    time_plot(time_nanogate_minimap2, time_mn_plot_filename)





def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([comparisons_plot
                            ])


if __name__ == '__main__':
    main()
