import seaborn as sns
from matplotlib import pyplot as plt
import argh
import pandas as pd
import stracolors.calls






def time_plot(summary_table, plot_filename):

    time_graph = sns.catplot(x="Parts", y="Computing time", col="Constructs number", data=summary_table,
                    kind="bar", height=5, color=stracolors.calls.call_palette_paired()[1], facet_kws={"gridspec_kws": {"wspace": 0.1}},
                     aspect=1, ci="sd")
    time_graph.set(xlabel='Parts per construct', ylabel='Time(sec)')
    time_graph.savefig(plot_filename, orientation="portrait")



def kmer_length(summary_table, plot_filename):

    time_graph = sns.catplot(x="Error rate", y="K-mer length", col="Parts", data=summary_table, row ="Constructs number",
                    kind="bar", height=5, color=stracolors.calls.call_palette_paired()[2], facet_kws={"gridspec_kws": {"wspace": 0.1}},
                     aspect=1, ci="sd")
    time_graph.set(xlabel='Error rate', ylabel='K-mer length')
    time_graph.savefig(plot_filename, orientation="portrait")



def average_reads_per_construct(summary_table, plot_filename):

    time_graph = sns.catplot(x="Constructs number", y="Analysed reads", col="Parts", data=summary_table,
                    row= 'Depth',kind="bar", height=5, color=stracolors.calls.call_palette_paired()[3],
                             facet_kws={"gridspec_kws": {"wspace": 0.1}}, aspect=1, ci="sd")
    time_graph.set(xlabel='Constructs number', ylabel='Average reads per construct')
    time_graph.savefig(plot_filename, orientation="portrait")





def run_stats_plot(summary_table_filename: 'input table',
                  time_plot_filename : 'time plot filename',
                  kmer_length_plot_filename: 'kmer length plot filename',
                  average_reads_per_construct_plot_filename: 'average reads per construct plot filename',
):

    """
    Summary plots and the and of the experiment
    :param summary_table_filename:
    :param time_plot_filename:
    :param kmer_length_plot_filename:
    :param average_reads_per_construct_plot_filename:
    :return:
    """
    summary_table = pd.read_csv(summary_table_filename)
    time_plot(summary_table, time_plot_filename)
    kmer_length(summary_table, kmer_length_plot_filename)
    average_reads_per_construct(summary_table, average_reads_per_construct_plot_filename)


def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([run_stats_plot
                            ])


if __name__ == '__main__':
    main()
