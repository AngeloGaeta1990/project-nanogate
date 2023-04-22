import pandas as pd
import argh
import numpy as np
import math
import screed
import os
import logging
from fnmatch import fnmatch
from pathlib import Path

#commit for bumpversion

class SummaryStats:
    def __init__(self):
        # number of constructs in the library
        self.constructs_number = None
        # number of parts in a construct
        self.parts = None
        # depht used
        self.depht = None
        # reads produced
        self.reads = None
        # average reads analysed per construct
        self.average_reads = None
        # good reads per construct
        self.good_reads = None
        # jaccards coefficient of all parts belonging to a constructs
        self.jaccards = []
        # average of jaccards accepted
        self.average_jaccards = None
        # reads filtered out
        self.filtered_out = None
        # computational time
        self.time = None
        # run analysed
        self.run = None
        # error rate
        self.error_rate = None
        # k-mer length
        self.kmer_length = None
        # precision
        self.precision = None
        # recall
        self.recall = None

    def to_dict(self):
        return {
            'Run': self.run,
            'Constructs number': self.constructs_number,
            'Parts': self.parts,
            'Depth': self.depht,
            'Error rate': self.error_rate,
            'K-mer length': self.kmer_length,
            'Analysed reads': self.reads,
            'Discarded reads': self.filtered_out,
            # 'Average reads per construct': self.average_reads,
            'Good reads per construct': self.good_reads,
            'Average jaccard': self.average_jaccards,
            'Computing time': self.time,
            'Precision': self.precision,
            'Recall': self.recall

        }
class Construct:
    def __init__(self):
        # id related to the construct
        self.id = None
        # Sequence of the construct
        self.sequence = None
        # Lenght of the sequence
        self.sequence_length = None
        # full length reads
        self.perfect_reads = None
        # run related to the construct
        self.run = None
        # number of constructs in the original file
        self.constructs_number = None
        # number of parts in the original file
        self.parts = None

class Read:
    def __init__(self):
        # id related to the construct
        self.id = None
        # Sequence of the construct
        self.sequence = None
        # Lenght of the sequence
        self.sequence_length = None
        # run related to the construct
        self.run = None
        # number of constructs in the original file
        self.read_number = None
        # number of parts in the original file
        self.parts = None

def get_parts_reads(positive_table):
    """
    method to read from positive table parts per construct , and reads analysed per jaccard
    """
    summary_statistics = []
    positive_table = pd.read_csv(positive_table)
    for construct_number in positive_table["constructs number"].unique():
        for parts in positive_table["parts"].unique():
            summary_stats = SummaryStats()
            summary_stats.constructs_number = construct_number
            summary_stats.parts = parts
            part_set_subtable = positive_table.loc[
                (positive_table["constructs number"] == summary_stats.constructs_number)
                & (positive_table["parts"] == summary_stats.parts)]
            if len(part_set_subtable) > 0:
                summary_stats.reads = part_set_subtable["reads"].values[0]
                summary_stats.error_rate = part_set_subtable["error_rate"].values[0]
                # summary_stats.average_reads = part_set_subtable["average reads"].values[0]
                summary_stats.depht = part_set_subtable["depth"].values[0]
                summary_stats.good_reads = part_set_subtable["good reads ratio"].values[0]
                summary_stats.precision = part_set_subtable["precision"].values[0]
                summary_stats.recall = part_set_subtable["recall"].values[0]
            else:
                summary_stats.reads = None
                summary_stats.error_rate = None
                # summary_stats.average_reads = None
                summary_stats.depht = None
                summary_stats.good_reads = None

            summary_statistics.append(summary_stats)
    return summary_statistics


def get_jaccard(summary_statistics, processed_tables_filename, run):
    """
    method to get all jaccard values obtained
    :param summary_statistics:
    :param processed_tables:
    :return:
    """
    processed_tables =[]
    for summary_stats in summary_statistics:
        for filename in processed_tables_filename:
            try:
                processed_table = pd.read_csv(filename)
                processed_tables.append(processed_table)
            except pd.errors.EmptyDataError:
                    pass
        processed_table = pd.concat(processed_tables)

        part_set_subtable = processed_table.loc[
            (processed_table["constructs number"] == summary_stats.constructs_number)
            & (processed_table["parts_number"] == summary_stats.parts) & (
                    processed_table["error rate"] == summary_stats.error_rate)]
        if len(part_set_subtable)> 0:
            summary_stats.filtered_out = part_set_subtable["junk reads"].values[0]
            summary_stats.time = part_set_subtable["time"].values[0]
            summary_stats.kmer_length =part_set_subtable["kmer length"].values[0]
            summary_stats.run = int(run)
        else:
            summary_stats.filtered_out = None
            summary_stats.time = None
            summary_stats.kmer_length = None
            summary_stats.run = int(run)

        jaccards = part_set_subtable["Identity"].values
        for jaccard in jaccards:
            if type(jaccard) == str:
                jaccard_list = jaccard.split(", ")
                for jaccard in jaccard_list:
                    jaccard = float(jaccard)
                    if not math.isnan(jaccard):
                        if type(jaccard) == float:
                            summary_stats.jaccards.append(jaccard)
    for summary_stats in summary_statistics:
        summary_stats.average_jaccards = np.mean(summary_stats.jaccards)
    return summary_statistics


def get_processed_tables_list(processed_table_main_folder):
    """
    qurying list of nanogate processed tables
    :param processed_table_main_folder:
    :return:
    """
    processed_tables_list = []
    pattern = "processed_*.csv"
    for path, subdirs, files in os.walk(processed_table_main_folder):
        for name in files:
            if fnmatch(name, pattern):
                file = (os.path.join(path, name))
                processed_tables_list.append(file)
    return processed_tables_list




def write_csv_file(summary_statistics, output_summary_stats_csv):
    """
    write txt summary stats output file
    :param summary_statistics:
    :param output_summary_stats_txt:
    :return:
    """
    run_table = pd.DataFrame.from_records([run.to_dict() for run in summary_statistics])
    run_table.to_csv(output_summary_stats_csv)




def run_stats(positive_table: 'positive table',
              output_summary_stats_txt: 'output .txt filename',
              processed_tables: 'processed table',
              run=1,
              ):
    """
    Method to plot summary statistics after nanogate simulation
    :param positive_table: table with positive
    :param output_summary_stats_txt: output filename
    :param processed_tables: processed tables from nanogate
    :return:  a .txt file with summary statistics
    """
    summary_statistics = get_parts_reads(positive_table)
    processed_tables = get_processed_tables_list(processed_tables)
    summary_statistics = get_jaccard(summary_statistics, processed_tables, run)

    write_csv_file(summary_statistics, output_summary_stats_txt)


def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([run_stats
                            ])
if __name__ == '__main__':
    main()
