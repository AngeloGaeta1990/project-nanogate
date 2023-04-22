import pandas as pd
import numpy as np
import argh
import screed
import os
import math
from fnmatch import fnmatch
"""
class to handle nanogate result
"""


class Badreads():
    def __init__(self):
        # parts
        # Jaccard threshold value used to filter uncorrect matches part-read
        self.threshold = None
        # number junk or random reads without assigned parts
        self.false_positives = 0
        # number junk or random reads with assigned parts
        self.false_positives_error= 0
        # number of reads evaluted
        self.reads_number = None
        # number of constructs in the constructs library
        self.constructs_number = None
        # depht sequencing
        self.depht = None
        # ratio correct results/number of reads
        self.false_positive_ratio = 0
        # parts per construct
        self.parts = None
        # error rate used to get kmer length
        self.error_rate = None
        # random reads
        self.random_reads = 0
        # random reads on all badreads
        self.random_reads_ratio = 0
        # chimera reads
        self.chimera_reads = 0
        # chimera reads on all badreads
        self.chimera_reads_ratio = 0
        # junk reads
        self.junk_reads = 0
        # junk reads on all badreads
        self.junk_reads_ratio = 0

    def to_dict(self):
        return {
            'reads': self.reads_number,
            'parts': self.parts,
            'constructs number': self.constructs_number,
            'depth': self.depht,
            'jaccard threshold': self.threshold,
            'false positive': self.false_positives,
            'false positive error': self.false_positives_error,
            'false positive ratio': self.false_positive_ratio,
            'error_rate': self.error_rate,
            'random reads ratio': self.random_reads_ratio,
            'chimeras ratio': self.chimera_reads_ratio,
            'junk reads ratio': self.junk_reads_ratio
        }

def read_error_rate(error_rate):
    """
    convert string "002" into float 0.02
    :param error rate:
    :return:
    """
    threshold_float = (float(error_rate) / 100)
    return threshold_float

def read_depht(depht):
    """
    convert string "002" into float 0.02
    :param error rate:
    :return:
    """
    depht = depht.replace("x","")
    depht = int(depht)
    return depht

def read_theshold(threshold):
    """
    convert string "02" into float 0.2
    :param threshold:
    :return:
    """
    threshold_float = 0
    threshold_float += float(threshold[0])
    threshold_float += (float(threshold[1])/10)
    if len(threshold)>2:
        threshold_float += (float(threshold[2])/100)
    return threshold_float




def get_badreads_stats(badreads_table, error_rate, depht ):
    """
    main method to fill all attributes
    """
    bad_reads = []

    for table in badreads_table:
        if os.stat(table).st_size > 0:
            table = pd.read_csv(table)
            if len(table) >0 :
                if table['error rate'].values[0] == error_rate and table['depht'].values[0]==depht:
                    bad_read = Badreads()
                    bad_read.error_rate = error_rate
                    bad_read.depht = depht
                    bad_read.threshold = table['jaccard threshold'].values[0]
                    bad_read.constructs_number = table['constructs number'].values[0]
                    bad_read.parts = 0
                    #Gets stats
                    false_positives_analysis(table, bad_read)
                    try:
                        bad_read.false_positive_ratio = bad_read.false_positives/(bad_read.junk_reads + bad_read.random_reads)
                    except ZeroDivisionError:
                        bad_read.false_positive_ratio = None
                    try:
                        bad_read.random_reads_ratio = bad_read.random_reads/(bad_read.junk_reads + bad_read.random_reads + bad_read.chimera_reads)
                    except ZeroDivisionError:
                        bad_read.random_reads_ratio = None
                    try:
                        bad_read.junk_reads_ratio = bad_read.junk_reads/(bad_read.junk_reads + bad_read.random_reads + bad_read.chimera_reads)
                    except ZeroDivisionError:
                        bad_read.junk_reads_ratio = None
                    try:
                        bad_read.chimera_reads_ratio =bad_read.chimera_reads/(bad_read.junk_reads + bad_read.random_reads + bad_read.chimera_reads)
                    except ZeroDivisionError:
                        bad_read.chimera_reads_ratio = None
                else:
                    bad_read = Badreads()
                    empty_table_evaluation(bad_read, error_rate, depht)
                #Filling list of object
                bad_reads.append(bad_read)
    return bad_reads


def empty_table_evaluation(bad_read, error_rate, depht ):
    """
    fill the badread object for empty tables
    :param bad_read:
    :return:
    """
    bad_read.error_rate = error_rate
    bad_read.depht = depht
    bad_read.threshold = None
    bad_read.constructs_number = None
    bad_read.parts = None
    bad_read.constructs_number = None
    bad_read.parts = None
    bad_read.false_positive_ratio = None
    bad_read.random_reads_ratio = None
    bad_read.junk_reads_ratio = None
    bad_read.chimera_reads_ratio = None




def false_positives_analysis(processed_table,bad_read):
    """
    counting chimeras, junk reads, random reads, and counting corretc false positive
    :return:
    """
    # Iterating trough nanogate rows
    for index, row in processed_table.iterrows():
        if "chimera" not in row['Quality']:
            if "junk" in row['Quality']:
                bad_read.junk_reads +=1
            elif "random" in row['Quality']:
                bad_read.random_reads +=1

            if type(row["Parts"]) == float:
                bad_read.false_positives += 1
            elif type(row["Parts"]) != float:
                bad_read.false_positives_error += 1
        else:
            bad_read.chimera_reads +=1



def get_processed_tables_list(processed_table_main_folder):
    """
    query subfolders from root folder for list of processed files to analyse
    :param processed_table_main_folder:
    :return:
    """
    processed_tables_list = []
    pattern = "bad_*.csv"
    for path, subdirs, files in os.walk(processed_table_main_folder):
        for name in files:
            if fnmatch(name, pattern):
                file = (os.path.join(path, name))
                processed_tables_list.append(file)
    return processed_tables_list


def to_integer(stringed_number_list):
    """
    :param stringed_number_list: a list of string in format  "n," e.g."4,"
    :return: a list of integers
    """
    integer_list = []
    for stringed_number in stringed_number_list:
        stringed_number = stringed_number.replace(",", "")
        number = int(stringed_number)
        integer_list.append(number)
    return integer_list


def write_output_csv(jaccard_results: list, output: str) -> None:
    """
    write output table
    """
    jaccards_dict = []
    for jaccard in jaccard_results:
        jaccard_pandas = jaccard.to_dict()
        jaccards_dict.append(jaccard_pandas)
    output_table = pd.DataFrame(jaccards_dict)
    output_table.to_csv(output)






def bad_reads_stats(bad_reads_table_filename: 'bad reads table filename',
                      badreads_tables_dir: 'processed table rootfolder',
                      error_rate = "010",
                      depht = "50x"):
    """
    Takes as input nanogate raw result and add:
     1) threshold jaccard value used to discard containments
     2) Expected parts taken from the table generated by testbed_generator.py
    :param table: raw result table generated by nanogate
    :param threshold: threshold value to accept or discard parts containment in constructs
    :param verification_table: table generated by testbed generator, contains the real mapping construct-parts
    :return: processed table including thresholds and expected results
    """

    error_rate = read_error_rate(error_rate)
    depht = read_depht(depht)
    processed_tables = get_processed_tables_list(badreads_tables_dir)
    badreads = get_badreads_stats(processed_tables, error_rate, depht)
    write_output_csv(badreads, bad_reads_table_filename)


def main():
    """
    parsing arguments
    """

    argh.dispatch_commands([bad_reads_stats
                       ])


if __name__ == '__main__':
    main()
