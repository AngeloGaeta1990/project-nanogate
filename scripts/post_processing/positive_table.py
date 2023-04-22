import pandas as pd
import numpy as np
import argh
import screed
import os
import math
import os.path
import ast
from pathlib import Path
from fnmatch import fnmatch
from collections import Counter

"""
class to handle nanogate result
"""


class JaccardResult():
    def __init__(self):
        # nanogate result_table
        self.table = None
        # Jaccard threshold value used to filter uncorrect matches part-read
        self.threshold = None
        # number of reads evaluted
        self.total_reads = None
        # number of constructs in the constructs library
        self.constructs_number = None
        # depht sequencing
        self.depht = None
        # run
        self.run = None
        # parts per construct
        self.parts = None
        # error rate used to get kmer length
        self.error_rate = None
        # when there is a mistake how many parts are wrong ? wrong_parts/wrong_reads
        self.error_within_error = None
        # number of reads matching the constructs length
        self.reads_analysed = None
        # reads matching the construct length / all reads
        self.good_read_ratio = None
        # percentage of constructs with at least 10 reads
        self.reads_10_constructs = None
        # number of reads whose number of expected parts match with the number of assigned parts
        self.positives = None
        # number of reads with a correct calling
        self.true_positives = None
        # number of reads matching a different construct in the library
        self.false_positives = None
        # number of reads not matching any constructs in the library
        self.false_negatives = None
        #precision
        self.precision = None
        # recall
        self.recall = None
        # number of parts called incorrectly
        self.misscalling = None
        # expected parts table related to the analysed library
        self.expected_parts_object = None
        # filename of naogate result
        self.filename = None
        #threshold_0
        self.threshold_0_result = None
        # ratio between constructs analysed correctly and constructs ordered correctly
        self.po_ratio = None
        # parts order positive
        self.po_positives = None
        # parts order true positives
        self.po_true_positives = None
        # number of reads matching a different construct in the library
        self.po_false_positives = None
        # number of reads not matching any constructs in the library
        self.po_false_negatives = None
        # precision
        self.po_precision = None
        # recall
        self.po_recall = None
        # number of parts called incorrectly
        self.po_misscalling = None
        # when there is a mistake how many parts are wrong ? wrong_parts/wrong_reads
        self.po_error_within_error = None

    def to_dict(self):
        return {
            'reads': self.total_reads,
            'parts': self.parts,
            'constructs number': self.constructs_number,
            'depth': self.depht,
            'jaccard threshold': self.threshold,
            'true positives': self.true_positives,
            'false positives': self.false_positives,
            'false negatives': self.false_negatives,
            'precision': self.precision,
            'recall': self.recall,
            'parts mismatch per read mismatch': self.error_within_error,
            'parts_order_ratio': self.po_ratio,
            'parts order true positives': self.po_true_positives,
            'parts order false positives': self.po_false_positives,
            'parts order false negatives': self.po_false_negatives,
            'parts order precision': self.po_precision,
            'parts order recall': self.po_recall,
            'parts order parts mismatch per read mismatch': self.po_error_within_error,
            'error_rate': self.error_rate,
            'reads analysed': self.reads_analysed,
            '10 reads constructs': self.reads_10_constructs,
            'good reads ratio': self.good_read_ratio
        }


class Threshold0():
    def __init__(self):
        #parts_per_construct
        self.parts_per_construct = None
        #library size
        self.constructs_library = None
        # error rate
        self.error_rate = None
        # depth
        self.depht = None
        # run
        self.run = None
        # table describing threshold 0 result
        self.table = None
        # filename
        self.filename = None
        #true positives
        self.true_positives_reads = None

    def build_threshold_0_object(self,filename):
        self.filename = filename
        try:
            self.table = pd.read_csv(self.filename)
            self.error_rate = self.table['error rate'].values[0]
            self.depht = self.table['depht'].values[0]
            self.threshold = self.table['jaccard threshold'].values[0]
            self.run = int(self.table["run"].values[0])
            self.constructs_library = self.table['constructs number'].values[0]
            self.parts_per_construct = self.table['parts_number'].values[0]
            # Gets stats
            self.find_true_positives_reads_id()
        except pd.errors.EmptyDataError:
            self.table = None


    def find_true_positives_reads_id(self):
        self.true_positives = 0
        pass_table = self.table[self.table['Positivity'] == "pass"]
        self.positives_reads = list(pass_table['Read id'].values)


#Class  to handle single Nanogate result row
class NanoResult:
    def __init__(self):
        #List of parts found by Nanogate
        self.nanogate_parts = []
        #List of parts expected
        self.expected_parts = []


"""
class describing constructs
"""
class Construct:
    def __init__(self):
        #construct id e.g. 1
        self.construct_id = None
        #size of the library
        self.construct_library = None
        #parts per constructs
        self.parts_per_construct = None
        #depth
        self.depth = None
        #error rate
        self.error_rate = None
        # run
        self.run = None
        #jaccard threshold
        self.threshold = None
        # table describing the construct, given by Processed_resut->pass->construct
        self.construct_table = None
        # processed nanogate result used to extract construct stats
        self.processed_table = None
        # number of reads for this construct, at this depth, error rate, run,
        self.reads_number = 0
        #positives
        self.positives = 0
        # true positives
        self.true_positives = 0
        # false positives
        self.false_positives = 0
        # false negatives
        self.false_negatives = 0
        # precision
        self.precision = 0
        # recall
        self.recall = 0
        # misscall per error
        self.error_within_error = 0
        # does the construct have more than 10 reads ? true or false
        self.reads_10 = None



    def to_dict(self):
        return {
            'construct id': self.construct_id,
            'constructs library': self.construct_library,
            'parts number': self.parts_per_construct,
            'depth': self.depth,
            'error_rate': self.error_rate,
            'run': self.run,
            'reads number': self.reads_number,
            '10 reads constructs': self.reads_10,
            'jaccard threshold': self.threshold,
            'true positives': self.true_positives,
            'false positives': self.false_positives,
            'false negatives': self.false_negatives,
            'precision': self.precision,
            'recall': self.recall,
            'parts mismatch per read mismatch': self.error_within_error
        }

#############################################################################################
#Reading parameters methods
#############################################################################################



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


def read_run(run):
    """
    convert string run "1" to integer run int 1
    :param run:
    :return:
    """
    run = int(run)
    return run

#############################################################################################################
#Jaccard result methods
###########################################################################################################
def get_jaccard_result(processed_tables, error_rate, depht,run,processed_table_main_folder,threshold_filename):
    """
    main method to fill all attributes

    """
    jaccard_results = []
    # Getting good read
    threshold_0_results = get_O_threshold_results(processed_table_main_folder,threshold_filename)
    for filename in processed_tables:
        try:
            table = pd.read_csv(filename)
            if error_rate in table['error rate'].values and depht in table['depht'].values and run in table['run'].values :
                jaccard_result = JaccardResult()
                jaccard_result.filename = filename
                jaccard_result.error_rate = table['error rate'].values[0]
                jaccard_result.depht = table["depht"].values[0]
                jaccard_result.run = table["run"].values[0]
                jaccard_result.threshold = table['jaccard threshold'].values[0]
                jaccard_result.constructs_number = table['constructs number'].values[0]
                jaccard_result.parts = table['parts_number'].values[0]
                jaccard_result.table = table
                # Gets stats
                positive_stats(jaccard_result,threshold_0_results)
                # Parts order confusion matrix stats
                parts_order_positive_stats(jaccard_result, threshold_0_results)
                #evaluate precision
                confusion_matrix_evalution(jaccard_result)
                # parts order confusion matrix
                parts_order_confusion_matrix_evalution(jaccard_result)
                # evalute reads quantity
                evalute_reads_analysed(jaccard_result,table)
                # count the amount of constructs with at least 10 reads
                count_10_reads_counstrcuts(jaccard_result, table)
                # Filling list of object
                jaccard_results.append(jaccard_result)
        except pd.errors.EmptyDataError:
            pass
    return jaccard_results


def confusion_matrix_evalution(jaccard_result):
    """
    evalute nanogate confusion matrix values
    :return:
    """
    try:
        jaccard_result.precision = jaccard_result.true_positives / (
                jaccard_result.true_positives + jaccard_result.false_positives)
    except ZeroDivisionError:
        jaccard_result.precision = None
        # evaluate bad calling per error
    try:
        jaccard_result.error_within_error = jaccard_result.misscalling
    except ZeroDivisionError:
        jaccard_result.error_within_error = None
        # evaluate reacall
    try:
        jaccard_result.recall = jaccard_result.true_positives / (
                jaccard_result.true_positives + jaccard_result.false_negatives)
    except ZeroDivisionError:
        jaccard_result.recall = None

    try:
        jaccard_result.po_ratio = jaccard_result.po_true_positives/jaccard_result.true_positives
    except ZeroDivisionError:
        jaccard_result.po_ratio = None

def find_theshold_false_negatives(jaccard_object, threshold_0_results):
    """
    add to the false negative counter, the false negative due to the threhsold
    :return:
       """
    threshold_negatives = None
    for threshold_0_result in threshold_0_results:
        if threshold_0_result.run == jaccard_object.run and \
            threshold_0_result.parts_per_construct == jaccard_object.parts and \
            threshold_0_result.constructs_library == jaccard_object.constructs_number and \
            threshold_0_result.error_rate == jaccard_object.error_rate and \
            threshold_0_result.depht == jaccard_object.depht:
            if jaccard_object.threshold > 0:
                threshold_negatives = 0
                fails_table = jaccard_object.table[jaccard_object.table["Positivity"] == "fail"]
                jaccard_object_fails = list(fails_table["Read id"].values)
                for fail_read in jaccard_object_fails:
                    if fail_read in threshold_0_result.positives_reads:
                        threshold_negatives += 1
            else:
                threshold_negatives = 0
    return threshold_negatives


def positive_stats( jaccard_result, threshold_0_results):
    """
    method to evaluted true and false positives
    :param processed_table:
    :param jaccard_result:
    :return:
    """
    jaccard_result.true_positives = 0
    jaccard_result.false_positives = 0
    jaccard_result.false_negatives = 0
    jaccard_result.misscalling = 0

    # positives-subtable

    positive_table = jaccard_result.table[jaccard_result.table['Positivity'] == "pass"]
    jaccard_result.positives = len(positive_table)

    # Iterating trough nanogate rows
    for index, row in positive_table.iterrows():
        expected_parts_list = row['expected_parts'].split()
        parts_list = row['Parts'].split()
        # Converting found and expected into an integer list
        int_parts_list = to_integer(parts_list)
        int_expected_parts_list = to_integer(expected_parts_list)
        # Measuring true and false positive
        if int_parts_list == int_expected_parts_list:
            jaccard_result.true_positives += 1
        else:
            status = evaluate_positive_or_negative(int_parts_list, positive_table)
            if status == "false positive":
                jaccard_result.false_positives += 1
            elif status == "false negative":
                jaccard_result.false_negatives += 1
            # Measuring bad parts per error
            jaccard_result.misscalling = evaluate_bad_parts_per_error(int_parts_list, int_expected_parts_list)
    # add false negative due to threshold
    jaccard_result.false_negatives += find_theshold_false_negatives(jaccard_result, threshold_0_results)



def evaluate_positive_or_negative(int_parts_list, positive_table):
    """
    if the result is not a true positive evalute if it is a false positive or false negative
    """
    status = None
    int_parts_set = set(int_parts_list)
    constructs_library = ast.literal_eval(positive_table["expected constructs library"].values[0])
    constructs_library_set = []
    for expected_construct_list in constructs_library:
        construct_set = set(expected_construct_list)
        constructs_library_set.append(construct_set)
    if int_parts_set in constructs_library_set:
             status = "false positive"
    else:
             status = "false negative"
    return status


def evalute_reads_analysed(jaccard_object, jaccard_table):
    """
    evalute the ratio pass reads/ all reads
    :param jaccard_object:
    :param jaccard_table:
    :return:
    """
    jaccard_object.total_reads = len(jaccard_table)
    positive_table = jaccard_table[jaccard_table['Positivity'] == "pass"]
    jaccard_object.reads_analysed = len(positive_table)
    jaccard_object.good_read_ratio = jaccard_object.reads_analysed/jaccard_object.total_reads


def count_10_reads_counstrcuts(jaccard_object, jaccard_table):
    """
    count in each nanogate result the amount of counstructs with at least 10 reads 
    """
    positive_table = jaccard_object.table[jaccard_object.table['Positivity'] == "pass"]
    pass_constructs_list =(list(positive_table["Construct"].values))
    constructs_frequency = dict(Counter(pass_constructs_list))
    frequencies = list(constructs_frequency.values())
    greater_10_frequencies = [x for x in frequencies if x >= 10]
    jaccard_object.reads_10_constructs = len(greater_10_frequencies)/jaccard_object.constructs_number



def evaluate_bad_parts_per_error(int_parts_list,int_expedced_parts_list):
    """
    given expected and called parts return how many bad called parts are present per error
    :return:
    """
    bad_part_calling = 0
    for part in int_parts_list:
        if part not in int_expedced_parts_list:
            bad_part_calling += 1
    bad_part_calling = bad_part_calling / len(int_expedced_parts_list)
    return bad_part_calling




def get_processed_tables_list(processed_table_main_folder, path):
    """
    query subfolders from root folder for list of processed files to analyse
    :param processed_table_main_folder:
    :return:
    """
    processed_tables_list = []
    pattern = path
    for path, subdirs, files in os.walk(processed_table_main_folder):
        for name in files:
            if fnmatch(name, pattern):
                file = (os.path.join(path, name))
                processed_tables_list.append(file)
    return processed_tables_list


def get_O_threshold_results(processed_table_main_folder,filename_standard):
    """
    for each nanogate result finds the related 0_trheshold result
    :return:
    """
    threshold_0_objects = []
    threshold_0_filenames = get_processed_tables_list(processed_table_main_folder,filename_standard)
    for threshold_0_filename in threshold_0_filenames:
            threshold_0 = Threshold0()
            threshold_0.build_threshold_0_object(threshold_0_filename)
            threshold_0_objects.append(threshold_0)
    return threshold_0_objects


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
    integer_list.sort()
    return integer_list

##########################################################################################
 #constructs methos
##########################################################################################
def load_constructs( processed_table, threshold_0_objects):
    """
    build constructs objects
    :param processed_table:
    :return:
    """
    constructs_objects = []
    pass_table = processed_table[processed_table['Positivity'] == "pass"]
    for id in pass_table['Construct'].unique():
        construct = Construct()
        construct.construct_id = id
        construct.processed_table = processed_table
        construct.construct_table = pass_table[pass_table['Construct']== construct.construct_id]
        construct.construct_library = construct.construct_table['constructs number'].values[0]
        construct.parts_per_construct = construct.construct_table['parts_number'].values[0]
        construct.depth = construct.construct_table['depht'].values[0]
        construct.error_rate = construct.construct_table['error rate'].values[0]
        construct.run = construct.construct_table['run'].values[0]
        construct.threshold = construct.construct_table['jaccard threshold'].values[0]
        construct.reads_number = len(construct.construct_table)
        find_10_reads_construct(construct)
        find_construct_positives(construct,threshold_0_objects)
        # measuring confusion matrix
        try:
            construct.precision = construct.true_positives / (
                    construct.true_positives + construct.false_positives)
        except ZeroDivisionError:
            construct.precision = None

        # evaluate recall
        try:
            construct.recall = construct.true_positives / (
                    construct.true_positives + construct.false_negatives)
        except ZeroDivisionError:
            construct.recall = None

        constructs_objects.append(construct)
    return constructs_objects


def find_construct_positives(construct_object,threshold_0_objects):
    """
    evaluate confusion matrix parameters for
    constructs
    :param construct_object:
    :return:
    """

    for index, row in construct_object.construct_table.iterrows():
        expected_parts_list = row['expected_parts'].split()
        parts_list = row['Parts'].split()
        # Converting found and expected into an integer list
        int_parts_list = to_integer(parts_list)
        int_expected_parts_list = to_integer(expected_parts_list)
        # Measuring true and false positive
        if int_parts_list == int_expected_parts_list:
            construct_object.true_positives += 1
        else:
            positive_table = construct_object.processed_table[construct_object.processed_table['Positivity'] == "pass"]
            status = evaluate_positive_or_negative(int_parts_list, positive_table)
            if status == "false positive":
                construct_object.false_positives += 1
            elif status == "false negative":
                construct_object.false_negatives += 1
            # Measuring bad parts per error
            construct_object.misscalling = evaluate_bad_parts_per_error(int_parts_list, int_expected_parts_list)
        construct_object.false_negatives += find_construct_theshold_false_negatives(construct_object, threshold_0_objects)

def  find_10_reads_construct(construct):
    """
    return True if the construct has more then 10 reads
    :param construct:
    :return:
    """
    if construct.reads_number > 10:
        construct.reads_10 = True
    else:
        construct.reads_10 = False

def find_construct_theshold_false_negatives(constructs_object, threshold_0_results):
    """
    add to the false negative counter, the false negative due to the threhsold
    :return:
       """
    threshold_negatives = None
    for threshold_0_result in threshold_0_results:
        if threshold_0_result.run == constructs_object.run and \
            threshold_0_result.parts_per_construct == constructs_object.parts_per_construct and \
            threshold_0_result.constructs_library == constructs_object.construct_library and \
            threshold_0_result.error_rate == constructs_object.error_rate and \
            threshold_0_result.depht == constructs_object.depth:
            if constructs_object.threshold > 0:
                threshold_negatives = 0
                fails_table = constructs_object.processed_table[constructs_object.processed_table["Positivity"] == "fail"]
                fails_table = fails_table[fails_table['Construct'] == constructs_object.construct_id]
                jaccard_object_fails = list(fails_table["Read id"].values)
                for fail_read in jaccard_object_fails:
                    if fail_read in threshold_0_result.positives_reads:
                        threshold_negatives += 1
            else:
                threshold_negatives = 0
    return threshold_negatives



def constucts_analysys(processed_tables,error_rate, depht, run,processed_table_main_folder, threshold_filename):
    """
    makes the constructs stats table
    :return:
    """
    threshold_0_results = get_O_threshold_results(processed_table_main_folder, threshold_filename)
    constructs_objects = []
    for table in processed_tables:
        try:
            table = pd.read_csv(table)
            if error_rate in table['error rate'].values and depht in table['depht'].values and run in table['run'].values:
                constructs_objects = constructs_objects + load_constructs(table,threshold_0_results)
        except pd.errors.EmptyDataError:
            pass
    return constructs_objects

######################################################################################################
# Parts order
########################################################################################
def parts_order_positive_stats( jaccard_result, threshold_0_results):
    """
    method to evaluted true and false positives
    :param processed_table:
    :param jaccard_result:
    :return:
    """
    jaccard_result.po_true_positives = 0
    jaccard_result.po_false_positives = 0
    jaccard_result.po_false_negatives = 0
    jaccard_result.po_misscalling = 0

    # positives-subtable

    positive_table = jaccard_result.table[jaccard_result.table['Positivity'] == "pass"]
    jaccard_result.positives = len(positive_table)

    # Iterating trough nanogate rows
    for index, row in positive_table.iterrows():
        expected_parts_list = row['expected_parts'].split()
        if type(row['Parts order']) == float:
            parts_list = []
        else:
            parts_list = row['Parts order'].split()
        # Converting found and expected into an integer list
        int_parts_list = to_integer_no_sorting(parts_list)
        int_expected_parts_list = to_integer_no_sorting(expected_parts_list)
        # Measuring true and false positive
        if int_parts_list == int_expected_parts_list:
            jaccard_result.po_true_positives += 1
        else:
            status = evaluate_positive_or_negative_parts_order(int_parts_list, positive_table)
            if status == "false positive":
                jaccard_result.po_false_positives += 1
            elif status == "false negative":
                jaccard_result.po_false_negatives += 1
            # Measuring bad parts per error
            jaccard_result.po_misscalling = evaluate_sorted_bad_parts_per_error(int_parts_list, int_expected_parts_list)


def to_integer_no_sorting(stringed_number_list):
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


def evaluate_positive_or_negative_parts_order(int_parts_list, positive_table):
    """
    if the result is not a true positive evalute if it is a false positive or false negative
    """
    status = None
    constructs_library = ast.literal_eval(positive_table["expected constructs library"].values[0])
    if int_parts_list in constructs_library:
             status = "false positive"
    else:
             status = "false negative"
    return status


def parts_order_confusion_matrix_evalution(jaccard_result):
    """
    evalute nanogate confusion matrix values for parts order
    :return:
    """
    try:
        jaccard_result.po_precision = jaccard_result.po_true_positives / (
                jaccard_result.po_true_positives + jaccard_result.po_false_positives)
    except ZeroDivisionError:
        jaccard_result.po_precision = None
        # evaluate bad calling per error
    try:
        jaccard_result.po_error_within_error = jaccard_result.po_misscalling
    except ZeroDivisionError:
        jaccard_result.po_error_within_error = None
        # evaluate recall
    try:
        jaccard_result.po_recall = jaccard_result.po_true_positives / (
                jaccard_result.po_true_positives + jaccard_result.po_false_negatives)
    except ZeroDivisionError:
        jaccard_result.po_recall = None


def evaluate_sorted_bad_parts_per_error(int_parts_list,int_expedced_parts_list):
    """
    given expected and called parts return how many bad called parts are present per error
    :return:
    """
    bad_part_calling = 0
    shortest_list =(min(len(int_parts_list),len(int_expedced_parts_list)))
    indexes = list(range(0, (shortest_list-1)))
    for index in indexes:
        if int_parts_list[index] != int_expedced_parts_list[index]:
            bad_part_calling += 1
    bad_part_calling = bad_part_calling / len(int_expedced_parts_list)
    return bad_part_calling
########################################################################################
#writing outputs
########################################################################################
def write_output_csv(jaccard_results: list, output: str):
    """
    write output table
    """
    jaccards_dict = []
    for jaccard in jaccard_results:
        jaccard_pandas = jaccard.to_dict()
        jaccards_dict.append(jaccard_pandas)
    output_table = pd.DataFrame(jaccards_dict)
    output_table.to_csv(output)

def write_construcs_csv(construcs_objects: list, output: str):
    """
    write output table
    """
    constructs_dict = []
    for construct in construcs_objects:
        construct_pandas = construct.to_dict()
        constructs_dict.append(construct_pandas)
    output_table = pd.DataFrame(constructs_dict)
    output_table.to_csv(output)


#####################################################################################################################
#MAIN
#######################################################################################################################à
def positive_analysis(jaccard_evaluation_table: 'output jaccard evaluation table filename',
                      constructs_stats_table: 'output constructs stats table filename',
                      processed_root_folder: 'processed table rootfolder',
                      error_rate = "010",
                      depht = "50x",
                      run = "1"):


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
    run = read_run(run)
    processed_tables = get_processed_tables_list(processed_root_folder,"processed_*.csv")
    jaccards = get_jaccard_result(processed_tables, error_rate, depht, run, processed_root_folder,"threshold_0_pro.csv")
    constructs = constucts_analysys(processed_tables,error_rate,depht,run, processed_root_folder,"threshold_0_pro.csv")
    write_output_csv(jaccards, jaccard_evaluation_table)
    write_construcs_csv(constructs, constructs_stats_table)


def main():
    """
    parsing arguments
    """

    argh.dispatch_commands([positive_analysis
                             ])


if __name__ == '__main__':
    main()
