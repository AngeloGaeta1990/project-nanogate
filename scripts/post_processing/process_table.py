import pandas as pd
import argh
import json
import math
import os


def add_theshold(table, threshold):
    """
    Method to add a column to nanogate raw result
    specifing the jaccard threshold applied
    """
    table = pd.read_csv(table)
    threshold_column_length = len(table)
    threshold_float = 0
    threshold_float += float(threshold[0])
    threshold_float += (float(threshold[1]) / 10)
    if len(threshold) > 2:
        threshold_float += (float(threshold[2]) / 100)
    threshold_column = [threshold_float] * threshold_column_length
    table['jaccard threshold'] = threshold_column
    return table


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
    method to convert the depht integer
    :param depht:
    :return:
    """
    depht =depht.replace("x","")
    depht = int(depht)
    return depht


def read_run(run):
    """
    method to convert the depht integer
    :param depht:
    :return:
    """
    depht = int(run)
    return depht




def add_construcsts_number(table, construcst_number,error_rate,depht, run):
    """
    add to processed table constructs number, error rate, depht
    :param table:
    :param construcst_number:
    :param error_rate:
    :return:
    """
    construcst_number_length = len(table)
    construcst_number_column =[construcst_number] * construcst_number_length
    error_rate = read_error_rate(error_rate)
    error_rate_column = [error_rate] * construcst_number_length
    depht = read_depht(depht)
    depht_column =[depht] * construcst_number_length
    run = read_run(run)
    run_column = [run] * construcst_number_length
    table['constructs number'] = construcst_number_column
    table['error rate'] = error_rate_column
    table['depht'] = depht_column
    table['run'] = run_column
    table['constructs number'] = construcst_number_column
    table['error rate'] = error_rate_column
    table['depht'] = depht_column
    return table


def add_expected_result(threshold_table, verification_table):
    """
    Method to add a column to nanogate raw result including
    compares the verification table to the nanogate result table

    """
    verification_table = pd.read_csv(verification_table)
    expected_parts = []
    for index_threshold, row_threshold in threshold_table.iterrows():
        construct_id_test = row_threshold['Construct']
        if math.isnan(construct_id_test) == True:
            expected_parts.append(None)
        else:
            for index_verification, row_verification in verification_table.iterrows():
                if row_threshold['Construct'] == row_verification['Construct']:
                    expected_parts.append(row_verification['Parts'])

    threshold_table["expected_parts"] = expected_parts
    return threshold_table

def add_expected_constructs_library(threshold_table, verification_table):
     """
     add the list of all constructs to each row
     :param threshold_table:
     :param verification_table:
     :return:
     """
     constructs_library = []
     verification_table = pd.read_csv(verification_table)
     for part_row in verification_table["Parts"].values:
         construct = from_expected_string_list_to_integer_list(part_row)
         constructs_library.append(construct)

     table_len = len(threshold_table)
     constructs_library_column = [constructs_library] * table_len
     threshold_table["expected constructs library"] = constructs_library_column
     return threshold_table

def from_expected_string_list_to_integer_list( part_row):
    """
    turn list of expected parts in string to list of expected parts in integer
    :return:
    """
    part_row = part_row.replace(",", "")
    part_row = part_row.split(" ")
    cleaned_part_row = [x for x in part_row if x.isdigit()]
    int_parts_row = list(map(int, cleaned_part_row))
    int_parts_row.sort()
    return int_parts_row


def add_expected_parts_number(threshold_table):
    """
    :return:table including the number of parts in each construct
    """
    parts_numbers = []
    for index_parts, row_parts in threshold_table.iterrows():
        expected_parts = row_parts['expected_parts']
        if expected_parts == None:
            parts_numbers.append(0)
        else:
            parts_number = expected_parts.split(", ")
            # delete extra space in the list
            del parts_number[-1]
            parts_numbers.append(len(parts_number))
    threshold_table["parts_number"] = parts_numbers
    return threshold_table


def create_badreads_table(table):
    """
    include junk reads, random reads and chimera from nanogate results
    :param table:
    :return:
    """
    bad_reads_table = table[table['Quality'] != "standard"]
    return bad_reads_table


def filter_good_reads(table):
    """
    remove junk reads, random reads and chimera from nanogate results
    :return:
    """
    filtered_table = table[table['Quality'] == "standard"]
    return filtered_table


def add_from_jsonlog(table, jsonlog):
    """
    add time, removed parts to table
    :param table:
    :param jsonlog:
    :return:
    """
    with open(jsonlog) as json_file:
        json_data = json.load(json_file)

    time = json_data["time"]
    outliars = json_data["junk reads"]
    kmer_length = json_data["kmer_length"]
    column_length = len(table)
    time_column = [time] * column_length
    outliars_column = [outliars] * column_length
    kmer_length_column = [kmer_length]* column_length
    table["time"] = time_column
    table["junk reads"] = outliars_column
    table["kmer length"] = kmer_length_column
    return table





def write_processed_table(expected_table, output_filename):
    """
    method to turn the pandas table into a .csv file
    """
    expected_table.to_csv(output_filename)


def process_table(table: 'raw result table',
                  verification_table: 'verification table',
                  jsonlog: 'json log file',
                  output_processed_table: 'output processed table',
                  junk_read_table: 'junk read table',
                  threshold="01",
                  constructs_number = 50,
                  error_rate ="001",
                  depht ="50x",
                  run = "1"
                  ):
    """
    Takes as input nanogate raw result and add:
     1) threshold jaccard value used to discard containments
     2) Expected parts taken from the table generated by testbed_generator.py


    :param table: raw result table generated by nanogate
    :param threshold: threshold value to accept or discard parts containment in constructs
    :param verification_table: table generated by testbed generator, containts the real mapping construct-parts
    :return: processed table including thresholds and expected results
    """
    table_test =  pd.read_csv(table)
    if len(table_test) > 0:
        threshold_table = add_theshold(table, threshold)
        constructs_number_table = add_construcsts_number(threshold_table,constructs_number, error_rate, depht, run)
        expected_table = add_expected_result(constructs_number_table, verification_table)
        expected_table = add_expected_parts_number(expected_table)
        expected_table = add_expected_constructs_library(expected_table, verification_table)
        constructs_number_table = add_construcsts_number(threshold_table,constructs_number, error_rate, depht, run)
        expected_table = add_expected_result(constructs_number_table, verification_table)
        expected_table = add_expected_parts_number(expected_table)
        time_filteredout_table = add_from_jsonlog(expected_table, jsonlog)
        bad_reads_table = create_badreads_table(time_filteredout_table)
        filtered_table = filter_good_reads(time_filteredout_table)
        write_processed_table(filtered_table, output_processed_table)
        write_processed_table(bad_reads_table, junk_read_table)
    else:
        with open(output_processed_table, "w") as my_empty_csv:
                pass
        with open(junk_read_table, "w") as my_empty_csv:
                pass




def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([process_table
                            ])


if __name__ == '__main__':
    main()
