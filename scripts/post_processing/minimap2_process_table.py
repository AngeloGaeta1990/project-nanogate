import pandas as pd
import argh
import screed
import re
import os

def read_depht(depht):
    """
    method to convert the depht to integer
    :param depht:
    :return:
    """
    depht =depht.replace("x","")
    depht = int(depht)
    return depht



def add_one_value_columns(table, construcst_number, parts_number, depht, logfile):
    """
    add to processed table constructs number, error rate, depht
    :param table:
    :param construcst_number:
    :param error_rate:
    :return:
    """
    #Constructs number
    if table is not None:
        table_length = len(table)
        construcst_number_column =[construcst_number] * table_length

        # # Run
        # run = read_run(run)
        # run_column = [run] * table_length

        #Depht
        depht = read_depht(depht)
        depht_column = [depht] * table_length

        #Parts number
        parts_number_column = [parts_number] * table_length

        #Time
        time = get_time(logfile)
        if time != None:
            real_time_column = [time["Real time"]] * table_length
            cpu_time_column = [time["CPU time"]] * table_length

        #Adding 1 value columns
        table['constructs number'] = construcst_number_column
        # table['run'] = run_column
        table['depht'] = depht_column
        table['parts_number'] = parts_number_column
        table['real time'] = real_time_column
        table['CPU time'] = cpu_time_column
        return table
    else:
        return None


def read_reads_file(input_reads):
    """
    read reads file and return a dictionary reads id : expected parts
    :return:
    """
    read_construct_dict = {}
    for record in screed.open(input_reads):
        features = record.name
        splitname = features.split(" ")
        read = splitname[0]
        doublesplit = splitname[1].split(",")
        construct= int(doublesplit[0])
        read_construct_dict[read] = construct
    return read_construct_dict

def read_run(run):
    """
    read run turning into int
    """
    run = int(run)
    return run

def add_expected_result(minimap_table, read_construct_dict):
    """
    Method to add a column to minimap2 result specifying expected construct
    """
    if os.stat(minimap_table).st_size > 0:
        minimap_table = pd.read_csv(minimap_table)
        expected_constructs = []
        for index_minimap_table, row_minimap_table in minimap_table.iterrows():
                if row_minimap_table['read id'] in read_construct_dict.keys():
                    expected_constructs.append(read_construct_dict[row_minimap_table['read id']])
        minimap_table["expected constructs"] = expected_constructs
        return minimap_table
    else:
        return None

def get_time(logfile = None):
    """
    Get time from logfile
    """
    if logfile != None:
        time_dict = {}
        logfile = open(logfile, "r")
        for line in logfile.readlines():
            if "Real time" in line:
                times = re.findall(r'\d+', line)
                times = [int(x) for x in times]
                time_dict["Real time"] = times[0] + (times[1]/1000)
                time_dict["CPU time"] =times[2] +(times[3]/1000)
        return time_dict
    else:
        return None


def write_processed_table(expected_table, output_filename):
    """
    method to turn the pandas table into a .csv file
    """
    expected_table.to_csv(output_filename)


def process_table(minimap_result: 'raw result table',
                  input_reads: 'input reads',
                  output_processed_table: 'output processed table',
                  logfile = None,
                  constructs_number = 50,
                  parts_number = 20,
                  depht ="50x"
                  # run = "1"
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


    read_id_dict = read_reads_file(input_reads)
    expected_table = add_expected_result(minimap_result, read_id_dict)
    expanded_table = add_one_value_columns(expected_table, constructs_number, parts_number, depht, logfile)
    if expanded_table is not None:
        write_processed_table(expanded_table, output_processed_table)
    else:
        with open(output_processed_table, 'w') as fp:
            pass


def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([process_table
                            ])


if __name__ == '__main__':
<<<<<<< HEAD
    main()
=======
    main()
>>>>>>> db58617d792aaffb4a97fe8bb31582316a9320e7
