import pandas as pd
import argh
import os
import numpy as np
import screed



"""
Class to analyse Minimap 2 results
"""


class Minimap2_Stats:
    def __init__(self):
        # table related to Minimap2 result
        self.table = None
        # number of constructs in the library
        self.constructs_library = None
        # number of parts per constructs
        self.parts_number = None
        # coverage eg. 5x, 100x
        self.coverage = None
        # Minimap2 real time
        self.real_time = None
        # Minimap2 CPU time
        self.cpu_time = None
        # Minimap2 positives
        self.positives = 0
        # Minimap2 true positives
        self.true_positives = 0
        # Minimap2 false positives
        self.false_positives = 0
        # Minimap2 false positives
        self.false_negatives = 0
        #Minimap2 Precision
        self.precision = 0
        #Minimap2 Recall
        self.recall = 0
        self.positives = None
        # Minimap2 true positives
        self.true_positives = None
        # Minimap2 total reads
        self.total_reads = None
        # Minimap2 alignment score
        self.alignment_score = None
        # Minimap2 mapping quality
        self.mapping_quality = None

    def to_dict(self):
        return {"Constructs library": self.constructs_library,
                "Parts number": self.parts_number,
                "Coverage": self.coverage,
                "Real Time": self.real_time,
                "CPU time": self.cpu_time,
                "True positives": self.true_positives,
                "recall": self.recall,
                "precision":self.precision,
                "Alignment score": self.alignment_score,
                "Mapping quality": self.mapping_quality
                }

"""
Class to analyse Minimap 2 results
"""
class Construct():
    def __init__(self):
        # construct id e.g. 1
        self.construct_id = None
        # size of the library
        self.construct_library = None
        # parts per constructs
        self.parts_per_construct = None
        # depth
        self.depth = None
        # run
        self.run = None
        # table describing the construct, given by Processed_resut->pass->construct
        self.construct_table = None
        # processed nanogate result used to extract construct stats
        self.minimap2_table = None
        # number of reads for this construct, at this depth, error rate, run,
        self.reads_number = 0
        # positives
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


    def to_dict(self):
        return {
            'construct id': self.construct_id,
            'Constructs library': self.construct_library,
            'Parts number': self.parts_per_construct,
            'Coverage': self.depth,
            'reads number': self.reads_number,
            'true positives': self.true_positives,
            'false positives': self.false_positives,
            'false negatives': self.false_negatives,
            'precision': self.precision,
            'recall': self.recall
        }


def evaluate_minimap2(output_filename: 'output filename',
                      constructs_filename: 'construct flename',
                      *processed_minimap2_tables: 'processed minimap2 tables'):
    """
    Main methods loop troughout all minimap2 processed tables to evalute stats
    :param output_filename:
    :param processed_minimap2_tables:
    :return:
    """
    minimap2_stats_objects = []
    constructs_objects = []
    # For each processed table get stats
    for processed_table in processed_minimap2_tables:
        if os.stat(processed_table).st_size > 0:

            table = pd.read_csv(processed_table)
            print(table)
            minimap2_stat = Minimap2_Stats()
            minimap2_stat.table = table
            # reads analysed
            minimap2_stat.total_reads = len(table)
            # Constructs library

            minimap2_stat.constructs_library = table["constructs number"].values[0]
            # Parts per contruct
            minimap2_stat.parts_number = table["parts_number"].values[0]
            # Coverage
            minimap2_stat.coverage = table["depht"].values[0]
            # Real time
            minimap2_stat.real_time = table["real time"].values[0]
            # CPU time
            minimap2_stat.cpu_time = table["CPU time"].values[0]
            #Alignment score
            minimap2_stat.alignment_score = np.mean(table["alignment score"].values)
            #Mapping quality
            minimap2_stat.mapping_quality = np.mean(table["mapping quality"].values)
            evalute_confusion_matrix_per_result(minimap2_stat)
            minimap2_stats_objects.append(minimap2_stat)
            #############
            #evaluate construct
            construct = minimap2_construct_analysis(table)
            constructs_objects.append(construct)
        # writing output
    write_output_table(minimap2_stats_objects, output_filename)
    write_construct_table(constructs_objects,constructs_filename)


def evalute_confusion_matrix_per_result(minimap2_object):
    """
    evalute the confusion matrix parameter per result
    """
    # Positives
    constructs_list = list(minimap2_object.table["expected constructs"].values)
    for index, row in minimap2_object.table.iterrows():
        if row["construct name"] == str(row["expected constructs"]):
            minimap2_object.true_positives += 1
        # True positives
        else:
            if row["construct name"] in constructs_list:
                minimap2_object.false_positives +=1
            else:
                minimap2_object.false_negatives +=1

    try:
        minimap2_object.precision = minimap2_object.true_positives / (minimap2_object.true_positives + minimap2_object.false_positives)
    except ZeroDivisionError:
        minimap2_object.precision = None
    try:
        minimap2_object.recall = minimap2_object.true_positives / (minimap2_object.true_positives + minimap2_object.false_negatives)
    except ZeroDivisionError:
        minimap2_object.recall = None

###################################################################################################à

#Construct analysis
###############################################################################################
def minimap2_construct_analysis(table):
    """"
    build construct object
    """
    for id in table['construct name'].unique():
        construct = Construct()
        construct.construct_id = id
        construct.minimap2_table = table
        construct.construct_table = table[table['construct name'] == construct.construct_id]
        construct.construct_library = construct.construct_table['constructs number'].values[0]
        construct.parts_per_construct = construct.construct_table['parts_number'].values[0]
        construct.depth = construct.construct_table['depht'].values[0]
        construct.reads_number = len(construct.construct_table)
        evaluate_constuct_confusion(construct)
        return construct


def evaluate_constuct_confusion(constuct_object):
    """
    evalute construct confusion matrix values
    :param constuct_object:
    :return:
    """
    # Positives
    constructs_list = list(constuct_object.minimap2_table["expected constructs"].values)
    for index, row in constuct_object.construct_table.iterrows():
        if row["construct name"] == str(row["expected constructs"]):
            constuct_object.true_positives += 1
        # True positives
        else:
            if row["construct name"] in constructs_list:
                constuct_object.false_positives +=1
            else:
                constuct_object.false_negatives +=1
    try:
        constuct_object.precision = constuct_object.true_positives / (constuct_object.true_positives + constuct_object.false_positives)
    except ZeroDivisionError:
        constuct_object.precision = None
    try:
        constuct_object.recall = constuct_object.true_positives / (constuct_object.true_positives + constuct_object.false_negatives)
    except:
        constuct_object.recall = None

###############################################################################################################


#writing output

#######################################################################################à

def write_output_table(minimap2_stats_objects,  output_filename):
    """
    Writes output table given minimap2 list objects and output filename
    :param minimap2_stats:
    :param output_filename:
    :return:
    """

    minimap_2_stat_table = pd.DataFrame.from_records([s.to_dict() for s in minimap2_stats_objects])
    minimap_2_stat_table.to_csv(output_filename)

def write_construct_table(minimap2_constructs_objects,  output_constructs_filename):
    """
    Writes output table given minimap2 list objects and output filename
    :param minimap2_stats:
    :param output_filename:
    :return:
    """

    minimap_2_stat_table = pd.DataFrame.from_records([s.to_dict() for s in minimap2_constructs_objects])
    minimap_2_stat_table.to_csv(output_constructs_filename)

#Testing input
# evaluate_minimap2("testing_output.csv", "minimap2_processed_alignment_table.csv")
def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([evaluate_minimap2
                            ])


if __name__ == '__main__':
    main()
