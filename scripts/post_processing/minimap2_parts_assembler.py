import screed
import re
import pandas as pd
import argh
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

"""
Abstract class contains elements shared by Parts and Constructs
"""


class DNAcomponent:
    def __init__(self):
        # integer representing the component
        self.id = None
        # DNA sequence
        self.sequence = None
        # Description  in input/output fasta
        self.description = None


"""
Class to handle to DNA parts
"""


class Part(DNAcomponent):
    def __init__(self):
        DNAcomponent.__init__(self)
        self.clone = None
        self.position = None

    # takes from description integer id representing the part
    def get_id(self):
        """
        get for sequence parts object
        """
        id = int(re.search(r'\d+', self.description).group())
        return id

    def get_position(self, position_parts_dict):
        """
        get position for parts objects
        """
        for key, value in position_parts_dict.items():
            if self.id in value:
                self.position = key

    # build parts object
    def build_part(self, record, position_parts_dict):
        """
        parts constructor
        """
        self.description = record.name
        self.sequence = record.sequence
        self.id = self.get_id()
        self.get_position(position_parts_dict)


class Construct(DNAcomponent):
    def __init__(self):
        DNAcomponent.__init__(self)
        # list of parts object in the construct
        self.parts = []
        self.parts_id = []
        self.parts_description = None

    def build_construct(self, combination, parts_list, construct_id):
        """
        construct constructor
        """
        self.get_parts_id(combination)
        self.get_parts(parts_list)
        self.get_sequence()
        self.get_id(construct_id)
        self.get_parts_description()

    def get_parts_id(self, combination):
        """
        get id of all parts in the construct
        """
        for part_id in combination:
            self.parts_id.append(part_id)

    def get_parts(self, parts_list):
        """
        get parts objects from parts object list
        """
        for part in parts_list:
            if part.id in self.parts_id:
                self.parts.append(part)

    def get_sequence(self):
        """
        build sequence from by summing parts sequences
        """
        self.parts.sort(key=lambda x: x.position, reverse=False)
        self.sequence = Seq("", generic_dna)
        for part in self.parts:
            self.sequence = self.sequence + part.sequence

    def get_id(self, construct_id):
        """
        get construct id
        """
        self.id = construct_id

    def get_parts_description(self):
        """
        get parts description to add to final fasta
        """
        parts_str = str(self.parts_id)
        parts_str = parts_str.replace("[", "")
        parts_str = parts_str.replace("]", "")
        parts_str = parts_str + ","
        self.parts_description = parts_str

    def to_dict(self):
        return {
            'Construct': self.id,
            'Parts': self.parts_description
        }


def build_position_parts_dictionary(input_constructs, input_parts_number):
    # reading from input constructs parts position
    postion_parts_dict = {}
    for part in range(input_parts_number):
        postion_parts_dict[part + 1] = []

    for record in screed.open(input_constructs):
        description_list = record.name.split(" ")
        description_list = description_list[4:]
        for part in range(len(description_list)):
            part_id = int(description_list[part].replace(",", ""))
            if part_id not in postion_parts_dict[part + 1]:
                postion_parts_dict[part + 1].append(part_id)

    return postion_parts_dict


# print(postion_parts_dict)
# function to prcombinations that contain
# one element from each of the given arrays
def get_combinations(arr):
    # number of arrays
    if len(arr) <= 10:
        n = len(arr)
        # to keep track of next element in each of the n arrays
        indices = [0 for i in range(n)]

        combinations = []
        while (1):

            # prcurrent combination
            new_combination = []
            for i in range(n):
                new_combination.append(arr[i][indices[i]])
                if len(new_combination) == n:
                    combinations.append(new_combination)
                    new_combination = []

            # find the rightmost array that has more elements left after the current element in that array
            next = n - 1
            while (next >= 0 and
                   (indices[next] + 1 >= len(arr[next]))):
                next -= 1

            # no such array is found so no more combinations left
            if (next < 0):
                return combinations
            # if found move to next element in that array
            indices[next] += 1

            # for all arrays to the right of this array current index again points to first element
            for i in range(next + 1, n):
                indices[i] = 0


def write_constructslibrary_fasta(constructs, filename):
    """
    Turns a list of constructs into a fasta file with all constructs
    """
    # writing a fasta file for constructs
    with open(filename, "w") as fastaseq:
        for construct in constructs:
            seq_record = SeqRecord(
                construct.sequence,
                id=str(construct.id),
                name="",
                description=">" + "Construct " + str(construct.id) + " " + "Parts " + str(
                    construct.parts_description) + "\n",
            )
            SeqIO.write(seq_record, fastaseq, "fasta")


def write_varification_csv(constructs, output_verifcation_csv):
    """
    method to write a csv file, with two columns:
    construct parts eg constuct =1 parts = 4 5 6
    """
    pandas_dict = []
    for construct in constructs:
        construct_dict = construct.to_dict()
        pandas_dict.append(construct_dict)
    verification_table = pd.DataFrame(pandas_dict)
    verification_table.to_csv(output_verifcation_csv)


# #####################################################################
# """
# DEFINING INPUTS/OUTPUTS
# """
# input_parts_number = 5
# input_parts = "20_parts_library.fasta"
# input_constructs = "5_parts_constructs.fasta"
# output_constrcuts_fasta = "prova_assembly.fasta"
# output_verification_csv = "prova_verification.csv"
# ########################################################################


"""
Main
"""


def parts_assembly(input_parts: 'input parts fasta filename',
                   input_constructs: 'input constructs filename',
                   output_combinatorial_assembly: 'output combinatorial assembly fasta filename',
                   output_assembly_mapping: 'output verification table csv filename ',
                   parts_number=5,
                   max_parts=10):

    if parts_number <= max_parts:
        # Initialsing parts position dictionary
        position_parts_dict = build_position_parts_dictionary(input_constructs, parts_number)

        # reading parts from input
        parts_list = []
        for record in screed.open(input_parts):
            part = Part()
            part.build_part(record, position_parts_dict)
            parts_list.append(part)

        # Initialising nested lists of all parts combinations
        combinations = get_combinations(list(position_parts_dict.values()))

        # Building constructs
        construct_id = 0
        construct_list = []
        for parts_id in combinations:
            construct = Construct()
            construct.build_construct(parts_id, parts_list, construct_id)
            construct_list.append(construct)
            construct_id += 1

        # Writing output files
        write_constructslibrary_fasta(construct_list, output_combinatorial_assembly)
        write_varification_csv(construct_list, output_assembly_mapping)
    else:
        with open(output_combinatorial_assembly, 'w') as fp:
            pass
        with open(output_assembly_mapping, 'w') as fp:
            pass

def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([parts_assembly
                            ])


if __name__ == '__main__':
    main()
