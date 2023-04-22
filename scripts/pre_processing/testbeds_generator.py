from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
import random
import argh
import pandas as pd
import json
import copy

"""
Class describing a DNA part 
"""


class DNApart:
    def __init__(self):
        # seq ="ATGATCTAGTAT...TGA"
        self.seq = None
        # number = integer describing the part e.g. 1 stands for part 1
        self.number = None
        # features read from the inpunt gebankfile
        self.features = None
        # length of the DNA part
        self.length = None
        # if different from None describe the template parts used for cloning
        self.clone = None
        # integer representing position in the construct
        self.position = None


"""
Class describing a construct, a set of parts
"""


class Construct:
    def __init__(self):
        # seq ="ATGATCTAGTAT...TGA"
        self.seq = None
        # number = integer describing the  construct e.g. 1 stands for construct 1
        self.number = None
        # features read from the inpunt gebankfile
        self.features = None
        # parts string describing which parts are included into the construct eg. 5 3 2
        self.parts_string = None
        # list of DNA parts objects to be included into the construct
        self.parts = []
        # number of parts included into the construct
        self.length = None
        #similarity
        self.clone = None

    def to_dict(self):
        return {
            'Construct': self.number,
            'Parts': self.parts_string,
            'Clone': self.clone,
        }



def parts_initialiser(genbank_file,cds_number):
    """
       Finding CDS
       and creating a list of all CDS, a number of CDS = cds_number is selected at random a tunred into a list of parts
       1 part = 1 random CDS taken from yeast chromosome XIV
    """
    # wildtype_cds_percentage = int((1 - clone_parts_amount) *100)
    # cds_to_download =  int((wildtype_cds_percentage * cds_number)/ 100)
    genbank_record = SeqIO.read(genbank_file, "genbank")
    cdss= []
    for feature in genbank_record.features:
        if feature.type == "CDS":
            cds = DNApart()
            cds.seq = genbank_record.seq[feature.location.start:feature.location.end]
            cds.length = len(cds.seq)
            cdss.append(cds)

    part_number = 0
    parts = random.sample(cdss, cds_number)
    for part in parts:
        part_number += 1
        part.number = part_number
    return parts


def read_cloned_parts(construct):
    """
    add a string describing witch parts is clone and which is not
    :param construct:
    :return:
    """
    construct.clone = ""
    for part in construct.parts:
        construct.clone = construct.clone + str(part.clone) + ", "
    return construct

def from_parts_to_construct(constructs_number, length_construct, parts):
    """
    Methods to create n constructs from a list of parts object return
    a list of construct objects
    """
    #Creating a construct joining at random length_construct CDS
    position_list = list(range(length_construct))
    parts_positions = create_parts_positions_groups(parts,position_list)

    constructs = []
    for construct_index in range(constructs_number):
        construct = Construct()
        construct.number = construct_index
        construct.length = length_construct
        for position in range(0,length_construct):
            selected_part = random.choice(parts_positions[position])
            construct.parts.append(selected_part)

        # Parts assembly
        construct.parts_string = ""
        construct = assembly(construct)
        construct = read_cloned_parts(construct)
        constructs.append(construct)
    return constructs

# def test_mutagenesis(constructs_number, length_construct, parts, parts_similarity):
#         #Initialising constructs list
#         constructs = []
        

#         #Reset parts clone
#         for part in parts:
#             part.clone = None
        
#         # Create a template construct
#         template_construct = Construct()
#         template_construct.number = "0"
#         template_construct.length = length_construct
#         for position in range(0,length_construct):
#                 selected_part = random.choice(parts)
#                 selected_part.position = position
#                 if selected_part not in template_construct.parts:
#                     template_construct.parts.append(selected_part)

#         # Parts assembly
#         template_construct.parts_string = ""
#         template_construct = assembly(template_construct)
#         template_construct = read_cloned_parts(template_construct)
#         constructs.append(template_construct)

#         #Copying template construct and editing last part
#         for construct_index in range(1,constructs_number):
#             construct = Construct()
#             construct.number = construct_index
#             construct.length = length_construct
#             construct.parts = copy.deepcopy(template_construct.parts)
#             part_clone = clone_part(template_construct.parts[-1],construct.parts[-1])
#             mutation = generate_random_mutation(part_clone, parts_similarity)
#             part_clone.seq = insert_mutation(part_clone, mutation)
#             part_clone.number = part_clone.number + construct_index
#             part_clone.clone =template_construct.parts[-1].number
#             construct.parts[-1] = part_clone
            
#             # Parts assembly
#             construct.parts_string = ""
#             construct = assembly(construct)
#             construct = read_cloned_parts(construct)
#             constructs.append(construct)
#         return constructs
       #Copy from template and clone only last part
        
# def mutagenesis_parts_library(constructs):
#     parts_library = []
#     part_names = []
#     for construct in constructs:
#         for part in construct.parts :
#             if part.number not in part_names:
#                 part_names.append(part.number)
#                 parts_library.append(part)
#     return parts_library 

def assembly(construct):
    """
    method to put parts into construct
    :param construct:
    :return:
    """
    # Parts assembly
    construct_parts_seq = []
    for part in construct.parts:
        construct_parts_seq.append(part.seq)
        construct.parts_string = construct.parts_string + str(part.number) + ", "
    construct.seq = sum(construct_parts_seq, Seq("", DNAAlphabet()))
    return construct


def clone_part(part,part_0):
    """
    method to create a clone of a part
    :param construct:
    :return:
    """
    part_clone = DNApart()
    part_clone.seq = part_0.seq
    part_clone.features = part_0.features
    part_clone.length = part_0.length
    part_clone.clone = part_0.number
    part_clone.position = part_0.position
    part_clone.number = part.number
    return part_clone



def generate_random_mutation(part_clone, parts_similarity):
    """
    Method to create a random sequence "mutation"
    to differentiate sequences
    :param clone:
    :param similarity:
    :return: mutation
    """
    mutation_length = ((1-parts_similarity) * 100) * len(part_clone.seq)/100

    mutation_length = int(mutation_length)
    mutation = ""
    for nucleotide in range(mutation_length):
        random_nucleotide = random.choice(["A","T","C","G"])
        mutation += random_nucleotide
    mutation = Seq(mutation, DNAAlphabet())
    return mutation

def insert_mutation(part_clone, mutation):
    """
    method to add the mutation to the construct
    :param clone:
    :param mutation:
    :return:
    """
    mutation_left_limit = len(part_clone.seq) - len(mutation)
    insertion_site = random.randrange(0,mutation_left_limit)
    part_clone.seq = part_clone.seq[:insertion_site] + mutation + part_clone.seq[insertion_site+len(mutation):]
    return part_clone.seq


def postions_to_mutate(construct_similarity,length_construct,method):
    """
    Method to determine position to mutate
    :param construct_similarity: % of cloned parts per construtcs
    :param length_construct: parts per constrcutd
    :param method: first, random or last, select the first , random or lost parts of constructs
    :return: list of parts postions to clone
    """
    construct_percentage = construct_similarity * 100
    parts_per_construct_to_mutate =int((construct_percentage * length_construct)/100)

    parts_indexes_to_clone = None
    if method == "first":
        parts_indexes_to_clone =list(range(0,parts_per_construct_to_mutate))
    elif method == "random":
        positions = list(range(length_construct))
        parts_indexes_to_clone = list(random.sample(positions,parts_per_construct_to_mutate))
    elif method == "last":
        parts_indexes_to_clone = list(range(length_construct-parts_per_construct_to_mutate,length_construct ))
    return parts_indexes_to_clone


def parts_cloning(parts,parts_indexes_to_clone,parts_similarity):
    """
    Takes a part_0 sample for each position where a clone is needed e.g. pos 2,3,4,
    all the other parts in the same position will be clone of parts 0
    :param parts: list of parts objects
    :param parts_indexes_to_clone: list of postion indexes , state which position mutate
    :param parts_similarity: similarity between parts
    :return: list of parts
    """
    if parts_similarity > 0:
        cloning_groups = create_parts_positions_groups(parts,parts_indexes_to_clone)
        clones = []
        for cloning_group in cloning_groups:
            part_0_index = random.randint(0,len(cloning_group)-1)
            part_0 = cloning_group[part_0_index]
            clones.append(part_0)
            for part in cloning_group:
                if cloning_group.index(part) != part_0_index:
                    part_clone = clone_part(part,part_0)
                    mutation = generate_random_mutation(part_clone, parts_similarity)
                    part_clone.seq = insert_mutation(part_clone, mutation)
                    clones.append(part_clone)

        # adding parts in non cloning positions
        for part in parts:
            if part.position not in parts_indexes_to_clone:
                clones.append(part)
        return clones
    else:
        return parts

def create_parts_positions_groups(parts,parts_indexes_to_clone):
    """
    split the parts list into a nested list including a list for each part to clone
    e.g. indexes = 2,3,4 return = [[all_parts_2], [all_parts_3], [all_parts_4]]
    :param parts: list of parts objects
    :param parts_indexes_to_clone: list of parts position to clone
    :return:  nested list including parts object to clone diveded by position
    """
    cloning_groups = []
    for position in range(max(parts_indexes_to_clone) + 1):
        empty_list = []
        cloning_groups.append(empty_list)

    for part in parts:
        if part.position in parts_indexes_to_clone:
            cloning_groups[part.position].append(part)
    cloning_groups = [x for x in cloning_groups if x != []]
    return cloning_groups

def assign_parts_position(parts,length_construct):
    """
    assign to each part which position must assume in a construct, and which
    must be cloned
    :param parts: list of parts objects
    :param parts_similarity:  % of parts with similarity in construct
    :param position_similarity: str choice first, last, random
    :return:
    """
   # splitting the part list into sublists of length length construct
    parts_groups = [parts[x:x + length_construct] for x in range(0, len(parts), length_construct)]
    # for each list recursively assign the index
    for part_group in parts_groups:
        part_index = 0
        for part in part_group:
            part.position = part_index
            part_index += 1
            if part_index == length_construct :
                part_index = 0

    #If there is spare sublist assign position at random
    if len(parts_groups[-1]) < length_construct:
        parts_indexes = list(range(0,length_construct))
        for part in parts_groups[-1]:
            part.position = random.choice(parts_indexes)


def write_constructslibrary_fasta(constructs, filename):
    """
    Turns a list of constructs into a fasta file with all constructs
    """
    # writing a fasta file for constructs
    with open(filename, "w") as fastaseq:
        for construct in constructs:
            seq_record = SeqRecord(
                construct.seq,
                id=str(construct.number),
                name="",
                description=">" + "Construct " + str(construct.number) + " " + "Parts " + str(
                    construct.parts_string) + "\n",
                )
            SeqIO.write(seq_record, fastaseq, "fasta")

def write_partslibrary_fasta(parts_library, filename):
    """
    turn a list of parts into a fasta file with all parts
    :return:
    """
    # writing a fasta file for parts
    with open(filename, "w") as fastaseq:
        seq_record = None
        for part in parts_library:
            if part.clone == None:
                seq_record = SeqRecord(
                    part.seq,
                    id=str(part.number),
                    name="",
                    description=">" + "Part " + str(part.number) + "\n",
                )
            elif part.clone != None:
                seq_record = SeqRecord(
                    part.seq,
                    id=str(part.number),
                    name="",
                    description=">" + "Part " + "Clone " + str(part.clone) + "\n",

                )
            SeqIO.write(seq_record, fastaseq, "fasta")


def write_csv(constructs, output_verifcation_csv):
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

def write_json_clone_parts(parts_library, json_filename):
    """
    Given a list of parts write as output a json file clone parts : 51,34,56,64
    :param parts:
    :param filename:
    :return:
    """
    json_dictionary = {}

    clone_id_parts = []
    for part in parts_library:
        if part.clone != None:
            clone_id_parts.append(part.number)

    json_dictionary["cloned parts id"] = clone_id_parts
    with open(json_filename, 'w') as outfile:
        json.dump(json_dictionary, outfile)



"""
main
"""


def generate_testbeds(input_genbank: 'input genbank file',
                      output_constructs_fasta: 'output constructs',
                      output_parts_fasta: 'output parts',
                      output_verifcation_csv: 'output csv',
                      output_cloned_parts_json: 'output cloned parts json',
                      # Number of CDS to download from the input DNA sequence
                      cds_number: 'number of cds' = 50,
                      # constructs library size
                      constructs_number: 'number of contructs' = 100,
                      # number of parts in each construct
                      length_construct: 'number of CDS in contructs' = 20,
                      constructs_similarity=0.4,
                      parts_similarity=0.7,
                      mutation_method = 'last'):

    """
    Takes as input -input-genbank file and reads all the CDS.
    The first - cds-number CDS are taken.
    Each of the CDS is turned into a part object, - length-constructs parts are randomly assembled into
    -constructs-number constructs
    Each construct is written into the fasta file - output_constructs_fasta
    Each part is written into the fasta file - output_parts_fasta


    :param input_genbank: genbank file to use as input
    :param output_constructs_fasta: filename for the output fasta file containing all constructs
    :param output_parts_fasta: filename for the output fasta file containing all parts
    :param output_parts_fasta: filename for the .csv file containing constructs and parts
    :param cds_number: pool of CDS to take from the input genbank file
    :param constructs_number: number of constructs to generate
    :param length_construct: number of parts to include in each construct
    :return: a fasta file including constructs, a fasta file including the library of parts
    a .csv file mapping constructs and parts
    """
    # Generating parts
    parts = parts_initialiser(input_genbank,cds_number)
    # Assigning positions to each part
    assign_parts_position(parts, length_construct)
    # Selecting positions to clone
    cloned_positions = postions_to_mutate(constructs_similarity, length_construct,mutation_method)
    # Adding clones
    parts_library = parts_cloning( parts,cloned_positions, parts_similarity)
    # Building Constructs
    constructs = from_parts_to_construct(constructs_number, length_construct, parts_library)
    # Writing outputs
    write_constructslibrary_fasta(constructs, output_constructs_fasta)
    write_partslibrary_fasta(parts_library, output_parts_fasta)
    write_csv(constructs, output_verifcation_csv)
    write_json_clone_parts(parts_library,output_cloned_parts_json)


def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([generate_testbeds
                            ])


if __name__ == '__main__':
    main()
