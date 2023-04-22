import pysam
import re
import pandas as pd
import argh
import os
import mappy as mp
import logging
import sys


class Alignment:
    """
    Class to analyse minimap2 results
    """

    def __init__(self, sam_dictionary):
        # dictionary uploaded from a single sam file line including all features
        self.sam_dictionary = sam_dictionary
        # id related to read query
        self.read_name = None
        # construct name reference
        self.construct_name = None
        # score related to the alignment
        self.alignment_score = None
        # mapping quality
        self.mapping_quality = None

    def to_dict(self):
        return {
            'read id': self.read_name,
            'construct name': self.construct_name,
            'alignment score': self.alignment_score,
            'mapping quality': self.mapping_quality
        }

    # fill the alignment object assigning read id, construct id and alignment score
    def build_alignment(self):
        self.read_name = self.sam_dictionary["name"]
        self.construct_name = self.sam_dictionary["ref_name"]
        self.mapping_quality = self.sam_dictionary["map_quality"]
        if len(self.sam_dictionary["tags"]) > 2 :
            alignment_score = self.sam_dictionary["tags"][2]
            alignment_score = int(re.sub("[^0-9]", "", alignment_score))
            self.alignment_score = alignment_score
        else:
            self.alignment_score = None

    # this method is called if a new alignment for the same read has a greater score
    # updates construct id and alignment score
    def update_alignment(self, alignment_score, sam_dictionary):
        self.alignment_score = alignment_score
        self.construct_name = sam_dictionary["ref_name"]


def minimap2_to_csv(sam_file_filename: 'input sam file',
                    output_filename: 'output csv table'):

    # Reading input and initialising keys variables
    try:
        samfile = pysam.AlignmentFile(sam_file_filename, "rb")
        reads_scores = {}
        alignments = []

        # For each line in the sam file look at the read
        for alignment in samfile:
            sam_dictionary = alignment.to_dict()
            read = sam_dictionary["name"]

            # if is the first time I see this read and there is an alignment score build an alignment object
            if read not in reads_scores.keys() :
                aln = Alignment(sam_dictionary)
                aln.build_alignment()
                alignments.append(aln)
                reads_scores[read] = aln.alignment_score

            # otherwise for the same read take the alignment with the highest alignment score
            elif read in reads_scores.keys() :
                alignment_score = sam_dictionary["tags"][2]
                alignment_score = int(re.sub("[^0-9]", "", alignment_score))
                if reads_scores[read] < alignment_score:
                    reads_scores[read] = alignment_score
                    for aln in alignments:
                        if aln.read_name == read:
                            aln.update_alignment(alignment_score,sam_dictionary)


       # Build output table
        minimap_2_result = pd.DataFrame.from_records([s.to_dict() for s in alignments])
        minimap_2_result.to_csv(output_filename)

    except ValueError:
        with open(output_filename, 'w') as fp:
            pass


def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([minimap2_to_csv
                            ])


if __name__ == '__main__':
    main()
