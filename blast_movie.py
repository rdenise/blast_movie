#! /usr/bin/env python3
# -*- coding: utf-8 -*-

##########################################################################################
##########################################################################################
##
##                                Library
##
##########################################################################################
##########################################################################################

import argparse
from textwrap import dedent
import os

import numpy as np
import pandas as pd
import tqdm

import multiprocessing

from common import network_creation
from common import utils
from common import utils_blast

##########################################################################################
##########################################################################################
##
##                                Main
##
##########################################################################################
##########################################################################################

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=dedent(
        """See the effect of alignement threshold based on annotation color"""
    ),
)

general_option = parser.add_argument_group(title="General input dataset options")
general_option.add_argument(
    "-b",
    "--blastfile",
    metavar="<file>",
    dest="blastFile",
    help="Blast all vs all file",
    required=True,
)
general_option.add_argument(
    "-s",
    "--selected_threshold",
    dest="thresholdName",
    choices=["score", "pident", "coverage", "evalue"],
    help="Choose the threshold you want to see evolve between ['score', 'pident', 'coverage', 'evalue']",
    required=True,
)
general_option.add_argument(
    "-a",
    "--annotation",
    metavar="<annotation>",
    dest="annotation",
    help="Tabulated file with the information about the sequence need to have at least, 'protein_id', 'length' and the columns to group the protein",
    required=True,
)
general_option.add_argument(
    "-c",
    "--column_name",
    metavar="<column_name>",
    dest="columnName",
    help="Name of the column in annotation file that contain the group you want to highligh in the figure",
    required=True,
)
general_option.add_argument(
    "-lcc",
    "--length_choice_cov",
    default="mean",
    dest="length_choice_cov",
    choices=["mean", "subject", "query", "shortest", "longest"],
    help=dedent(
        """
                            Length used for percentage overlap calculation 
                                       between 2 sequences:
                                       'mean'=mean of the 2 lengths (default),
                                       'subject'=subject length, 'query'=query length,
                                       'shortest'=shortest length, 'longest'=longest length
                            """
    ),
)
general_option.add_argument(
    "-id",
    "--length_choice_id",
    default="mean",
    dest="length_choice_id",
    choices=["mean", "subject", "query", "shortest", "longest", "HSP"],
    help=dedent(
        """
                            Length used for percentage identity calculation 
                                       between 2 sequences:
                                       'mean'=mean of the 2 lengths (default),
                                       'subject'=subject length, 'query'=query length,
                                       'shortest'=shortest length, 'longest'=longest length
                                       'HSP'=HSP length
                            """
    ),
)
general_option.add_argument(
    "-o",
    "--output",
    default=None,
    dest="output",
    metavar="<OUTPUT>",
    help="Name of the output file (default: In the same folder as the blast output)",
)
general_option.add_argument(
    "-t",
    "--threads",
    metavar="<num_threads>",
    dest="threads",
    help="Number of threads to use (default:1)",
    required=True,
    default=1,
    type=int,
    choices=range(1, multiprocessing.cpu_count()),
)

##########################################################################################

args = parser.parse_args()

##########################################################################################

if args.output:
    OUTPUT = args.output
else:
    OUTPUT = os.path.join(os.getcwd(), os.path.basename(args.blast_tbl))

##########################################################################################

utils.create_folder(OUTPUT)

BLAST_TBL = args.blastFile
ANNOT = args.annotation

##########################################################################################

annot_dtypes = {
    "protein_id": "string",
    "length": "int",
    args.columnName: "string",
}

annot_df = pd.read_table(
    ANNOT, usecols=["protein_id", "length", args.columnName], dtype=annot_dtypes
).set_index("protein_id")

# Opening blast_out and preparation
blast_names = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


# Get the types of te columns for multiple HSPs dataframe
blast_dtypes = [
    ("qseqid", "S100"),
    ("sseqid", "S100"),
    ("pident", np.float64),
    ("length", np.int32),
    ("mismatch", np.int32),
    ("gapopen", np.int32),
    ("qstart", np.int32),
    ("qend", np.int32),
    ("sstart", np.int32),
    ("send", np.int32),
    ("evalue", np.float64),
    ("bitscore", np.float64),
]

dict_prot = annot_df.length.to_dict()

blast_df = utils_blast.read_blast(
    blast_out=BLAST_TBL,
    blast_dtypes=blast_dtypes,
    protein_dict=dict_prot,
    option_pid=args.length_choice_id,
    option_cov=args.length_choice_cov,
)

all_score = blast_df[args.thresholdName].sort_values().unique()
num_score = all_score.shape[0]
max_str_len_score = len(str(max(all_score)))


G = network_creation.create_graph(blast_tbl=blast_df, thresholdName=args.thresholdName)


print()
print("----------------------")
print("Plotting network files")
print("----------------------")

if args.threads == 1:
    results = []

    for score_index in tqdm.tqdm(range(num_score)):
        score = all_score[score_index]

        score_txt = str(int(score)).zfill(max_str_len_score)

        name_output = os.path.join(
            OUTPUT,
            os.path.basename(BLAST_TBL).replace(
                os.path.splitext(BLAST_TBL)[-1], f".{score_txt}.png"
            ),
        )

        tmp_results = network_creation.visu_graph(
            graph=G,
            output=name_output,
            threshold=score,
            annot_dict=annot_df.loc[
                ~annot_df[args.columnName].isna(), args.columnName
            ].to_dict(),
            thresholdName=args.thresholdName,
        )

        results.append(tmp_results)
else:
    args_func = []
    for score_index in range(num_score):
        score = all_score[score_index]
        score_txt = str(int(score)).zfill(max_str_len_score)

        name_output = os.path.join(
            OUTPUT,
            os.path.basename(BLAST_TBL).replace(
                os.path.splitext(BLAST_TBL)[-1], f".{score_txt}.png"
            ),
        )

        args_func.append(
            (
                G,
                name_output,
                score,
                annot_df.loc[
                    ~annot_df[args.columnName].isna(), args.columnName
                ].to_dict(),
                args.thresholdName,
            )
        )

    pool = multiprocessing.Pool(args.threads)
    results = list(
        tqdm.tqdm(
            pool.imap(network_creation.star_visu_graph, args_func), total=num_score
        )
    )
    pool.close()

df = pd.DataFrame(results, columns=["threshold_value", "red_edges"])

df.to_csv(
    os.path.join(OUTPUT, f"number_of_red_edges_for_{args.thresholdName}.tsv"),
    index=False,
    sep="\t",
)
