#!/usr/bin/env python

# Ryan Groussman
# Armbrust Lab 2016

"""
make_seq_info_all.py
Takes as input a trimmed FASTA file with taxid included in a special field in the defline:
Ex:
>Thaps3_262753_tax296543
>Heterosigma_akashiwo_0194584208_MMETSP0292_tax39354
# "seq_info_all.csv must contain at least two columns,
# one column containing the names of the fasta sequences in your tree (seqname) and one with the corresponding tax ids (tax_id)"
"""

import argparse
import re


def write_line_to_csv(line):
    """Given a FASTA defline, outputs the original sequence description
    along with the tax id, separated by a comma.
    If it doesn't find a taxid, writes out an error message for the seq."""

    taxid_found = False
    line_elts = line.split("_")
    for elt in line_elts:
        if tax_pattern.match(elt):
            taxid_found = True
            tax_id = tax_pattern.findall(elt)[0][3:]
            outline = line[1:] + "," + tax_id + "\n"
            output_file.write(outline)

    if taxid_found == False:
        print("Tax ID not found for seq:", line)
        outline = line[1:] + "," + "\n"
        output_file.write(outline)

# parse incoming argument
parser = argparse.ArgumentParser()
parser.add_argument("fasta_file", help="FASTA file with tax_id in defline")
parser.add_argument("-o", "--output", help="Output file name (Default: seq_info_all.csv)", type=str)
args = parser.parse_args()

# here's our regular expression for finding taxids:
tax_pattern = re.compile("(tax[0-9]+)")

# open up the fasta
fasta_file = open(args.fasta_file, 'r')
if args.output == None:
    print("Output file: seq_info_all.csv")
    output_handle = "seq_info_all.csv"
else:
    print("Output file:", args.output)
    output_handle = args.output

output_file = open(output_handle, 'w')

# first, write out a standard header:
out_header='seqname,tax_id' + "\n"
output_file.write(out_header)

# for each sequence in the fasta file;
for line in fasta_file:
    if line.startswith(">"):
        write_line_to_csv(line.strip())


fasta_file.close()
output_file.close()
