#!/usr/bin/env python

"""
Splitting the fasta file from FluPipe into seperate files by subtypes and segments
"""

import argparse
from typing import Tuple

from Bio import SeqIO

def get_substrings(
    mode: str,
    seq_id: str,
) -> Tuple[str, str]:
    """
    Performs the alignment of query on reference using a pairwise alignment.

    Args:
        mode (str): Specify source of data for right splitting of ID
        seq_id (str): ID for the Influenza sequence

    Returns:
        Tuple[str, str]: Tuple containing string with subtype, segment
    """
    if mode == 'FluPipe' or mode == 'flupipe':
        split_id = seq_id.split("|")
        subtype = split_id[2].split("/")[1].strip()
        segment = split_id[3].strip()
    elif mode == 'GISAID' or mode == 'gisaid':
        split_id = seq_id.split("|")
        subtype = split_id[-1].split("/")[1].strip()
        segment = split_id[1].strip()
    elif mode == 'OpenFlu' or mode == 'openflu':
        split_id = seq_id.split("|")
        subtype = split_id[-1].strip()
        segment = split_id[1].strip()

    return subtype, segment
    

def main():
    """
    Main function

    Returns
    ---------
    None
    """
    parser = argparse.ArgumentParser(
        description="Splitting: Into subtypes and segments",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Set of Influenza genome sequences (fasta format)",
    )
    parser.add_argument(
        "-m",
        "--mode",
        required=True,
        help="The source of the sequences: FluPipe, GISAID or OpenFlu",
    )
    parser.add_argument(
        "--h1n1",
        default="h1n1_ha.fasta",
        help="Output fasta file with H1N1 HA sequences",
    )
    parser.add_argument(
        "--h3n2",
        default="h3n2_ha.fasta",
        help="Output fasta file with H3N2 HA sequences",
    )
    parser.add_argument(
        "--vic",
        default="vic_ha.fasta",
        help="Output fasta file with Victoria HA sequences",
    )
    args = parser.parse_args()

    fastain = args.input
    mode = args.mode
    h1n1 = args.h1n1
    h3n2 = args.h3n2
    vic = args.vic

    for record in SeqIO.parse(fastain, "fasta"):
        subtype, segment = get_substrings(mode, record.description)

        if subtype == 'H1N1':
            if segment == 'HA':
                with open(h1n1, 'a') as output:
                    SeqIO.write(record, output, 'fasta')
            #elif segment == 'NA':
            else:
                print(
                    f"Found H1N1 sequence of segment {segment}, the sequence will not be processed"
                )
                SeqIO.write(record, 'h1n1_other_seg.fasta', 'fasta')
        elif subtype == 'H3N2':
            if segment == 'HA':
                with open(h3n2, 'a') as output:
                    SeqIO.write(record, output, 'fasta')
            #elif segment == 'NA':
            else:
                print(
                    f"Found H3N2 sequence of segment {segment}, the sequence will not be processed"
                )
                SeqIO.write(record, 'h3n2_other_seg.fasta', 'fasta')
        elif subtype == 'H0N0' or subtype == "HN":
            if segment == 'HA':
                with open(vic, 'a') as output:
                    SeqIO.write(record, output, 'fasta')
            #elif segment == 'NA':
            else:
                print(
                    f"Found Victoria sequence of segment {segment}, the sequence will not be processed"
                )
                SeqIO.write(record, 'vic_other_seg.fasta', 'fasta')
        else:
            print(
                f"Found sequence of unknown subtype {subtype}, the sequence will not be processed"
            )
            SeqIO.write(record, 'other_subtype.fasta', 'fasta')

if __name__ == "__main__":
    main()
