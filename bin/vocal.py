#!/usr/bin/env python

"""
Vocal is a simple Variants of Concern Alert System Script

@hrichard @StefanFrankBio @hoelzer
"""

import argparse
import time
from typing import Tuple

from Bio import pairwise2
from Bio import SeqIO
import pandas as pd

import aligner
import AnnotateVariants
from data_loader import D_GENEPOS
from data_loader import GAP_CHAR
import PSL_helper

# Temp constants
ref_id = "HA"
ref_start, ref_end = D_GENEPOS[ref_id][1:]


def Pairwise_align_query_ref(
    ref_seq: str,
    query_seq: str,
    align_method: str = "parasail",
) -> Tuple[str, str]:
    """
    Performs the alignment of query on reference using a pairwise alignment.

    Args:
        query_seq (str): Query sequence to align.
        ref_seq (str): Reference sequence to align against.
        align_method (str, optional): Aligner to use. Defaults to parasail.

    Returns:
        Tuple[str, str]: Tuple containing the aligned reference sequence and query sequence
    """
    aligner_func = None
    if align_method == "parasail":
        print("== Will be using parasail-python alignment Module")
        aligner_func = aligner.parasail_align
    else:
        raise ValueError("Unknown aligner function!")

    print(
        f"Finding Alignment of the {ref_id} protein (length: {len(ref_seq)}) in the consensus"
    )
    
    s_time = time.time()
    ref_al, query_al = aligner_func(ref_seq, query_seq)
    print("{:.1f} seconds needed".format(time.time() - s_time))
    print(
        f"Found alignment of length {len(ref_al)}"
    )
    return ref_al, query_al


def PSL_align_query_ref(
    ref_seq,
    query_record,
    PSL_df,
    gene_nt_seq,
    gene_nt_start=ref_start,
    gene_nt_end=ref_end,
    k_size=18,
):
    """
    Returns the alignment of the reference to the query for the gene gene_nt_seq using the
    information in the PSL alignment dataframe
    the parameter k_size is the size of the seed used to find the position of the
    protein in the aligned reference sequence. So in theory gene_nt_seq occurs in ref_seq
    (but this is not tested for)
    """
    try:
        PSL_row = PSL_df.loc[query_record.id]
    except KeyError:
        raise KeyError()

    (
        ref_al_complete,
        query_al_complete,
        r_coord_corresp,
        q_coord_corresp,
    ) = PSL_helper.PSLpretty(ref_seq, str(query_record.seq), PSL_row)

    # We need a converter of positions from the PSLpretty function
    # How many nucleotides to find the beginning / end of the spike?
    ankor_start = gene_nt_seq[:k_size]
    ankor_end = gene_nt_seq[-k_size:]
    # Find the beginning of the spike on the ref/query alignment now.
    # gene_start_ankor = ref_al_complete.find(ankor_start)
    # gene_end_ankor = ref_al_complete.find(ankor_end) + k_size #exclusive
    gene_start_corresp = PSL_helper.posconverter(
        gene_nt_start - 1, r_coord_corresp
    )  # position is 1 based
    gene_end_corresp = PSL_helper.posconverter(
        gene_nt_end, r_coord_corresp
    )  # position is exclusive
    # if gene_start_ankor == -1 and gene_end_ankor == -1:
    #     print(f"Start (anchor, corresp): ({gene_start_ankor}, {gene_start_corresp})")
    #     print(f"End (anchor, corresp): ({gene_end_ankor}, {gene_end_corresp})")
    #     warnings.warn( f"The seed of size {k_size} could not be found, sequence: {query_record.id}")
    ref_al = ref_al_complete[gene_start_corresp:gene_end_corresp]
    query_al = query_al_complete[gene_start_corresp:gene_end_corresp]
    offset = (
        len(query_al_complete[:gene_start_corresp].replace(GAP_CHAR, ""))
        + PSL_row["tStart"]
        + 1
    )
    return (offset, ref_al, query_al)


def main():
    """
    Main function

    Returns
    ---------
    None
    """
    parser = argparse.ArgumentParser(
        description="Vocal: Alert system for variants of concerns on the Spike protein",
        epilog="Usage: python vocal.py -i sequences.fasta [-o variant_table.tsv]",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Set of complete consensus genome sequences (fasta format)",
    )
    parser.add_argument(
        "-r",
        "--ref_nt",
        required=True,
        help="Reference sequence (DNA version)",
    )
    parser.add_argument(
        "--control",
        required=True,
        help="Sequence with MOC and ROI",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="variant_table.tsv",
        help="output table with the set of variants detected over all sequences, with translation",
    )
    parser.add_argument(
        "--PSL",
        type=str,
        default="",
        help="PSL file with the alignment of the sequences (deactivates pairwise alignment)",
    )
    parser.add_argument(
        "--president",
        action="store_true",
        help="indicates that the PSL file is given as a president report (default: False)",
    )
    parser.add_argument(
        "--complete",
        required=True,
        help="Only consider complete genomes",
    )
    parser.add_argument(
        "-n",
        "--number",
        required=True,
        help="Number of nucleotides a sequence can differ from the length of the reference sequence to be considered as complete",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="verbose mode (prints alignments on stdout)",
    )
    args = parser.parse_args()

    # Setting up the alignment strategy
    al_module = "parasail"

    PSL_parse = False
    if len(args.PSL) > 0:
        PSL_parse = True
        if args.president:
            PSL_df = PSL_helper.readPSLfromPresident(args.PSL)
        else:
            PSL_df = PSL_helper.readPSL(args.PSL)
    
    ref_nt = [str(rec.seq) for rec in SeqIO.parse(open(args.ref_nt),'fasta')][0]

    nrecords = 0
    nskipped = 0
    nincomplete = 0
    list_of_dfs = []
    fastain = args.input
    controlin = args.control
    verbose = args.verbose

    n = int(args.number)
    if args.complete == 'yes':
        complete = True
        print(
            "FluWarnSystem running only for COMPLETE SEQUENCES"
        )
        print(
            f"The reference sequence is {len(ref_nt)} nucleotides long"
        )
        print(
            f"A sequence is considered as complete with a length between {len(ref_nt)-n} and {len(ref_nt)+n}"
        )
    elif args.complete == 'no':
        complete = False
        print(
            "FluWarnSystem running for ALL SEQUENCES"
        )
    else:
        raise ValueError("The parameter complete only takes the options yes and no!")

    # Add a sequence with a MOC and ROI to catch the error from prepare_report.R 
    # that occurs when there is no MOC / ROI in input fasta
    # The sequence is removed in report.Rmd and is not in the final report
    with open(fastain, 'a') as in_file:
        with open(controlin, 'r') as control_file:
            for line in control_file:
                in_file.write(line)

    print(f"Opening {fastain}")
    for record in SeqIO.parse(fastain, "fasta"):
        print("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
        print(record.description)

        ##Quick hack to have the record IDs with the spaces masked (compatible with president)
        record.id = record.description.replace(" ", "%space%")
        query_seq = str(record.seq.upper())

        if complete == True and record.description == "OFL3839 | HA | A/Germany/767/95 | OFL_ISL_789 | AF008754 | H3N2":
            print(
                "MOC and ROI sequence is included, even though it is incomplete"
            )
            print(
                "It is removed before the report is generated"
            )
        elif complete == True and (len(record.seq) < (len(ref_nt)-n) or len(record.seq) > (len(ref_nt)+n)):
            SeqIO.write(record, 'incomplete_seq.fasta', 'fasta')
            print(
                f"Found incomplete sequence of length {len(record.seq)}, the sequence will not be processed"
            )
            nincomplete += 1
            continue
        else:
            print(
                f"The sequence is {len(record.seq)} nucleotides long"
            )

        if not set(query_seq) <= {"A", "T", "G", "C", "N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}:
            print(record.id)
            print(
                "Found non IUPAC nucleotide characters in query, the sequence will not be processed"
            )
            print(
                f"Characters found {set(query_seq) - set(['A', 'T', 'G', 'C', 'N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'])}"
            )
            nskipped += 1
            continue

        if 'R' in query_seq: # IUPAC A or G
            query_seq = query_seq.replace('R', 'N')
        if 'Y' in query_seq: # IUPAC C or T
            query_seq = query_seq.replace('Y', 'N')
        if 'S' in query_seq: # IUPAC G or C
            query_seq = query_seq.replace('S', 'N')
        if 'W' in query_seq: # IUPAC A or T
            query_seq = query_seq.replace('W', 'N')
        if 'K' in query_seq: # IUPAC G or T
            query_seq = query_seq.replace('K', 'N')
        if 'M' in query_seq: # IUPAC A or C
            query_seq = query_seq.replace('M', 'N')
        if 'B' in query_seq: # IUPAC G or C or T
            query_seq = query_seq.replace('B', 'N')
        if 'D' in query_seq: # IUPAC A or G or T
            query_seq = query_seq.replace('D', 'N')
        if 'H' in query_seq: # IUPAC A or C or T
            query_seq = query_seq.replace('H', 'N')
        if 'V' in query_seq: # IUPAC A or G or C
            query_seq = query_seq.replace('V', 'N')

        if nrecords % 100 == 0:
            print(f"Sequence number {nrecords+1}: {record.id}")
        try:
            if PSL_parse:
                offset, ref_al, query_al = PSL_align_query_ref(
                    ref_nt, 
                    query_seq, 
                    PSL_df,
                )
            else:
                ref_al, query_al = Pairwise_align_query_ref(
                    ref_nt,
                    query_seq,
                    al_module,
                )
        except KeyError:
            print("This key doesn't exist>", record.id)
            nskipped += 1
            continue

        if verbose:
            print("Alignment of query to ref:")
            print(pairwise2.format_alignment(ref_al, query_al, 0, 1, len(ref_al)))
        
        # We only process sequences with ACGT- characters
        if not set(query_al) <= {"A", "T", "G", "C", "-", "N"}:
            print(record.id)
            print(
                "Found non ACGTN- characters in alignment of query, the sequence will not be processed"
            )
            print(
                f"Characters found {set(query_al) - set(['A', 'T', 'G', 'C', '-', 'N'])}"
            )
            nskipped += 1
            continue
        nrecords += 1

        df = AnnotateVariants.alignedCDSvariants(ref_al, query_al, verbose=verbose)
        
        if not df.empty:
            df["ID"] = record.id
            df["target_gene"] = ref_id
            list_of_dfs.append(df)

    print(f"Processed {nrecords} sequences, skipped {nskipped}, incomplete {nincomplete}")
    df_all_mutations = pd.concat(list_of_dfs, ignore_index=True, sort=False)

    df_all_mutations.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    main()
