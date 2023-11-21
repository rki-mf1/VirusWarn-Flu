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
# ref_start = 21563
# ref_end   = 25384
ref_start, ref_end = D_GENEPOS[ref_id][1:]
### We fixed an divergence of max 5% with the reference
### This is a very conservative estimate
#max_percent_shift = 0.05
### Influenza genome is about 13600 bases long
#max_shift_allowed = int(13600 * max_percent_shift)
QUERY_START_RESTRICTED = 1 #ref_start - max_shift_allowed - 1
QUERY_END_RESTRICTED = 1750 #ref_end + max_shift_allowed


def Pairwise_align_query_ref(
    ref_seq: str,
    query_seq: str,
    align_method: str = "parasail",
    align_start_on_ref: int = QUERY_START_RESTRICTED,
    align_end_on_ref: int = QUERY_END_RESTRICTED,
) -> Tuple[str, str]:
    """
    Performs the alignment of query on reference using a pairwise alignment.

    Args:
        query_seq (str): Query sequence to align.
        ref_seq (str): Reference sequence to align against.
        align_method (str, optional): Aligner to use. Defaults to parasail.
        align_start_on_ref (int, optional): Start of the query to align. Defaults to QUERY_START_RESTRICTED.
        align_end_on_ref (int, optional): End of the query to align. Defaults to QUERY_END_RESTRICTED.

    Returns:
        Tuple[str, str]: Tuple containing the aligned reference sequence and query sequence
    """
    aligner_func = None
    if align_method == "parasail":
        print("== Will be using parasail-python alignment Module")
        aligner_func = aligner.parasail_align
    else:
        raise ValueError("Aligner function option not known")

    print(
        f"Finding Alignment of the {ref_id} protein (length: {len(ref_seq)}) to the subsequence ({align_start_on_ref}-{align_end_on_ref}) in the consensus"
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
        "--ref_nt",
        required=True,
        help="Reference sequence (DNA version)",
    )
    parser.add_argument(
        "--ref_aa",
        required=True,
        help="Reference sequence (Protein version)",
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
        "--restrict-start",
        type=int,
        default=QUERY_START_RESTRICTED,
        help="start coordinate for restricted alignment of the query sequences (estimated position of the spike)",
    )
    parser.add_argument(
        "--restrict-end",
        type=int,
        default=QUERY_END_RESTRICTED,
        help="END coordinate for restricted alignment of the query sequences (estimated position of the spike)",
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
    
    # 1 is for HA, 0 is NA
    ref_nt = [str(rec.seq) for rec in SeqIO.parse(open(args.ref_nt),'fasta')][1]
    ref_aa = [str(rec.seq) for rec in SeqIO.parse(open(args.ref_aa),'fasta')][1]

    nrecords = 0
    nskipped = 0
    list_of_dfs = []
    fastain = args.input
    verbose = args.verbose

    print(f"Opening {fastain}")
    for record in SeqIO.parse(fastain, "fasta"):
        ##Quick hack to have the record IDs with the spaces masked (compatible with president)
        record.id = record.description.replace(" ", "%space%")
        query_seq = str(record.seq.upper())

        if 'R' in query_seq: # IUPAC A or G
            query_seq = query_seq.replace('R', 'N')
        elif 'S' in query_seq: # IUPAC G or C
            query_seq = query_seq.replace('S', 'N')
        
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
                    args.restrict_start,
                    args.restrict_end,
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
        # print(f"Searching for Amino acids differences, found a total of {len(df)} variations")
        
        if not df.empty:
            df["ID"] = record.id
            df["target_gene"] = ref_id
            #df["target_gene_start"] = offset
            # merge with the other dataframes
            # print(df)
            list_of_dfs.append(df)
        # if nrecords >= 1000:
        #    break
        # clean up the output dataframe and order the columns

    print(f"Processed {nrecords} sequences, skipped {nskipped}")
    df_all_mutations = pd.concat(list_of_dfs, ignore_index=True, sort=False)

    df_all_mutations.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    main()
