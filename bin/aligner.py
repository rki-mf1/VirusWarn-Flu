#!/usr/bin/env python

"""
Historical pairwise global aligner using the needleman_wunsch algorithm and
helper functions for aligning sequences from the BioPython pairwise2 library

@StefanFrankBio @hrichard
"""

from typing import Tuple
import parasail


"""
Parasail part
Slightly modified version taken from dev/covsonar2-alpha covsonar/src/covsonar/align.py

@stephan-fuchs
"""

def parasail_align(
    query_seq: str, 
    ref_seq: str, 
    gapopen: int = 16, 
    gapextend: int = 4, 
) -> Tuple[str, str]:
    """
    Use semi-global Smith-Waterman algorithm from parasail-python to align the sequence to the reference

    Args:
        query_seq (str): Query sequence to align.
        ref_seq (str): Reference sequence to align against.
        gapopen (int, optional): Gap open penalty. Defaults to 16.
        gapextend (int, optional): Gap extension penalty. Defaults to 4.

    Returns:
        Tuple[str, str]: Tuple containing the aligned reference sequence and query sequence
    """
    matrix = parasail.dnafull
    alignment = parasail.sg_trace_striped_32(
            query_seq, ref_seq, gapopen, gapextend, matrix
    )

    aligned_ref = alignment.traceback.ref
    aligned_query = alignment.traceback.query

    return aligned_ref, aligned_query