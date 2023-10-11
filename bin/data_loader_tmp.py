#!/usr/bin/env python

"""
the data to load, temporary (to merge with data_loader from the other script)
"""

import pandas as pd
from Bio import SeqIO

aminoacids = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]

GAP_CHAR = "-"
AASTOP = "*"
codons_stop = ["TAG", "TGA", "TAA"]

CODONDICT = dict(
    A=["GCT", "GCC", "GCA", "GCG"],
    C=["TGT", "TGC"],
    D=["GAT", "GAC"],
    E=["GAA", "GAG"],
    F=["TTT", "TTC"],
    G=["GGT", "GGC", "GGA", "GGG"],
    H=["CAT", "CAC"],
    I=["ATT", "ATC", "ATA"],
    K=["AAA", "AAG"],
    L=["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    M=["ATG"],
    N=["AAT", "AAC"],
    P=["CCT", "CCC", "CCA", "CCG"],
    Q=["CAA", "CAG"],
    R=["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    S=["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    T=["ACT", "ACC", "ACA", "ACG"],
    V=["GTT", "GTC", "GTA", "GTG"],
    W=["TGG"],
    Y=["TAT", "TAC"],
    aastop=codons_stop,
)

AADICT = {}
for k, v in CODONDICT.items():
    for codon in v:
        AADICT[codon] = k

D_GENEPOS = {
    #"PB2": ("A/Brisbane/2/2018", 1, 759),   # Segment 1 - RNA polymerase subunit
    #"PB1": ("A/Brisbane/2/2018", 1, 757),   # Segment 2 - RNA polymerase subunit & PB1-F2 protein
    #"PA": ("A/Brisbane/2/2018", 1, 716),    # Segment 3 - RNA polymerase subunit & PA-X protein
    #"PA-X": ("A/Brisbane/2/2018", 1, 232),  # Segment 3 - RNA polymerase subunit & PA-X protein
    "HA": ("A/Brisbane/2/2018", 1, 566),    # Segment 4 - hemagglutinin
    #"NP": ("A/Brisbane/2/2018", 1, 498),    # Segment 5 - nucleoprotein
    #"NA": ("A/Brisbane/2/2018", 1, 469),    # Segment 6 - neuraminidase
    #"MS1": ("A/Brisbane/2/2018", 1, 252),   # Segment 7 - two matrix proteins
    #"MS2": ("A/Brisbane/2/2018", 1, 97),    # Segment 7 - two matrix proteins
    #"NS1": ("A/Brisbane/2/2018", 1, 219),   # Segment 8 - two non-structural proteins
    #"NS2": ("A/Brisbane/2/2018", 1, 121),   # Segment 8 - two non-structural proteins
}


L_GENES_ROI = [
    ### ------ Glycosylation sites ------ ###
    ### -------- Paper Goka 2014 -------- ###
    (
        "HA",
        213,
        215,
        71,
        71,
        "RegionOfInterest",
        "Glycosylation site",
    ),
    (
        "HA",
        312,
        314,
        104,
        104,
        "RegionOfInterest",
        "Glycosylation site",
    ),
    (
        "HA",
        426,
        428,
        142,
        142,
        "RegionOfInterest",
        "Glycosylation site",
    ),
    (
        "HA",
        432,
        434,
        144,
        144,
        "RegionOfInterest",
        "Glycosylation site (variant)",
    ),
    (
        "HA",
        531,
        533,
        177,
        177,
        "RegionOfInterest",
        "Glycosylation site",
    ),
    (
        "HA",
        516,
        518,
        172,
        172,
        "RegionOfInterest",
        "Glycosylation site (variant)",
    ),
    (
        "HA",
        537,
        539,
        179,
        179,
        "RegionOfInterest",
        "Glycosylation site (variant)",
    ),
    (
        "HA",
        858,
        860,
        286,
        286,
        "RegionOfInterest",
        "Glycosylation site",
    ),
]


def extend_tuple(t):
    name, nt_start, nt_end, aa_start, aa_end, type, fun = t
    return [(name, i, type, fun) for i in range(aa_start, aa_end)]


L_GENES_ROI_EXTENDED = [x for t in L_GENES_ROI for x in extend_tuple(t)]

##Would be easier to expand the table for the following merges
DF_SPIKE_ROI = pd.DataFrame(
    L_GENES_ROI_EXTENDED,
    columns=["gene", "aa_position", "type", "functional domain"],
)

### Load HA from reference file
HA_PROTEIN = [
    str(rec.seq) for rec in SeqIO.parse(open('test-data/A(H1N1)pdm09/A-Brisbane-2-2018_protein.fa'),'fasta')
][1]

HA_NUCLEOTIDE = [
    str(rec.seq) for rec in SeqIO.parse(open('test-data/A(H1N1)pdm09/A-Brisbane-2-2018_nucleotide.fa'),'fasta')
][1]

### Load NA from reference file
NA_PROTEIN = [
    str(rec.seq) for rec in SeqIO.parse(open('test-data/A(H1N1)pdm09/A-Brisbane-2-2018_protein.fa'),'fasta')
][0]

NA_NUCLEOTIDE = [
    str(rec.seq) for rec in SeqIO.parse(open('test-data/A(H1N1)pdm09/A-Brisbane-2-2018_nucleotide.fa'),'fasta')
][0]

def get_codon_dict():
    """
    Gets a dictionary of amino acids (keys) and synonymous codon lists (values).

    Returns
    -------
    CODONDICT : dict
        A dictionary of amino acids (keys) and synonymous codon lists (values).

    @MelaniaNowicka, @JakubBartoszewicz
    """
    return CODONDICT


def get_aa_dict():
    """
    Gets a dictionary of codons (keys) and amino acids (values).

    Returns
    -------
    AADICT : dict
        A dictionary of codons (keys) and amino acids (values).

    @MelaniaNowicka, @JakubBartoszewicz
    """
    return AADICT


def mutation_pattern(mutation_type, pos_start, pos_end, AAfrom, AAto):
    """
    str * str * str * str - > str
    return the mutation pattern associated with a
    """
    if mutation_type == "M":
        return f"{AAfrom}{pos_start}{AAto}"
    elif mutation_type == "T":
        return f"{AAfrom}{pos_start}{AASTOP}"
    elif mutation_type == "D":
        return f"{AAfrom}{pos_start}-{pos_end}del"
    elif mutation_type == "I":
        return f"{AAto}{pos_start}-{pos_end}ins"
    else:
        infos = (pos_start, pos_end, AAfrom, AAto)
        raise ValueError(
            f"Error in mutation_pattern, mutation '{mutation_type}' is not known {infos}"
        )
