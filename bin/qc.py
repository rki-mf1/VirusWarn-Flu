#!/usr/bin/env python

import argparse
import pandas as pd

def main():
    """
    Main function

    Returns
    ---------
    None
    """
    parser = argparse.ArgumentParser(
        description="Quality Control: Generates table with QC info",
        epilog="Usage: python qc.py -i nextclade.csv [-o qc_table.tsv]",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="CSV with mutations generated with Nextclade",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="qc_table.tsv",
        help="Output table with the Quality Control stats",
    )
    args = parser.parse_args()

    table_in = args.input
    df_in = pd.read_csv(table_in, sep=";")

    subset_cols = [
        'seqName', 'qc.overallScore', 'qc.overallStatus', 
        'alignmentScore', 'alignmentStart', 'alignmentEnd', 'coverage', 'isReverseComplement', 
        'qc.missingData.missingDataThreshold', 'qc.missingData.score', 'qc.missingData.status', 'qc.missingData.totalMissing', 
        'qc.mixedSites.mixedSitesThreshold', 'qc.mixedSites.score', 'qc.mixedSites.status', 'qc.mixedSites.totalMixedSites', 
        'qc.privateMutations.cutoff', 'qc.privateMutations.excess', 'qc.privateMutations.score', 
        'qc.privateMutations.status', 'qc.privateMutations.total', 
        'qc.snpClusters.clusteredSNPs', 'qc.snpClusters.score', 'qc.snpClusters.status', 'qc.snpClusters.totalSNPs', 
        'qc.frameShifts.frameShifts', 'qc.frameShifts.totalFrameShifts', 'qc.frameShifts.frameShiftsIgnored', 
        'qc.frameShifts.totalFrameShiftsIgnored', 'qc.frameShifts.score', 'qc.frameShifts.status', 
        'qc.stopCodons.stopCodons', 'qc.stopCodons.totalStopCodons', 'qc.stopCodons.score', 'qc.stopCodons.status',
        'failedCdses', 'warnings', 'errors'
    ]

    qc_table = df_in[subset_cols]

    qc_table.rename(columns={"seqName": "ID"}, inplace=True)
    
    qc_table.to_csv(args.output, sep="\t")

if __name__ == "__main__":
    main()