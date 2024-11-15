#!/usr/bin/env python

"""
Annotate the list of mutations found in a dataframe

@hrichard
"""

import sys
import argparse
from typing import List
import pandas as pd


def aggregate_mutation_table(
    df_annot: pd.DataFrame,
    group_cols: List = ["gene", "amino acid", "type"],
    ID_col: str = "ID",
    comment_col: str = "comment",
    sep: str = ";",
) -> pd.DataFrame:
    """
    Aggregates information from a table of mutations according to group_cols and merges information from infos_cols

    Args:
        df_annot (DataFrame): DataFrame with the Mutations of Concern
        group_cols (List, optional): Columns from df_annot to aggregate. Defaults to ["gene", "amino acid", "type"].
        ID_col (str, optional): Name of the column including IDs. Defaults to ID.
        comment_col (str, optional): Name of the column including comments. Defaults to comment.
        sep (str, optional): The used seperator. Defaults to ;.

    Returns:
        DataFrame: Containing the aggregated mutations of annot_file
    """
    df_aggregated = df_annot.groupby(group_cols).agg(
        ID_list=pd.NamedAgg(
            column=ID_col, aggfunc=lambda x: sep.join([str(a) for a in x.tolist()])
        ),
        comment_list=pd.NamedAgg(
            column=comment_col, aggfunc=lambda x: sep.join([str(a) for a in x.tolist()])
        ),
    )
    return df_aggregated.reset_index()


def merge_variants_annotation(
    df_mutations: pd.DataFrame,
    df_annot: pd.DataFrame,
    roi_table: str,
    mut_merge: List = ["target_gene", "aa_pattern"],
    annot_merge: List = ["gene", "amino acid"],
) -> pd.DataFrame:
    """
    Merges DataFrames with result from Step 1 and with MOCs. Reports as well if the mutation is in a ROI.

    Args:
        df_mutations (DataFrame): Results from variant_table.tsv in a DataFrame
        df_annot (DataFrame): DataFrame with the Mutations of Concern
        roi_table (str): The path to the csv table including the Regions of Interest.
        mut_merge (List, optional): Columns from df_mutations for merge. Defaults to ["target_gene", "aa_pattern"].
        annot_merge (List, optional): Columns from df_annot for merge. Defaults to ["gene", "amino acid"].

    Returns:
        DataFrame: Tuple containing the aligned reference sequence and query sequence
    """
    if 'variant_type' in df_mutations.columns and 'type' in df_annot.columns:
        print("There is a type column in each of the DataFrame. Continue.")
    else:
        sys.exit("ERROR: There has to be a type column in each of the DataFrame!")

    df_merge_position = pd.merge(
        df_mutations,
        df_annot,
        left_on=mut_merge,
        right_on=annot_merge,
        how="left",
        validate="many_to_many",
    )

    df_merge_position = df_merge_position.fillna(value={"type": "NotAnnotated"}).rename(
        columns={"ID_list": "infos"}
    )

    ROI_temp = pd.read_csv(roi_table)
    #ROI_temp['aa_start'] = ROI_temp['aa_start'].astype(str) + ":" + ROI_temp['aa_end'].astype(str)
    ROI_temp = ROI_temp.drop(columns=['nt_start', 'nt_end', 'aa_end'])

    ROI = ROI_temp.rename(columns={'name': 'gene', 'aa_start': 'aa_position', 'fun': 'functional domain'})

    df_merge_region = pd.merge(
        df_mutations,
        ROI,
        left_on=["target_gene", "aa_pos_ref_start"],
        right_on=["gene", "aa_position"],
        how="inner",
    ).rename(columns={"functional domain": "infos"})

    return pd.concat([df_merge_position, df_merge_region], ignore_index=True)


def annotate_variant_table(
        df_variants: pd.DataFrame, 
        annot_file: str, 
        roi_table: str,
    ) -> pd.DataFrame:
    """
    Annotate the variants in the DataFrame using the table of annotations given in annot_file.
    The mutations in annot_file are aggregated.

    Args:
        df_variants (DataFrame): The DataFrame from variant_table.tsv covSonar built in Step 1.
        annot_file (str): The path to the tsv table including the Mutations of Concern.
        roi_table (str): The path to the csv table including the Regions of Interest.

    Returns:
        DataFrame: Containing the annotated variants 
    """
    if annot_file.endswith(".csv"):
        df_annot = pd.read_csv(annot_file)
    elif annot_file.endswith(".tsv"):
        df_annot = pd.read_csv(annot_file, sep="\t")
    else:
        raise TypeError("Not recognized file type")
    
    df_agg_annot = aggregate_mutation_table(df_annot)
    df_variant_with_annot = merge_variants_annotation(df_variants, df_agg_annot, roi_table)
    return df_variant_with_annot


def main():
    """
    Main function

    Returns
    ---------
    None
    """
    parser = argparse.ArgumentParser(
        description="Mutations2Functions: Intersect variant table with a set of sequences",
        epilog="Usage: python Mutations2Functions.py -i variant_table.tsv -a moc_table.tsv -r roi_table.csv -o [variants_with_phenotypes.tsv]",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", 
        "--input", 
        required=True, 
        help="Table of mutations (produced by Vocal)"
    )
    parser.add_argument(
        "-a",
        "--annotation",
        default="table_h1n1_mutations_annotation.tsv",
        help="Table with information about Mutations of Concern for the subtype",
    )
    parser.add_argument(
        "-r",
        "--roi_table",
        default="table_h1n1_roi.csv",
        help="Table with Regions of Interest for the subtype",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="variants_with_phenotypes.tsv",
        help="Output table with the set of variants detected over all sequences, with annotation",
    )
    parser.add_argument(
        "-L",
        "--largetable",
        action="store_true",
        help="Write all columns for the result table (there may be redundancies)",
    )

    args = parser.parse_args()

    tab_file = args.input
    annot_file = args.annotation
    roi_table = args.roi_table
    out_table = args.output

    # Reading in the files
    df_variants = pd.read_csv(tab_file, sep="\t")
    df_variant_with_annot = annotate_variant_table(df_variants, annot_file, roi_table).sort_values(
        [
            "ID", 
            "aa_pos_ref_start", 
            "variant_type", 
            "type",
        ]
    )

    columns_save = [
        "ID",
        "target_gene",
        "aa_pattern",
        "nt_pattern",
        "aa_pos_ref_start",
        "variant_type",
        "variant_size",
        "type",
        "infos",
    ] 

    if args.largetable:
        df_variant_with_annot.to_csv(out_table, sep="\t", index=False)
    else:
        df_variant_with_annot_small = df_variant_with_annot[columns_save]

        # Adds a row with type RegionOfInterest to prevent error in following R script 
        df_variant_with_annot_small.loc[len(df_variant_with_annot.index)] = [
            'ROI_ERROR_CATCH',	
            '',
            '',
            '',
            0,	
            '',	
            0,	
            'RegionOfInterest', 
            '',
        ] 
        df_variant_with_annot_small.loc[len(df_variant_with_annot.index)+1] = [
            'NotAnnotated_ERROR_CATCH',	
            '',
            '',
            '',
            0,	
            '',	
            0,	
            'NotAnnotated', 
            '',
        ] 
        
        df_variant_with_annot_small.to_csv(out_table, sep="\t", index=False)


if __name__ == "__main__":
    main()
