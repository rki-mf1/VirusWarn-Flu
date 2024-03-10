#!/usr/bin/env python

from typing import Tuple

import re
import argparse
import numpy as np
import pandas as pd

def get_pos_start(
    text: str
) -> int:
    numbers = re.findall(r'\d+', text)
    return int(numbers[0]) if numbers else None

#def get_pos_end(
#        text: str
#    ) -> int:
#    numbers = re.findall(r'\d+', text)
#    return int(numbers[1]) if numbers else None

def change_numbering_h1n1(
    pattern: str,
) -> Tuple[str, str]:
    specified_gene = pattern.split(":")[0]
    mut = pattern.split(":")[1]

    if specified_gene == "HA1":
        # HA1 has a length of 327 in Nexclade
        return pattern, mut
    elif specified_gene == "HA2":
        # HA2 has a length of 222 in Nexclade
        numbering = mut[0] + str(int(mut[1:-1]) + 327) + mut[-1]
        return pattern, numbering
    elif specified_gene == "SigPep":
        # SigPep has a length of 17 in Nexclade
        return pattern, specified_gene
    else:
        raise KeyError("Unknown specified target gene!")
    
def change_numbering_h3n2(
    pattern: str,
) -> Tuple[str, str]:
    specified_gene = pattern.split(":")[0]
    mut = pattern.split(":")[1]

    if specified_gene == "HA1":
        # HA1 has a length of 329 in Nexclade
        return pattern, mut
    elif specified_gene == "HA2":
        # HA2 has a length of 221 in Nexclade
        numbering = mut[0] + str(int(mut[1:-1]) + 329) + mut[-1]
        return pattern, numbering
    elif specified_gene == "SigPep":
        # SigPep has a length of 16 in Nexclade
        return pattern, specified_gene
    else:
        raise KeyError("Unknown specified target gene!")
    
def change_numbering_vic(
    pattern: str,
) -> Tuple[str, str]:
    specified_gene = pattern.split(":")[0]
    mut = pattern.split(":")[1]

    if specified_gene == "HA1":
        # HA1 has a length of 183 in Nexclade
        return pattern, mut
    elif specified_gene == "HA2":
        # HA2 has a length of 224 in Nexclade
        numbering = mut[0] + str(int(mut[1:-1]) + 183) + mut[-1]
        return pattern, numbering
    elif specified_gene == "SigPep":
        # SigPep has a length of 15 in Nexclade
        return pattern, specified_gene
    else:
        raise KeyError("Unknown specified target gene!")

def build_dataframe(
    df_in: pd.DataFrame,
    column: str,
    subtype: str,
) -> pd.DataFrame:
    df_in[column] = df_in[column].astype(str)
    df_in[column] = df_in[column].apply(lambda x: x.split(","))
    df_temp = df_in.explode(column, ignore_index=True)
    df_temp.rename(columns={column: "aa_pattern", "seqName": "ID"}, inplace=True)

    df_temp = df_temp[df_temp['aa_pattern'] != 'nan']

    df_temp["target_gene"] = "HA"

    nextclade_aa_pattern = []
    aa_pattern = []
    if subtype == 'h1n1':
        for i in df_temp["aa_pattern"]:
            res1, res2 = change_numbering_h1n1(i)
            nextclade_aa_pattern.append(res1)
            aa_pattern.append(res2)
    elif subtype == 'h3n2':
        for i in df_temp["aa_pattern"]:
            res1, res2 = change_numbering_h3n2(i)
            nextclade_aa_pattern.append(res1)
            aa_pattern.append(res2)
    elif subtype == 'vic':
        for i in df_temp["aa_pattern"]:
            res1, res2 = change_numbering_vic(i)
            nextclade_aa_pattern.append(res1)
            aa_pattern.append(res2)
    df_temp["nextclade_aa_pattern"] = nextclade_aa_pattern
    df_temp["aa_pattern"] = aa_pattern

    if column == "aaSubstitutions":
        df_temp["variant_type"] = "M"
        df_temp["aa_pos_ref_start"] = df_temp['aa_pattern'].apply(get_pos_start)
        df_temp["aa_pos_ref_end"] = df_temp["aa_pos_ref_start"]
        df_temp["variant_size"] = 1
    elif column == "aaDeletions":
        df_temp["variant_type"] = "D"
        df_temp["aa_pos_ref_start"] = df_temp['aa_pattern'].apply(get_pos_start)
        #df_temp["aa_pos_ref_end"] = df_temp['aa_pattern'].apply(get_pos_end)
        df_temp["aa_pos_ref_end"] = np.NaN
        df_temp["variant_size"] = np.NaN
    elif column == "aaInsertions":
        df_temp["variant_type"] = "I"
        df_temp["aa_pos_ref_start"] = df_temp['aa_pattern'].apply(get_pos_start)
        #df_temp["aa_pos_ref_end"] = df_temp['aa_pattern'].apply(get_pos_end)
        df_temp["aa_pos_ref_end"] = np.NaN
        df_temp["variant_size"] = np.NaN

    #df_temp["variant_size"] = (df_temp["aa_pos_ref_end"] - df_temp["aa_pos_ref_start"]) + 1

    # Set columns that can't be reproduced to NaN
    df_temp["variant_start_aligned"] = np.NaN
    df_temp["variant_end_aligned"] = np.NaN
    df_temp["aa_ref"] = np.NaN
    df_temp["aa_variant"] = np.NaN
    df_temp["aa_pos_query_start"] = np.NaN
    df_temp["aa_pos_query_end"] = np.NaN
    df_temp["nt_pos_ref_start"] = np.NaN
    df_temp["nt_pos_ref_end"] = np.NaN
    df_temp["nt_pos_query_start"] = np.NaN
    df_temp["nt_pos_query_end"] = np.NaN
    df_temp["nt_pattern"] = np.NaN
        
    # Row 1-3: VOCAL variant_table.tsv cols
    # Row 4: Added in FluWarnSystem for additional information cols
    # Row 5-6: Nextclade nextclade.csv kept cols
    if subtype == 'vic':
        subset_cols = [
            "variant_start_aligned", "variant_end_aligned", "variant_type", "aa_ref", "aa_variant", "variant_size",
            "aa_pos_ref_start", "aa_pos_ref_end", "aa_pos_query_start", "aa_pos_query_end", "nt_pos_ref_start", "nt_pos_ref_end", 
            "nt_pos_query_start", "nt_pos_query_end", "nt_pattern", "aa_pattern", "ID", "target_gene", 
            "nextclade_aa_pattern",
            "clade", "subclade", "totalAminoacidSubstitutions", "totalAminoacidDeletions",
            "totalAminoacidInsertions", "glycosylation", "nonACGTNs", "qc.overallStatus", "failedCdses", "warnings", "errors"
        ]
    else:
        subset_cols = [
            "variant_start_aligned", "variant_end_aligned", "variant_type", "aa_ref", "aa_variant", "variant_size",
            "aa_pos_ref_start", "aa_pos_ref_end", "aa_pos_query_start", "aa_pos_query_end", "nt_pos_ref_start", "nt_pos_ref_end", 
            "nt_pos_query_start", "nt_pos_query_end", "nt_pattern", "aa_pattern", "ID", "target_gene", 
            "nextclade_aa_pattern",
            "clade", "short-clade", "subclade", "totalAminoacidSubstitutions", "totalAminoacidDeletions",
            "totalAminoacidInsertions", "glycosylation", "nonACGTNs", "qc.overallStatus", "failedCdses", "warnings", "errors"
        ]

    df_temp_subset = df_temp[subset_cols]
    
    return df_temp_subset

def rearrange_table(
    table_in: str,
    subtype: str,
) -> pd.DataFrame:
    df_in = pd.read_csv(table_in, sep=";", index_col=0)

    # Column aaSubstitutions from nextclade.csv
    mut_df = build_dataframe(df_in, "aaSubstitutions", subtype)

    # Column aaDeletions from nextclade.csv
    del_df = build_dataframe(df_in, "aaDeletions", subtype)

    # Column aaInsertions from nextclade.csv
    in_df = build_dataframe(df_in, "aaInsertions", subtype)

    df_out = pd.concat([mut_df, del_df, in_df])
    
    return df_out


def main():
    """
    Main function

    Returns
    ---------
    None
    """
    parser = argparse.ArgumentParser(
        description="Nextclade: Rearrange Nextclade results for FluWarnSystem",
        epilog="Usage: python nextclade.py -i nextclade.csv [-o variant_table.tsv]",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="CSV with mutations generated with Nextclade",
    )
    parser.add_argument(
        "--subtype",
        required=True,
        help="Subtype of the input data",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="variant_table.tsv",
        help="Output table with the set of variants detected over all sequences",
    )
    parser.add_argument(
        "--sigpep",
        default="sigpep_table.tsv",
        help="Output table with the set of signal peptide mutations detected over all sequences",
    )
    args = parser.parse_args()

    table_in = args.input
    subtype = args.subtype
    table_out = rearrange_table(table_in, subtype)

    table_out_sigpep = table_out[table_out['aa_pattern'] == 'SigPep']
    table_out_variant = table_out[table_out['aa_pattern'] != 'SigPep']

    table_out_sigpep.to_csv(args.sigpep, sep="\t")
    table_out_variant.to_csv(args.output, sep="\t")

if __name__ == "__main__":
    main()