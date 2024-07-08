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
        description="Remove Fixed Substitutions: Delete substitutions that are fixed in the population from the table",
        epilog="Usage: python remove_fixed.py -i variants_with_phenotypes.tsv -f fixed_2019-2023_Wisconsin.csv [-o variants_with_phenotypes_without_fixed.tsv]",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="TSV with the set of variants detected over all sequences, with annotation",
    )
    parser.add_argument(
        "-f",
        "--fixed",
        required=True,
        help="CSV with mutations fixed in the population for the season",
    )
    parser.add_argument(
        "-s",
        "--season",
        default="22/23",
        help="CSV with mutations fixed in the population for the season",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="variants_with_phenotypes_without_fixed.tsv",
        help="Output table without substitutions fixed in the population",
    )
    args = parser.parse_args()

    in_file = pd.read_csv(args.input, sep='\t')
    fixed_file = pd.read_csv(args.fixed, sep=',')
    season = args.season
    out_file = args.output

    subset_fixed = fixed_file[fixed_file['Season'].isin([season])]

    res = in_file[~in_file['aa_pattern'].isin([list(subset_fixed['Mutation'])])]

    res.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    main()