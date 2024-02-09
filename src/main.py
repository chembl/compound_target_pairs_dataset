import argparse
import chembl_downloader
import logging
import sqlite3

import get_dataset

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract the compound-target pairs dataset from ChEMBL. \
            The full dataset plus filtering columns for binding vs. binding+functional data will always be written to csv. \
            Additional outputs and output types can be chosen with the parameters below."
    )

    parser.add_argument(
        "--chembl",
        "-v",
        metavar="<version>",
        type=str,
        default=None,
        help="ChEMBL version. \
                            Latest version if None. Required if a path to a SQLite database is provided, i.e., if --sqlite is set. (default: None)",
    )
    parser.add_argument(
        "--sqlite",
        "-s",
        metavar="<path>",
        type=str,
        default=None,
        help="Path to SQLite database. \
                            ChEMBL is downloaded as an SQLite database and handled by chembl_downloader if None. (default: None)",
    )
    parser.add_argument(
        "--output",
        "-o",
        metavar="<path>",
        type=str,
        required=True,
        help="Path to write the output file(s) to. (required)",
    )
    parser.add_argument(
        "--delimiter",
        "-d",
        metavar="<delimiter>",
        type=str,
        default=";",
        help="Delimiter in output csv-files.  (default: ;)",
    )
    parser.add_argument(
        "--all_sources",
        action="store_true",
        help="If this is set, the dataset is calculated based on all sources in ChEMBL. \
                            This includes data from BindingDB which may skew the results. \
                            Default (not set): the dataset is calculated based on only literature data.",
    )
    parser.add_argument(
        "--rdkit",
        action="store_true",
        help="Calculate RDKit-based compound properties.",
    )
    parser.add_argument(
        "--excel",
        action="store_true",
        help="Write the results to excel. Note: this may fail if the output is too large.",
    )
    parser.add_argument(
        "--BF", action="store_true", help="Write binding+functional data subsets."
    )
    parser.add_argument("--B", action="store_true", help="Write binding data subsets.")
    parser.add_argument(
        "--debug", action="store_true", help="Log additional debugging information."
    )
    args = parser.parse_args()

    # Set arguments that are always true.
    # Write the results to csv.
    csv = True
    # Write the full dataset plus filtering columns for binding vs. binding+functional data.
    full_df = True

    if args.debug:
        log_level = "DEBUG"
    else:
        log_level = "INFO"
    numeric_log_level = getattr(logging, log_level, None)
    assert isinstance(numeric_log_level, int), f"Invalid log level: %{args.log_level}"
    logging.basicConfig(level=numeric_log_level)

    if args.sqlite:
        logging.info(
            f"Using provided sqlite3 path ({args.sqlite}) to connect to ChEMBL."
        )
        assert args.chembl, "Please provide a ChEMBL version."
        with sqlite3.connect(args.sqlite) as chembl_con:
            get_dataset.get_ct_pair_dataset(
                chembl_con,
                args.chembl,
                args.output,
                not args.all_sources,
                args.rdkit,
                csv,
                args.excel,
                args.delimiter,
                full_df,
                args.BF,
                args.B,
            )
    else:
        logging.info("Using chembl_downloader to connect to ChEMBL.")
        if args.chembl is None:
            args.chembl = chembl_downloader.latest()

        with chembl_downloader.connect(version=args.chembl) as chembl_con:
            get_dataset.get_ct_pair_dataset(
                chembl_con,
                args.chembl,
                args.output,
                not args.all_sources,
                args.rdkit,
                csv,
                args.excel,
                args.delimiter,
                full_df,
                args.BF,
                args.B,
            )
