"""
Get the compound-target pairs dataset from ChEMBL using the given arguments.
"""

import logging
import sqlite3

import chembl_downloader

import arguments
import get_dataset


def main():
    """
    Call get_ct_pair_dataset to get the compound-target dataset using the given arguments.
    """
    args, calc_args, output_args = arguments.get_args()

    log_level = "DEBUG" if args.debug else "INFO"
    numeric_log_level = getattr(logging, log_level, None)
    assert isinstance(numeric_log_level, int), f"Invalid log level: %{args.log_level}"
    logging.basicConfig(level=numeric_log_level)

    if args.sqlite:
        logging.info(
            "Using provided sqlite3 path (%s) to connect to ChEMBL.", args.sqlite
        )
        assert args.chembl, "Please provide a ChEMBL version."
        with sqlite3.connect(args.sqlite) as chembl_con:
            get_dataset.get_ct_pair_dataset(
                chembl_con,
                calc_args,
                output_args,
            )
    else:
        logging.info("Using chembl_downloader to connect to ChEMBL.")
        if args.chembl_version is None:
            args.chembl_version = chembl_downloader.latest()

        with chembl_downloader.connect(version=args.chembl_version) as chembl_con:
            get_dataset.get_ct_pair_dataset(
                chembl_con,
                calc_args,
                output_args,
            )


if __name__ == "__main__":
    main()
