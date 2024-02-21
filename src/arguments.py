"""
Dataclasses related to handling arguments, 
specifically arguments related to how to calculate or output the dataset. 
"""

import argparse

from dataclasses import dataclass


@dataclass(frozen=True)
class CalculationArgs:
    """
    Collection of arguments related to how to calculate the dataset.

    - chembl_version:         Version of ChEMBL for output file names
    - calculate_rdkit:        True if RDKit-based compound properties should be calculated
    - limit_to_literature:    Include only literature sources if True
    - limited_flag:           String version of limit_to_literature used in file names
    - min_nof_cpds_bf:        Minimum number of compounds per target for the BF subset
    - min_nof_cpds_b:         Minimum number of compounds per target for the B subset
    """

    chembl_version: str
    calculate_rdkit: bool
    limit_to_literature: bool
    limited_flag: str
    min_nof_cpds_bf: int
    min_nof_cpds_b: int


@dataclass(frozen=True)
class OutputArgs:
    """
    Collection of arguments related to how to output the dataset.

    - output_path:        Path to write output files to
    - delimiter:          Delimiter in csv-output
    - write_to_csv:       True if output should be written to csv
    - write_to_excel:     True if output should be written to excel
    - write_full_dataset: True if the full dataset should be written to output
    - write_bf:           True if subsets based on binding+functional data should be written to output
    - write_b:            True if subsets based on binding data only should be written to output
    """

    output_path: str
    delimiter: str
    write_to_csv: bool
    write_to_excel: bool
    write_full_dataset: bool
    write_bf: bool
    write_b: bool


def parse_args() -> argparse.Namespace:
    """
    Get arguments with argparse.

    :return: Populated argparse.Namespace
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Extract the compound-target pairs dataset from ChEMBL. \
            The full dataset plus filtering columns for binding vs. binding+functional data \
            will always be written to csv. \
            Additional outputs and output types can be chosen with the parameters below."
    )

    parser.add_argument(
        "--chembl",
        "-v",
        dest="chembl_version",
        metavar="<version>",
        type=str,
        default=None,
        help="ChEMBL version. \
            Latest version if None. \
            Required if a path to a SQLite database is provided, \
            i.e., if --sqlite is set. (default: None)",
    )
    parser.add_argument(
        "--sqlite",
        "-s",
        metavar="<path>",
        type=str,
        default=None,
        help="Path to SQLite database. \
            ChEMBL is downloaded as an SQLite database \
            and handled by chembl_downloader if None. (default: None)",
    )
    parser.add_argument(
        "--output",
        "-o",
        dest="output_path",
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
        dest="calculate_rdkit",
        action="store_true",
        help="Calculate RDKit-based compound properties.",
    )
    parser.add_argument(
        "--excel",
        dest="write_to_excel",
        action="store_true",
        help="Write the results to excel. Note: this may fail if the output is too large.",
    )
    parser.add_argument(
        "--BF",
        dest="write_bf",
        action="store_true",
        help="Write binding+functional data subsets.",
    )
    parser.add_argument(
        "--B", dest="write_b", action="store_true", help="Write binding data subsets."
    )
    parser.add_argument(
        "--debug", action="store_true", help="Log additional debugging information."
    )
    args = parser.parse_args()

    return args


def get_args() -> tuple[argparse.Namespace, CalculationArgs, OutputArgs]:
    """
    Get parsed and default arguments.

    :return: parserd arguments,
        arguments related to how to calculate the dataset as CalculationArgs,
        arguments related to how to output the dataset as OutputArgs
    :rtype: tuple[argparse.Namespace, CalculationArgs, OutputArgs]
    """
    args = parse_args()

    calc_args = CalculationArgs(
        chembl_version=args.chembl_version,
        calculate_rdkit=args.calculate_rdkit,
        limit_to_literature=not args.all_sources,
        # used in file names
        limited_flag="literature_only" if not args.all_sources else "all_sources",
        min_nof_cpds_bf=100,
        min_nof_cpds_b=100,
    )

    output_args = OutputArgs(
        output_path=args.output_path,
        delimiter=args.delimiter,
        # Always write the results to csv.
        write_to_csv=True,
        write_to_excel=args.write_to_excel,
        # Always write the full dataset plus filtering columns
        # for binding vs. binding+functional data.
        write_full_dataset=True,
        write_bf=args.write_bf,
        write_b=args.write_b,
    )

    return args, calc_args, output_args
