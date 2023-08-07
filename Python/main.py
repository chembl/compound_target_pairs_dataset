import argparse
import chembl_downloader
import logging
import sqlite3

import get_dataset

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Extract the compound-target dataset from ChEMBL.')
    parser.add_argument('--chembl_version', '-v',
                        metavar='<chembl_version>',
                        type=str,
                        default="",
                        help='ChEMBL version. Latest version if empty.')
    parser.add_argument('--sqlite_path', '-s',
                        metavar='<sqlite_path>',
                        type=str,
                        default="",
                        help='Path to sqlite database. Downloaded and handled by chembl_downloader if empty.')
    parser.add_argument('--output_path', '-o',
                        metavar='<output_path>',
                        type=str,
                        required=True,
                        help='Path to write output file(s) to.')
    parser.add_argument('--limit_to_literature', '-l',
                        metavar='<limit_to_literature>',
                        type=bool,
                        default=True,
                        help='Limit dataset to literature data. Note that average pchembl_values, ligand efficiencies,\
                            first_publication_cpd_target_pair, first_publication_cpd_target_pair_w_pchembl and first_publication_cpd \
                            will be based on literature data only if this is set to True.')
    parser.add_argument('--calculate_RDKit', '-r',
                        metavar='<calculate_RDKit>',
                        type=bool,
                        default=True,
                        help='calculate RDKit-based compound properties')
    parser.add_argument('--write_to_csv', '-c',
                        metavar='<write_to_csv>',
                        type=bool,
                        default=True,
                        help='write results to csv')
    parser.add_argument('--delimiter', '-d',
                        metavar='<delimiter>',
                        type=str,
                        default=";",
                        help='Delimiter in output csv-files.')
    # TODO: make this not fail even if the output is too large
    parser.add_argument('--write_to_excel', '-e',
                        metavar='<write_to_excel>',
                        type=bool,
                        default=False,
                        help='Write results to excel. Note: this may fail if the output is too large.')
    parser.add_argument('--write_full_dataset',
                        metavar='<write_full_dataset>',
                        type=bool,
                        default=True,
                        help='write full dataset plus filtering columns for binding vs. binding+functional data')
    parser.add_argument('--write_BF',
                        metavar='<write_BF>',
                        type=bool,
                        default=False,
                        help='write binding+functional data subsets')
    parser.add_argument('--write_B',
                        metavar='<write_B>',
                        type=bool,
                        default=False,
                        help='write binding data subsets')
    parser.add_argument('--log_level', 
                        metavar='<log_level>',
                        type=str,
                        default="INFO",
                        help='Level of detail of output. INFO: basic progress information; DEBUG: debugging information')
    args = parser.parse_args()

    assert(args.write_to_csv or args.write_to_excel), \
        "Please set either write_to_csv or write_to_excel to True. \
        Otherwise the results will not be saved."
    
    assert(args.write_full_dataset or args.write_BF or args.write_B), \
        "Please set either write_full_dataset, write_BF or write_B to True. \
        Otherwise the results will not be saved."

    numeric_log_level = getattr(logging, args.log_level.upper(), None)
    assert(isinstance(numeric_log_level, int)), f"Invalid log level: %{args.log_level}"
    logging.basicConfig(level=numeric_log_level)

    if args.chembl_version == "":
        args.chembl_version = chembl_downloader.latest()

    if args.sqlite_path == "":
        logging.info("Using chembl_downloader to connect to ChEMBL.")
        with chembl_downloader.connect(version=args.chembl_version) as chembl_con:
            get_dataset.get_ct_pair_dataset(chembl_con,
                                            args.chembl_version,
                                            args.output_path,
                                            args.limit_to_literature,
                                            args.calculate_RDKit,
                                            args.write_to_csv,
                                            args.write_to_excel,
                                            args.delimiter,
                                            args.write_full_dataset, args.write_BF, args.write_B)
    else:
        logging.info("Using provided sqlite3 path ({args.sqlite_path}) to connect to ChEMBL.")
        with sqlite3.connect(args.sqlite_path) as chembl_con:
            get_dataset.get_ct_pair_dataset(chembl_con,
                                            args.chembl_version,
                                            args.output_path,
                                            args.limit_to_literature,
                                            args.calculate_RDKit,
                                            args.write_to_csv,
                                            args.write_to_excel,
                                            args.delimiter,
                                            args.write_full_dataset, args.write_BF, args.write_B)
