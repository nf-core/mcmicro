#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_IMAGE_FORMATS = (".tiff", ".tif")

    VALID_MARKER_FORMATS = ".csv"

    def __init__(
        self,
        sample_col="sample",
        first_col="cycle_number",
        second_col="channel_count",
        third_col="image_tiles",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            first_col (str): The name of the column that contains the cycles number.
            second_col (str): The name of the column that contains the number of channels.
            third_col (str): The name of the column that contains the image tiles file
                path (default "tiff").

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._first_col = first_col
        self._second_col = second_col
        self._third_col = third_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """

        '''
        self._validate_sample(row)
        print('*** done validating sample ***')
        self._validate_first(row)
        print('*** done validating first ***')
        self._validate_second(row)
        print('*** done validating second ***')
        self._validate_third(row)
        print('*** done validating third ***')
        self._seen.add((row[self._sample_col], row[self._first_col]))
        self.modified.append(row)
        '''

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("Sample input is required.")
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_first(self, row):
        """Assert that the cycle entry has the right format and exists"""
        if len(row[self._first_col]) <= 0:
            raise AssertionError("cycle required.")
        self._validate_cycle_format(row[self._first_col])

    def _validate_second(self, row):
        """Assert that the channel_count entry has the right format if it exists."""
        if len(row[self._second_col]) <= 0:
            raise AssertionError("channel_count required.")
        self._validate_channel_count_format(row[self._second_col])

    def _validate_third(self, row):
        """Assert that the image entry has the right format if it exists."""
        if len(row[self._first_col]) <= 0:
            raise AssertionError("Image required.")
        self._validate_image_format(row[self._third_col])

    def _validate_image_format(self, filename):
        """Assert that a given filename has image extension."""
        if not any(filename.endswith(extension) for extension in self.VALID_IMAGE_FORMATS):
            raise AssertionError(
                f"The image file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_IMAGE_FORMATS)}"
            )

    def _validate_cycle_format(self, cycle):
        """Assert that the cycle is an integer."""
        print(f'cycle is {cycle}')
        try:
            cycle = int(cycle)
        except Exception as err:
            print(err)
            print("cycle must be an integer")
            sys.exit(1)

    def _validate_channel_count_format(self, channel_count):
        """Assert that the channel_count is an integer."""
        print(f'channel_count is {channel_count}')
        try:
            channel_count = int(channel_count)
        except Exception as err:
            print(err)
            print("channel_count must be an integer")
            sys.exit(1)

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and image filename is unique.

        In addition to the validation, also rename all samples to have a suffix of _T{n}, where n is the
        number of times the same sample exist, but with different image files, e.g., multiple runs per experiment.

        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of channel and image must be unique.")
        seen = Counter()
        for row in self.modified:
            sample = row[self._sample_col]
            seen[sample] += 1
            # row[self._sample_col] = f"{sample}_T{seen[sample]}"


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect

def check_marker_sheet(file_in, file_out):
    import csv
    import collections

    marker_dict = collections.defaultdict(list)
    input_file = csv.DictReader(open(file_in))
    for row in input_file:
        for key,value in row.items():
            marker_dict[key].append(value)

    # uniqueness of marker name in marker sheet

    tmp_name_list = []
    for name in marker_dict['marker_name']:
        if name in tmp_name_list:
            raise Exception('Duplicate marker_name in marker sheet!')
        else:
            tmp_name_list.append(name)

    # uniqueness of (channel, cycle) tuple in marker sheet

    tmp_tup_list = []
    for i in range(len(marker_dict[list(marker_dict.keys())[0]])):
        curr_tup = (marker_dict['channel_number'][i], marker_dict['cycle_number'][i])
        if curr_tup in tmp_tup_list:
            raise Exception('Duplicate (channel_number, cycle_number) tuple in marker sheet!')
        else:
            tmp_tup_list.append(curr_tup)

    # cycle and channel are 1-based so 0 should throw an exception
    # cycle and channel cannot have skips and must be in order

    if int(marker_dict['channel_number'][0]) <= 0 or int(marker_dict['cycle_number'][0]) <= 0:
        raise Exception('channel_number and cycle number in the marker sheet are 1-based, so cannot be 0 or negative!')

    for i in range(1, len(marker_dict[list(marker_dict.keys())[0]])):
        if ( (marker_dict['channel_number'][i] != marker_dict['channel_number'][i-1]) and
             (int(marker_dict['channel_number'][i]) != int(marker_dict['channel_number'][i-1])+1) ):
            raise Exception('channel_number must be incresing without any gaps')
        if ( (marker_dict['cycle_number'][i] != marker_dict['cycle_number'][i-1]) and
             (int(marker_dict['cycle_number'][i]) != int(marker_dict['cycle_number'][i-1])+1) ):
            raise Exception('cycle_number must be incresing without any gaps')

    # TODO: this could be simplified to just returning the file_in atm, but leaving this here
    #   in case we want to make changes to the values in the block above
    with open(file_out, 'w') as fout:
        fout.write(','.join(list(marker_dict.keys())))
        fout.write("\n")
        # TODO: figure out a more pythonic way to get the following
        for i in range(len(marker_dict[list(marker_dict.keys())[0]])):
            curr_row_list = []
            for k in marker_dict:
                curr_row_list.append(marker_dict[k][i])
            curr_row_str = ','.join(curr_row_list) + "\n"
            fout.write(curr_row_str)


    '''
    required_columns = {"channel_number", "cycle_number", "marker_name", "excitation_wavelength", "emission_wavelength"}

    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        print('*** reader fieldnames ***')
        print(reader.fieldnames)
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()
    header = list(reader.fieldnames)
    # header.insert(1, "single_end")
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)
    '''


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_marker_sheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
