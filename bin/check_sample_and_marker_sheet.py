#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import collections
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

def check_sample_and_marker_sheet(input_path, marker_sheet_path):

    sample_dict = collections.defaultdict(list)
    input_file = csv.DictReader(open(input_path))
    for row in input_file:
        for key,value in row.items():
            sample_dict[key].append(value)
    
    if 'cycle_number' not in list(sample_dict.keys()):
        # no cycle_number in sample_sheet, so no additional validation
        return
    
    marker_dict = collections.defaultdict(list)
    input_file = csv.DictReader(open(marker_sheet_path))
    for row in input_file:
        for key,value in row.items():
            marker_dict[key].append(value)

    if set(sample_dict['cycle_number']) != set(marker_dict['cycle_number']):
        raise Exception('cycle_number values in sample and marker sheets must be 1:1 match.')

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "input",
        metavar="INPUT",
        type=Path,
        help="Tabular input sample sheet in CSV format.",
    )
    parser.add_argument(
        "marker_sheet",
        metavar="MARKER_SHEET",
        type=Path,
        help="Tablular input marker sheet in CSV format.",
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
    print("*** entering check_sample_and_marker_sheet: main ***")
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.input.is_file():
        logger.error(f"The given input file {args.input} was not found!")
        sys.exit(2)
    args.marker_sheet.parent.mkdir(parents=True, exist_ok=True)
    print('*** args ***')
    print(args)
    check_sample_and_marker_sheet(args.input, args.marker_sheet)


if __name__ == "__main__":
    sys.exit(main())
