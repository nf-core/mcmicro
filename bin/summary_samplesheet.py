#!/usr/bin/env python

import pandas
import sys
from io import StringIO


def pipeline():
    check              = '\u2705'
    cross              = '\u274C'
    warning            = '\u26A0'
    output_samplesheet = [["row_id", "variable_name", "value", "expected", "check"]]
    counter            = 0

    meta = None

    try:
        meta = pandas.read_json(sys.argv[1])
    except ValueError:
        meta = pandas.read_json(sys.argv[1], nrows=1, lines=True)

    for index, row in meta.iterrows():
        for column, value in row.items():
            temp = [str(counter), str(column), str(value), None, None]
            if column in ["pixel_size",
                    "channel_count",
                    "tile_count",
                    "tile_size_x",
                    "tile_size_y",
                    "cycle_number"]:
                temp[3] = "Number"
                temp[4] = cross if (pandas.isna(value) or
                                    not (isinstance(value, (int, float, complex)) and not isinstance(value, bool))
                                    ) else check
            elif column in ["pixel_size_unit", "pixel_datatype"]:
                temp[3] = "Unit"
                temp[4] = cross if (pandas.isna(value) or not isinstance(value, str)) else check
            elif column in ["id"]:
                temp[3] = "String"
                temp[4] = cross if (pandas.isna(value) or not isinstance(value, str)) else check
            else:
                temp[3] = "?"
                temp[4] = warning
            counter += 1
            output_samplesheet.append(temp)

    for column in ["pixel_size", "channel_count", "tile_count", "tile_size_x", "tile_size_y",
                    "cycle_number", "pixel_size_unit", "pixel_datatype"]:
        if column not in meta.columns:
            output_samplesheet.append([str(counter), column, "MISSING", "?", warning])
            counter += 1

    with open("{}_samplesheet_mqc.tsv".format(sys.argv[2]), "w") as f:
        f.write("\n".join(["\t".join(x) for x in output_samplesheet]))
        f.write("\n")


if __name__ == "__main__":
    pipeline()
