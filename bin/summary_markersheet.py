#!/usr/bin/env python

import pandas
import sys
from io import StringIO

def pipeline():
    check              = '\u2705'
    cross              = '\u274C'
    warning            = '\u26A0'
    output_markersheet = [["row_id", "variable_name", "value", "expected", "check"]]
    counter = 0

    markersheet = pandas.read_json(sys.argv[1])

    for index, row in markersheet.iterrows():
        print(row)
        counter += 1
        for column, value in row.items():
            temp = [str(counter), str(column), str(value), None, None]

            if column in ["channel_number",
                        "cycle_number",
                        "exposure_time"]:
                temp[3] = "Number"
                temp[4] = cross if (pandas.isna(value) or
                                    not isinstance(value, (int, float, complex)) and
                                    not isinstance(value, bool)
                                    ) else check
            elif column in["excitation_wavelength",
                        "emission_wavelength"]:
                temp[3] = "Number"
                temp[4] = warning if (pandas.isna(value) or
                                    not isinstance(value, (int, float, complex)) and
                                    not isinstance(value, bool)
                                    ) else check

            elif column in ["marker_name"]:
                temp[3] = "Uppercase marker name"
                temp[4] = cross if (pandas.isna(value) or not value.isupper() ) else check

            elif column in ["filter",
                        "exposure",
                        "background",
                        "remove"]:
                temp[3] = "Boolean"
                temp[4] = warning if (pandas.isna(value) or not isinstance(value, bool)) else check

            elif column in ["exposure_time_unit"]:
                temp[3] = "Time unit"
                temp[4] = cross if (pandas.isna(value) or not isinstance(value, str) ) else check

            else:
                temp[3] = "?"
                temp[4] = cross

            output_markersheet.append(temp)

    for column in ["channel_number", "cycle_number", "marker_name"]:
        if column not in markersheet.columns:
            output_markersheet.append([str(counter), column, "MISSING", "?", cross])
            counter += 1

    for column in ["excitation_wavelength", "emission_wavelength", "exposure_time", "filter", "background",
                   "remove", "exposure_time_unit"]:
        if column not in markersheet.columns:
            output_markersheet.append([str(counter), column, "MISSING", "?", warning])
            counter += 1

    with open("{}_markersheet_mqc.tsv".format(sys.argv[2]), "w") as f:
        f.write("\n".join(["\t".join(x) for x in output_markersheet]))
        f.write("\n")


if __name__ == "__main__":
    pipeline()
