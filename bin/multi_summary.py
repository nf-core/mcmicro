#!/usr/bin/env python

import argparse
import pandas
import sys


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("xmls", nargs="+", type=str, help="Xml summary files, (more than 1).")
    parser.add_argument("--output", type=str, help="Output file.", default="multi_summary")
    parser.add_argument("--errors", type=str, help="Output file with detailed errors.", default="multi_summary")

    return parser.parse_args(sys.argv[1:])


def check_consistency(data, variable_name):
    return len(data[data["variable_name"] == variable_name]["value"].value_counts()) == 1


def check_consistency_cycle(data, variable_name):
    for name, subset in data.groupby("cycle"):
        if not len(subset[subset["variable_name"] == variable_name].value_counts()) == 1:
            return False
    return True


def extract_errors(data, variable_name, is_cycle=False):
    output = list()
    if variable_name not in data["variable_name"]:
        for name in data["name"].unique():
            output.append((name, variable_name, "", "NOT PRESENT", "VALUE"))
    elif is_cycle:
        for cycle, subset in data[data["variable_name"] == variable_name].groupby("cycle"):
            reference = subset["value"]

            if len(reference) == 0:
                continue

            reference = reference.iat[0]

            for index, error in subset[subset["value"] != reference].iterrows():
                print(subset, reference, error)
                output.append((error["name"], variable_name, str(cycle), error["value"], reference))
    else:
        reference = data[data["variable_name"] == variable_name]["value"].iat[0]
        for index, error in data[(data["variable_name"] == variable_name) & (data["value"] != reference)].iterrows():
            output.append((error["name"], variable_name, str(error["cycle"]), error["value"], reference))

    return output


def pipeline():
    args = parse_args()
    data = load_summaries(args)

    output = [["Variable Name", "Expected", "Pass/Fail", "Comment"]]
    errors = [["Sample Name", "Variable Name", "Cycle", "Actual Value", "Expected Value"]]
    check = '\u2705'
    cross = '\u274C'
    warning = '\u26A0'

    if check_consistency_cycle(data, "ExposureTime"):
        output.append(["ExposureTime", "Consistent intra sample values", check, ""])
    else:
        output.append(["ExposureTime",
                       "Consistent intra sample values",
                       warning,
                       ""])

    if check_consistency_cycle(data, "ExposureTimeUnit"):
        output.append(["ExposureTimeUnit", "Consistent intra sample values", check, ""])
    else:
        output.append(["ExposureTimeUnit",
                       "Consistent intra sample values",
                       warning,
                       ""])

    if check_consistency(data, "SizeX"):
        output.append(["SizeX", "Consistent intra sample values", check, ""])
    else:
        output.append(["SizeX", "Consistent intra sample values", cross, ""])
        errors.extend(extract_errors(data, "SizeX"))

    if check_consistency(data, "SizeY"):
        output.append(["SizeY", "Consistent intra sample values", check, ""])
    else:
        output.append(["SizeY", "Consistent intra sample values", cross, ""])
        errors.extend(extract_errors(data, "SizeY"))

    if check_consistency(data, "PhysicalSizeX"):
        output.append(["PhysicalSizeX", "Consistent intra sample values", check, ""])
    else:
        output.append(["PhysicalSizeX", "Consistent intra sample values", cross, ""])
        errors.extend(extract_errors(data, "PhysicalSizeX"))

    if check_consistency(data, "PhysicalSizeY"):
        output.append(["PhysicalSizeY", "Consistent intra sample values", check, ""])
    else:
        output.append(["PhysicalSizeY", "Consistent intra sample values", cross, ""])
        errors.extend(extract_errors(data, "PhysicalSizeY"))

    if check_consistency(data, "PhysicalSizeXUnit"):
        output.append(["PhysicalSizeXUnit", "Consistent intra sample values", check, ""])
    else:
        output.append(["PhysicalSizeXUnit", "Consistent intra sample values", cross, ""])
        errors.extend(extract_errors(data, "PhysicalSizeXUnit"))

    if check_consistency(data, "PhysicalSizeYUnit"):
        output.append(["PhysicalSizeYUnit", "Consistent intra sample values", check, ""])
    else:
        output.append(["PhysicalSizeYUnit", "Consistent intra sample values", cross, ""])
        errors.extend(extract_errors(data, "PhysicalSizeYUnit"))

    if check_consistency(data, "Type"):
        output.append(["Type", "Consistent intra sample values", check, ""])
    else:
        output.append(["Type", "Consistent intra sample values", cross, ""])
        errors.extend(extract_errors(data, "Type"))

    for i, samples in data.groupby("name"):
        if (samples[samples["variable_name"] == "PhysicalSizeX"]["value"].values !=
                samples[samples["variable_name"] == "PhysicalSizeY"]["value"].values).any():
            errors.extend([[str(i), "PhysicalSizeX|PhysicalSizeY", "", "Non square pixels", "Non square pixels"]])

    with open(args.output + "_mqc.tsv", 'w') as f:
        f.write(
            '\n'.join(['\t'.join(x) for x in output])
        )
        f.write("\n")

    with open(args.errors + "_errors.tsv", 'w') as f:
        f.write(
            '\n'.join(['\t'.join(x) for x in errors])
        )
        f.write("\n")


def load_summaries(args):
    data = list()

    for summary in args.xmls:
        temp = pandas.read_csv(summary, sep="\t")
        temp = temp[["variable_name", "cycle", "value"]]  # could have more
        temp["name"] = summary.split("/")[-1]

        data.append(temp)

    return pandas.concat(data)


if __name__ == "__main__":
    pipeline()
