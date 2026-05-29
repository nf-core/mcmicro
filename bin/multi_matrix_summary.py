#!/usr/bin/env python

import pandas
import argparse
import sys

VALUE_CHECK = '\u2705'
VALUE_CROSS = '\u274C'
VALUE_WARN  = '\u26A0'


def parse_args():
    output = argparse.ArgumentParser()

    output.add_argument("--xml", type=str, nargs="+", required=True)
    #output.add_argument("--markersheet", type=str, nargs="+", required=True)
    output.add_argument("--samplesheet", type=str, nargs="+", required=True)
    output.add_argument("--merged", type=str, nargs="+", required=True)

    output.add_argument("--prefix", type=str, default="multi_matrix")

    return output.parse_args(sys.argv[1:])


def join_tests(data):
    # takes list(tuple(str, check|warn|err))
    global VALUE_WARN, VALUE_CROSS, VALUE_CHECK
    output = dict()

    entries = pandas.DataFrame(data)
    for i, entry in entries.groupby(0):
        if "ExposureTime" in i:  # reported per cycle
            continue

        if (entry[1] == VALUE_CROSS).any():
            output[i] = VALUE_CROSS
        elif (entry[1] == VALUE_WARN).any():
            output[i] = VALUE_WARN
        else:
            output[i] = VALUE_CHECK

    # handle exposuretime/units to be single entry and not per cycle
    if (entries[entries[0].str.contains("ExposureTime")][1] == VALUE_CROSS).any():
        output["ExposureTime|ExposureTimeUnits"] = VALUE_CROSS
    elif (entries[entries[0].str.contains("ExposureTime")][1] == VALUE_WARN).any():
        output["ExposureTime|ExposureTimeUnits"] = VALUE_WARN
    else:
        output["ExposureTime|ExposureTimeUnits"] = VALUE_CHECK

    return output


def extract_from_xml(file_paths):
    output = dict()
    for path in file_paths:
        data = pandas.read_csv(path, sep='\t')
        name, cycle = path.split('/')[-1].split("_")[:2]
        if name in output:
            output[name].extend([(row["variable_name"], row["check"]) for i, row in data.iterrows()])
        else:
            output[name] = [(row["variable_name"], row["check"]) for i, row in data.iterrows()]

    for key in output:
        output[key] = join_tests(output[key])

    return output


def extract_from_samplesheet(file_paths):
    output = dict()

    for path in file_paths:
        temp = pandas.read_csv(path, sep='\t')
        name = temp[temp['variable_name'] == 'id']['value']
        output[name.iat[0]] = { row["variable_name"]: row["check"] for i, row in temp.iterrows() }

    return output


def pipeline():
    args = parse_args()

    xml_data = extract_from_xml(args.xml)
    samplesheet_data = extract_from_samplesheet(args.samplesheet)
    #markersheet_data = extract_from_markersheet(args)

    merged = xml_data

    for key in samplesheet_data:
        print(xml_data, samplesheet_data, key)
        if key in merged:
            merged[key].update(samplesheet_data[key])
        else:
            merged[key] = samplesheet_data[key]

    merged = pandas.DataFrame().from_dict(merged, orient='index').T.fillna(VALUE_WARN)

    merged.to_csv("{}_matrix_summary_mqc.tsv".format(args.prefix), sep="\t")


if __name__ == "__main__":
    pipeline()
