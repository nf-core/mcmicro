#!/usr/bin/env python

import xml.etree.ElementTree as ET
import sys


def checkExposureTimeWithinCycles(node, cycle):
    planes = node.findall(".//{*}Plane"+"[@TheC='{}'][@ExposureTime][@ExposureTimeUnits]".format(cycle))

    if len(set([x.attrib["ExposureTime"] for x in planes])) == 1 and \
            len(set([x.attrib["ExposureTimeUnits"] for x in planes])) == 1:
        return True
    else:
        return False

def elementInTagConsistent(node, tag, attribute):
    if len(node.findall('.//{*}' + tag)) != len(node.findall('.//{*}' + tag + "[@{}]".format(attribute))) or \
            len(node.findall('.//{*}' + tag + "[@{}]".format(attribute))) == 0:
        return False  # some elements missing tag

    if len(set([e.attrib[attribute] for e in node.findall('.//{*}' + tag + "[@{}]".format(attribute))])) != 1:
        return False  # different values
    return True

def getAllValuesFromAttrib(node, tag, attribute):
    return str([e.attrib[attribute] if attribute in e.attrib else '' for e in node.findall('.//{*}' + tag)])

def getAllValuesFrom2Attrib(node, tag, attribute1, attribute2):
    return str([
        [e.attrib[attribute1] if attribute1 in e.attrib else '', e.attrib[attribute2] if attribute2 in e.attrib else '']
        for e in node.findall('.//{*}' + tag)
    ])

def pipeline():
    check         = '\u2705'
    cross         = '\u274C'
    warn          = '\u26A0'
    res           = None
    MAX_TILE_SIZE = 2048

    data = [["variable_name", "value", "expected", "check"]]
    variables = [["variable_name", "cycle", "value"]]

    root = ET.parse(sys.argv[1]).getroot()

    if not elementInTagConsistent(root, 'Pixels', 'SizeX') or \
            not elementInTagConsistent(root, 'Pixels', 'SizeY'):
        res = cross
    else:
        res = check

    data.append(
        [
            'SizeX|SizeY',
            getAllValuesFrom2Attrib(root, 'Pixels', 'SizeX', 'SizeY'),
            'Same Integer',
            res
        ]
    )

    for i, sizex in enumerate(root.findall('.//{*}Pixels[@SizeX]')):
        variables.append(["SizeX", str(i), sizex.attrib["SizeX"]])

    for i, sizey in enumerate(root.findall('.//{*}Pixels[@SizeY]')):
        variables.append(["SizeY", str(i), sizey.attrib["SizeY"]])

    if len(root.findall('.//{*}Pixels[@SizeX]')) > 0 and \
       int(root.findall('.//{*}Pixels[@SizeX]')[0].attrib['SizeX']) > MAX_TILE_SIZE:
        res = warn
    else:
        res = check

    data.append(
        [
            'Max tile size',
            root.findall('.//{*}Pixels[@SizeX]')[0].attrib['SizeX'],
            'Number must be smaller than {}'.format(MAX_TILE_SIZE),
            res
        ]
    )

    if not elementInTagConsistent(root, 'Pixels', 'PhysicalSizeX') or \
            not elementInTagConsistent(root, 'Pixels', 'PhysicalSizeY') or \
            int(float(root.findall('.//{*}Pixels[@PhysicalSizeX]')[0].attrib['PhysicalSizeX']) * 1000) != \
            int(float(root.findall('.//{*}Pixels[@PhysicalSizeY]')[0].attrib['PhysicalSizeY']) * 1000):
        res = cross
    else:
        res = check

    data.append(
        [
            "PhysicalSizeX|PhysicalSizeY",
            getAllValuesFrom2Attrib(root, 'Pixels', 'PhysicalSizeX', 'PhysicalSizeY'),
            "Numbers that are equal within 3 DP",
            res
        ]
    )

    for i, physicalsizex in enumerate(root.findall('.//{*}Pixels[@PhysicalSizeX]')):
        variables.append(["PhysicalSizeX", str(i), str(int(float(physicalsizex.attrib["PhysicalSizeX"]) * 1000)/1000)])
    for i, physicalsizey in enumerate(root.findall('.//{*}Pixels[@PhysicalSizeY]')):
        variables.append(["PhysicalSizeY", str(i), str(int(float(physicalsizey.attrib["PhysicalSizeY"]) * 1000)/1000)])

    if not elementInTagConsistent(root, 'Pixels', 'SizeC') or \
            root.findall('.//{*}Pixels[@SizeC]')[0].attrib['SizeC'] == 0:
        res = cross
    else:
        res = check

    data.append(
        [
            'SizeC',
            getAllValuesFromAttrib(root, 'Pixels', 'SizeC'),
            'Consistent > 0 numbers',
            res
        ]
    )

    for i, sizec in enumerate(root.findall('.//{*}Pixels[@SizeC]')):
        variables.append(["SizeC", str(i), sizec.attrib["SizeC"]])

    valid_physical_units = ["mm", "cm", "um", "µm"]  # , "reference_frame"]

    if not elementInTagConsistent(root, 'Pixels', 'PhysicalSizeXUnit') or \
            not elementInTagConsistent(root, 'Pixels', 'PhysicalSizeYUnit') or \
            root.findall('.//{*}Pixels[@PhysicalSizeXUnit]')[0].attrib['PhysicalSizeXUnit'] != \
            root.findall('.//{*}Pixels[@PhysicalSizeYUnit]')[0].attrib['PhysicalSizeYUnit'] or \
            root.findall('.//{*}Pixels[@PhysicalSizeXUnit]')[0].attrib['PhysicalSizeXUnit'] not in valid_physical_units:
        res = cross
    else:
        res = check

    data.append(
        [
            "PhysicalSizeXUnit|PhysicalSizeYUnit",
            getAllValuesFrom2Attrib(root, 'Pixels', 'PhysicalSizeXUnit', 'PhysicalSizeYUnit'),
            "Consistent units (mm, cm, um, µm)",
            res
        ]
    )

    for i, physicalsizexunit in enumerate(root.findall('.//{*}Pixels[@PhysicalSizeXUnit]')):
        variables.append(["PhysicalSizeXUnit", str(i), physicalsizexunit.attrib["PhysicalSizeXUnit"]])
    for i, physicalsizeyunit in enumerate(root.findall('.//{*}Pixels[@PhysicalSizeYUnit]')):
        variables.append(["PhysicalSizeYUnit", str(i), physicalsizeyunit.attrib["PhysicalSizeYUnit"]])

    # other pipelines only compatible with uint8/16 we may add more later
    valid_datatypes = ['uint8', 'uint16'] #, 'uint32', 'int8', 'int16', 'int32', 'float', 'double']

    if not elementInTagConsistent(root, 'Pixels', 'Type') or \
            root.findall('.//{*}Pixels[@Type]')[0].attrib['Type'] not in valid_datatypes:
        res = cross
    else:
        res = check

    data.append(
        [
            "Type",
            getAllValuesFromAttrib(root, 'Pixels', 'Type'),
            "Consistent valid datatypes (uint8, float16...)",
            res
        ]
    )

    for i, datatype in enumerate(root.findall('.//{*}Pixels[@Type]')):
        variables.append(["Type", str(i), datatype.attrib["Type"]])

    for c in sorted(set([x.attrib["TheC"] for x in root.findall('.//{*}Plane[@TheC]')])):
        if not checkExposureTimeWithinCycles(root, c):
            res = warn
        else:
            res = check

        data.append(
            [
                "ExposureTime|ExposureTimeUnit cycle {}".format(c),
                str(
                    [
                        (x.attrib["ExposureTime"] if "ExposureTime" in x.attrib else "",
                         x.attrib["ExposureTimeUnits"] if "ExposureTimeUnits" in x.attrib else "")
                        for x in root.findall(".//{*}Plane" + "[@TheC='{}']".format(c))
                    ]
                ),
                "Consistent valid exposure time and units",
                res
            ]
        )

    for planes in root.findall('.//{*}Plane[@ExposureTime]'):
        variables.append(["ExposureTime", planes.attrib["TheC"], planes.attrib["ExposureTime"]])
    for planes in root.findall('.//{*}Plane[@ExposureTimeUnits]'):
        variables.append(["ExposureTimeUnits", planes.attrib["TheC"], planes.attrib["ExposureTimeUnits"]])

    with open("{}_xml_mqc.tsv".format(sys.argv[2]), 'w') as f:
        f.write(
            '\n'.join(['\t'.join(x) for x in data])
        )
        f.write("\n")

    with open("{}_variables.tsv".format(sys.argv[2]), 'w') as f:
        f.write(
            '\n'.join(['\t'.join(x) for x in variables])
        )
        f.write("\n")


if __name__ == "__main__":
    pipeline()
