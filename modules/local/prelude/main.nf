import groovy.xml.XmlSlurper


process SUMMARY_XML {
    tag "${meta.id}_${meta.cycle_number}"
    label 'process_single'

    container "quay.io/biocontainers/python:3.13"

    input:
    tuple val(meta), path(xml)

    output:
    path "*.tsv", emit: output

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.cycle_number}"
    """
#! /usr/local/bin/python
import xml.etree.ElementTree as ET

def elementInTagConsistent(node, tag, attribute):
  if len(node.findall('.//{*}' + tag)) != len(node.findall('.//{*}' + tag + "[@{}]".format(attribute))):
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

check = '\u2705'
cross = '\u274C'
res = None

data = [["variable_name", "value", "expected", "check"]]

root = ET.parse('${xml}').getroot()

if not elementInTagConsistent(root, 'Pixels', 'SizeX') or not elementInTagConsistent(root, 'Pixels', 'SizeY') or \
root.findall('.//{*}Pixels[@SizeX]')[0].attrib['SizeX'] != root.findall('.//{*}Pixels[@SizeY]')[0].attrib['SizeY']:
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

valid_physical_units = ["mm", "cm", "um", "µm", "reference_frame"]

if not elementInTagConsistent(root, 'Pixels', 'PhysicalSizeXUnit') or \
not elementInTagConsistent(root, 'Pixels', 'PhysicalSizeYUnit') or \
root.findall('.//{*}Pixels[@PhysicalSizeXUnit]')[0].attrib['PhysicalSizeX'] != \
root.findall('.//{*}Pixels[@PhysicalSizeYUnit]')[0].attrib['PhysicalSizeY'] or \
root.findall('.//{*}Pixels[@PhysicalSizeXUnit]')[0].attrib['PhysicalSizeX'] not in valid_physical_units:
  res = cross
else:
  res = check

data.append(
[
  "PhysicalSizeXUnit|PhysicalSizeYUnit",
  getAllValuesFrom2Attrib(root, 'Pixels', 'PhysicalSizeXUnit', 'PhysicalSizeYUnit'),
  "Consistent units (mm, cm, um, µm, reference_frame)",
  res
]
)

valid_datatypes = ['uint8', 'uint16', 'uint32', 'int8', 'int16', 'int32', 'float', 'double']

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

if not elementInTagConsistent(root, 'Plane', 'ExposureTime') or \
not elementInTagConsistent(root, 'Plane', 'ExposureTimeUnit'):
  res = cross
else:
  res = check

data.append(
[
  "ExposureTime|ExposureTimeUnit",
  getAllValuesFrom2Attrib(root, 'Plane', 'ExposureTime', 'ExposureTimeUnit'),
  "Consistent valid exposure time and units",
  res
]
)

print(data)

with open("${prefix}" + "_xml_mqc.tsv", 'w') as f:
  f.write(
    '\\n'.join(['\\t'.join(x) for x in data])
  )
    """
}

process SUMMARY_MARKERSHEET_LITERAL {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), val(markersheet)

    output:
    path "*.tsv", emit: output

    when:
    task.ext.when == null || task.ext.when

    exec:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"

    def header = [
            "channel_number",
            "cycle_number",
            "excitation_wavelength",
            "emission_wavelength",
            "exposure_time",
            "marker_name",
            "filter",
            "exposure",
            "background",
            "remove",
            "exposure_time_unit"]

    def output = [header]

    markersheet.collect { m -> output.add( header.collect{ h -> m[h] ?: "" } ) }

    def output_file_markersheet = prefix + "_markersheet_mqc.tsv"
    def f1                      = task.workDir.resolve(output_file_markersheet)
    f1.text                     = output*.join("\t").join("\n")
}

process SUMMARY_SAMPLESHEET {
    tag "${meta.id}_${meta.cycle_number}"
    label 'process_single'

    input:
    tuple val(meta), val(samplesheet)

    output:
    path "*.tsv", emit: output

    when:
    task.ext.when == null || task.ext.when

    exec:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}_${meta.cycle_number}"

    def check              = '\u2705'
    def cross              = '\u274C'
    def output_samplesheet = [["row_id", "variable_name", "value", "expected", "check"]]
    def counter            = 0

    meta
        .each {
            key, value ->
            temp = [counter, key, value, "", ""]
            counter++

            if(key in [
                "pixel_size",
                "channel_count",
                "tile_count",
                "pixel_size_x",
                "pixel_size_y",
                "cycle_number"
            ]) {
                temp[3] = "Number"
                temp[4] = (value == null || !(value instanceof Number)) ? cross : check
            }
            else if (key in [
                "pixel_size_unit",
                "pixel_datatype"
            ]) {
                temp[3] = "Unit"
                temp[4] = (value == null || !(value instanceof String)) ? cross : check
            }
            else if (key in ["id"]) {
                temp[3] = "String"
                temp[4] = (value == null || !(value instanceof String)) ? cross : check
            }
            else {
                temp[3] = "?"
                temp[4] = cross
            }

            output_samplesheet.add(temp)
        }

    def output_file_samplesheet = prefix + "_samplesheet_mqc.tsv"
    def f1                  = task.workDir.resolve(output_file_samplesheet)
    f1.text                 = output_samplesheet*.join("\t").join("\n")
}

process SUMMARY_MARKERSHEET {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), val(markersheet)

    output:
    path "*.tsv", emit: output

    when:
    task.ext.when == null || task.ext.when

    exec:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"

    def check              = '\u2705'
    def cross              = '\u274C'
    def output_markersheet = [["row_id", "variable_name", "value", "expected", "check"]]
    def counter = 0
    markersheet
        .each { map ->
            map.each{ key, value ->
                temp = [counter, key, value, "", ""]
                counter++

                if (key in [
                    "channel_number",
                    "cycle_number",
                    "excitation_wavelength",
                    "emission_wavelength",
                    "exposure_time"
                ]) {
                    temp[3] = "Number"
                    temp[4] = (value == null || !(value instanceof Number)) ? cross : check
                }
                else if (key in [
                    "marker_name"
                ]){
                    temp[3] = "Uppercase marker name"
                    temp[4] = (value == null || value.toUpperCase() != value) ? cross : check
                }
                else if (key in [
                    "filter",
                    "exposure",
                    "background",
                    "remove"
                ]) {
                    temp[3] = "Boolean"
                    temp[4] = (value == null || !(value instanceof Boolean)) ? cross : check
                }
                else if (key in ["exposure_time_unit"]) {
                    temp[3] = "Time unit"
                    temp[4] = (value == null || !(value instanceof String)) ? cross : check
                }
                else {
                    temp[3] = "?"
                    temp[4] = cross
                }

                output_markersheet.add(temp)
            }
        }

    def output_file_markersheet = prefix + "_markersheet_mqc.tsv"
    def f1                      = task.workDir.resolve(output_file_markersheet)
    f1.text                     = output_markersheet*.join("\t").join("\n")
}
