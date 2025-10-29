import groovy.xml.XmlSlurper

process SUMMARY_XML {
    tag "${meta.id}_${meta.cycle_number}"
    label 'process_single'
    container 'community.wave.seqera.io/library/groovy:4_0_24--bd5a0401545c33a6'

    input:
    tuple val(meta), path(xml)

    output:
    path "*.tsv", emit: output

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.cycle_number}"
    """
    #!/usr/bin/env groovy
    import groovy.xml.XmlSlurper

    def check              = '\u2705'
    def cross              = '\u274C'
    def output_xml         = [["variable_name", "value", "expected", "check"]]

    def xs = new XmlSlurper().parse(new File("${xml}"))

    def tile_size = xs.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@SizeX != '' && node.@SizeY != ''
        }
        .collect {
            node -> [node.@SizeX.toInteger(), node.@SizeY.toInteger()]
        }.toSet()

    if (tile_size == null || tile_size[0] == null || tile_size[1] == null){
        output_xml.add(
            ["SizeX|SizeY", tile_size.toString(), "Same Integer", cross]
        )
    }
    else{
        output_xml.add(
            ["SizeX|SizeY", tile_size.toString(), "Same Integer", check]
        )
    }

    def pixels = xs.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@PhysicalSizeX != '' && node.@PhysicalSizeY != ''
        }
        .collect {
            node -> [node.@PhysicalSizeX.toDouble(), node.@PhysicalSizeY.toDouble()]
        }

    if (pixels == null || pixels.toSet().size() != 1 || pixels.toSet()[0][0] == null || pixels.toSet()[0][1] == null
        || (pixels[0][0]).round(3) != (pixels[0][1]).round(3)) {
        output_xml.add(
            ["PhysicalSizeX|PhysicalSizeY", pixels.toString(), "Numbers that are equal within 3 DP", cross]
        )
    }
    else {
        output_xml.add(
            ["PhysicalSizeX|PhysicalSizeY", pixels.toString(), "Numbers that are equal within 3 DP", check]
        )
    }

    def n_channels = xs.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@SizeC != ''
        }
        .collect { node -> node.@SizeC.toInteger() }

    if (n_channels == null || n_channels.toSet().size() != 1 || n_channels[0] == 0) {
        output_xml.add(
            ["SizeC", n_channels.toString(), "Consistent > 0 numbers", cross]
        )
    }
    else {
        output_xml.add(
            ["SizeC", n_channels.toString(), "Consistent > 0 numbers", check]
        )
    }

    def size_units = xs.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@PhysicalSizeXUnit != '' && node.@PhysicalSizeYUnit != ''
        }
        .collect { node -> [node.@PhysicalSizeXUnit.toString(), node.@PhysicalSizeYUnit.toString()] }

    if (size_units == null || size_units.toSet().size() != 1 || size_units.flatten().toSet().size() != 1 ||
        size_units.flatten().toSet()[0] == null ||
        !(size_units.flatten().toSet()[0] in ["mm", "cm", "um", "µm", "reference_frame"])) {
        output_xml.add(
            ["PhysicalSizeXUnit|PhysicalSizeYUnit", size_units.toString(), "Consistent units (mm, cm, um, µm, reference_frame)", cross]
        )
    }
    else {
        output_xml.add(
            ["PhysicalSizeXUnit|PhysicalSizeYUnit", size_units.toString(), "Consistent units (mm, cm, um, µm, reference_frame)", check]
        )
    }

    def pixel_datatype = xs.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@Type
        }
        .collect { node -> node.@Type.toString() }

    if (pixel_datatype == null || pixel_datatype.toSet().size() != 1 ||
        !(pixel_datatype.toSet()[0] ==~ /[u]?int(8|16|32)|float|double/)  // There are more bit|complex|double-complex
        ) {
        output_xml.add(
            ["Type", pixel_datatype.toString(), "Consistent valid datatypes (uint8, float16...)", cross]
        )
    }
    else {
        output_xml.add(
            ["Type", pixel_datatype.toString(), "Consistent valid datatypes (uint8, float16...)", check]
        )
    }

    def exposure_time = xs.'**'.findAll {
            node -> node.name() == 'Plane'
        }
        .collect {
            node ->
                [
                    node.@ExposureTime.toDouble(),
                    node.@ExposureTimeUnit.toString()
                ]
        }
        .toSet()

    if (exposure_time == null || exposure_time.size() != 1 || exposure_time[0] == null || exposure_time[1] == null || exposure_time[1] == "") {
        output_xml.add(
            ["ExposureTime|ExposureTimeUnit", exposure_time.toString(),
            "Consistent valid exposure time and units", cross]
        )
    }
    else {
        output_xml.add(
            ["ExposureTime|ExposureTimeUnit", exposure_time.toString(),
            "Consistent valid exposure time and units", check]
        )
    }

    def output_filename = "${prefix}_xml_mqc.tsv"
    new File(output_filename).text = output_xml*.join("\\t").join("\\n")
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

    def output_filename = prefix + "_markersheet_mqc.tsv"
    def f1              = task.workDir.resolve(output_filename)
    f1.text             = output*.join("\t").join("\n")
}

process SUMMARY_SAMPLESHEET {
    tag "$meta.id"
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

    def output_filename = prefix + "_samplesheet_mqc.tsv"
    def f1              = task.workDir.resolve(output_filename)
    f1.text             = output_samplesheet*.join("\t").join("\n")
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

    def output_filename = prefix + "_markersheet_mqc.tsv"
    def f1              = task.workDir.resolve(output_filename)
    f1.text             = output_markersheet*.join("\t").join("\n")
}
