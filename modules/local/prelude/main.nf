import groovy.xml.XmlSlurper

process PRELUDE {
    tag "$meta.id"
    label 'process_single'

    exec:
    val meta
    val markersheet
    val samplesheet
    val xml

    output:
    val meta           , emit: meta
    path output_file   , emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "summary"

    def check = '\u2705'
    def cross = '\u274C'
    def output = [["source", "variable_name", "value", "expected", "check"]]

    xml = new XmlSlurper().parse(new File(xml.toString()))

    tile_size = xml.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@SizeX != '' && node.@SizeY != ''
        }
        .collect {
            node -> [node.@SizeX.toInteger(), node.@SizeY.toInteger()]
        }.toSet()

    if (tile_size == null || tile_size[0] == null || tile_size[1] == null){
        output.append(
            ["xml", "SizeX|SizeY", tile_size.toString(), "Integer", cross]
        )
    }
    else{
        output.append(
            ["xml", "SizeX|SizeY", tile_size.toString(), "Integer", check]
        )
    }

    pixels = xml.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@PhysicalSizeX != '' && node.@PhysicalSizeY != ''
        }
        .collect {
            node -> [node.@PhysicalSizeX.toDouble(), node.@PhysicalSizeY.toDouble()]
        }

    if (pixels == null || pixels.toSet().size() != 1 || pixels.toSet()[0][0] == null || pixels.toSet()[0][1] == null
        || (pixels[0][0]).round(3) != (pixels[0][1]).round(3)) {
        output.append(
            ["xml", "PhysicalSizeX|PhysicalSizeY", pixels.toString(), "Numbers that are equal within 3 DP", cross]
        )
    }
    else {
        output.append(
            ["xml", "PhysicalSizeX|PhysicalSizeY", pixels.toString(), "Numbers that are equal within 3 DP", check]
        )
    }

    n_channels = xml.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@SizeC != ''
        }
        .collect { node -> node.@SizeC.toInteger() }

    if (n_channels == null || n_channels.toSet().size() != 1 || n_channels[0] == 0) {
        output.append(
            ["xml", "SizeC", n_channels.toString(), "Consistent > 0 numbers", cross]
        )
    }
    else {
        output.append(
            ["xml", "SizeC", n_channels.toString(), "Consistent > 0 numbers", check]
        )
    }

    size_units = xml.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@PhysicalSizeXUnit != '' && node.@PhysicalSizeYUnit != ''
        }
        .collect { node -> [node.@PhysicalSizeXUnit.toString(), node.@PhysicalSizeYUnit.toString()] }

    if (size_units == null || size_units.toSet().size() != 1 || size_units.flatten().toSet().size() != 1 ||
        size_units.flatten().toSet()[0] == null ||
        !(size_units.flatten().toSet()[0] in ["mm", "cm", "um", "µm", "reference_frame"])) {
        output.append(
            ["xml", "PhysicalSizeXUnit|PhysicalSizeYUnit", size_units.toString(), "Consistent units (mm, cm, um, µm, reference_frame)", cross]
        )
    }
    else {
        output.append(
            ["xml", "PhysicalSizeXUnit|PhysicalSizeYUnit", size_units.toString(), "Consistent units (mm, cm, um, µm, reference_frame)", check]
        )
    }

    pixel_datatype = xml.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@Type
        }
        .collect { node -> node.@Type.toString() }

    if (pixel_datatype == null || pixel_datatype.toSet().size() != 1 ||
        !(pixel_datatype.toSet()[0] ==~ /[u]?(int|float)(8|16|32)(_t)?/)
        ) {
        output.append(
            ["xml", "Type", pixel_datatype.toString(), "Consistent valid datatypes (uint8, float16...)", cross]
        )
    }
    else {
        output.append(
            ["xml", "Type", pixel_datatype.toString(), "Consistent valid datatypes (uint8, float16...)", check]
        )
    }

    exposure_time = xml.'**'.findAll {
            node -> node.name() == 'Plane'
        }
        .collect {
            node ->
                [
                    'exposure_time': node.@ExposureTime.toDouble(),
                    'exposure_time_unit': node.@ExposureTimeUnit.toString()
                ]
        }
        .toSet()

    if (exposure_time == null || exposure_time.size() != 1 || exposure_time[0] == null || exposure_time[1] == null) {
        output.append(
            ["xml", "ExposureTime|ExposureTimeUnit", exposure_time.toString(),
            "Consistent valid exposure time and units", cross]
        )
    }
    else {
        output.append(
            ["xml", "ExposureTime|ExposureTimeUnit", exposure_time.toString(),
            "Consistent valid exposure time and units", check]
        )
    }

    markersheet.map{
        entry ->
        entry.map{
            key, value ->
            temp = ["markersheet", key, value, "", ""]

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
                temp[4] = (value == null || value.toUpperCase() == value) ? cross : check
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

            output.append(temp)
        }
    }

    samplesheet.map {
            meta, _, _, _ -> meta
        }
        .map {
            key, value ->
            temp = ["samplesheet", key, value, "", ""]

            if(key in [
                "pixel_size",
                "channel_count",
                "tile_count",
                "pixel_size_x",
                "pixel_size_y",
                "cycle_number"
            ]) {
                temp[3] = "Number"
                temp[4] = (value == null || !(value instalceof Number)) ? cross : check
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

            output.append(temp)
        }

    def output_file = new File(${prefix}".csv")
    output_file.withWriter {
        w -> new CSVPrinter(w, CSVFormat.DEFAULT).printRecords(output)
    }
    //output_file.text = output*.join(",").join(System.lineSeparator())

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "summary"

    """
    touch "${prefix}.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prelude: \$(prelude --version)
    END_VERSIONS
    """
}
