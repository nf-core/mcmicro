import groovy.xml.XmlSlurper
import groovy.json.JsonOutput

process OMEVALIDATION {
    tag "$meta.id"
    label 'process_single'
    container 'community.wave.seqera.io/library/groovy:4_0_24--bd5a0401545c33a6'

    input:
    tuple val(meta), path(xml)

    output:
    tuple val(meta), path("sample_meta.json"), path("marker_meta.json")

    script:
    def prefix = "${meta.id}_${meta.cycle_number}"
    """
    #!/usr/bin/env groovy
    import groovy.xml.XmlSlurper
    import groovy.json.JsonOutput

    def xmlData = new XmlSlurper().parse(new File("${xml}"))

    /*
    SAMPLESHEET DATA ----------------------------------------------------------------------------------------------
    */

    def tile_count = xmlData.'**'.findAll { node -> node.name() == "Image"}.size()

    // TODO check if pre-stiched
    //if (tile_count < 2) {
    //   error 'Single image found in OMEXML metadata'
    //}

    def tile_size = xmlData.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@SizeX != '' && node.@SizeY != ''
        }
        .collect {
            node -> [node.@SizeX.toInteger(), node.@SizeY.toInteger()]
        }.toSet()

    if (tile_size.size() != 1) {
        error "Inconsistent tile sizes in images."
    }

    tile_size = tile_size[0]

    def pixels = xmlData.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@PhysicalSizeX != '' && node.@PhysicalSizeY != ''
        }
        .collect {
            node -> [node.@PhysicalSizeX.toDouble(), node.@PhysicalSizeY.toDouble()]
        }

    if (pixels.size() == 0) {
        error 'Images are missing pixel physical size metadata.'
    }
    if (pixels.toSet().size() != 1 || (pixels[0][0]).round(3) != (pixels[0][1]).round(3)) {
        error "Found non consistent pixels sizes in images."
    }

    def n_channels = xmlData.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@SizeC != ''
        }
        .collect { node -> node.@SizeC.toInteger() }

    if (n_channels.toSet().size() != 1 || n_channels[0] == 0) {
        error "Found inconsistent number of channels in images."
    }

    def size_units = xmlData.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@PhysicalSizeXUnit != '' && node.@PhysicalSizeYUnit != ''
        }
        .collect { node -> [node.@PhysicalSizeXUnit.toString(), node.@PhysicalSizeYUnit.toString()] }

    if (size_units.toSet().size() != 1 || size_units.flatten().toSet().size() != 1) {
        error "Inconsistent pixels size unit in images."
    }

    // transform pixel size to microns
    def s_units = size_units[0][0]

    if (s_units == 'mm'){
      pixels = pixels[0][0] / 1000
    }
    else if (s_units == 'cm'){
      pixels = pixels[0][0] / 10000
    }
    else if (s_units == 'um' || s_units == 'µm' || s_units == 'reference frame'){
      pixels = pixels[0][0]
    }
    else{
        error "Invalid pixel size unit found."
    }

    pixels = pixels.round(3)

    def pixel_datatype = xmlData.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@Type
        }
        .collect { node -> node.@Type.toString() }
        .toSet()

    if (pixel_datatype.size() != 1) {
        error "Inconsistent pixels datatype in images."
    }

    pixel_datatype = pixel_datatype[0]

    /*
    MARKERSHEET DATA ----------------------------------------------------------------------------------------------
    this needs to have channel specific data stored along so we can merge it reliably downstream
    */

    def exposure_time = xmlData.'**'.findAll {
            node -> node.name() == 'Plane'
        }
        .collect {
            node ->
            if (node.@ExposureTime == '' || node.@ExposureTimeUnit == '') {
                return [
                    'cycle_number': meta.cycle_number,
                    'channel_number': node.@TheC.toInteger() + 1,
                    'exposure_time': null,
                    'exposure_time_unit': null
                ]
            }
            else {
                return [
                    'cycle_number': meta.cycle_number,
                    'channel_number': node.@TheC.toInteger() + 1,
                    'exposure_time': node.@ExposureTime.toDouble(),
                    'exposure_time_unit': node.@ExposureTimeUnit.toString()
                ]
            }
        }.toSet()

    if (exposure_time.size() == 0) { //Plane is optional entry
        exposure_time = []

        for(int i = 0; i < n_channels[0]; i++)
            exposure_time.add([
                'cycle_number': meta.cycle_number,
                'channel_number': i+1,
                'exposure_time': null,
                'exposure_time_unit': null
            ])
    }

    def sample_meta = [
                    'pixel_size': pixels,
                    'channel_count':n_channels[0],
                    'pixel_size_unit':'µm',
                    'pixel_datatype':pixel_datatype,
                    'tile_count': tile_count,
                    'tile_size_x': tile_size[0],
                    'tile_size_y': tile_size[1]
                ]

    def marker_meta = exposure_time

    // Write outputs to JSON files
    new File("sample_meta.json").text = JsonOutput.toJson(sample_meta)
    new File("marker_meta.json").text = JsonOutput.toJson(marker_meta)
    """
}
