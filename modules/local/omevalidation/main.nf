import groovy.xml.XmlSlurper
process OMEVALIDATION {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), val(xmlPath)

    output:
    tuple val(meta), val(sample_meta), val(marker_meta)

    exec:
    xml = new XmlSlurper().parse(new File(xmlPath.toString()))

    /*
    SAMPLESHEET DATA ----------------------------------------------------------------------------------------------
    */

    tile_count = xml.'**'.findAll { node -> node.name() == "Image"}.size()

    if (tile_count < 2) {
        error 'Single image found in OMEXML metadata'
    }

    tile_size = xml.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@SizeX != '' && node.@SizeY != ''
        }
        .collect {
            node -> [node.@SizeX.toInteger(), node.@SizeY.toInteger()]
        }.toSet()

    if (tile_size.size() != 1) {
        error "Inconsistent tile sizes in images."
    }

    tile_size = tile_size[0]

    pixels = xml.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@PhysicalSizeX != '' && node.@PhysicalSizeY != ''
        }
        .collect {
            node -> [node.@PhysicalSizeX.toDouble(), node.@PhysicalSizeY.toDouble()]
        }

    if (pixels.toSet().size() != 1 || (pixels[0][0]).round(3) != (pixels[0][1]).round(3)) {
        error "Found non consistent pixels sizes in images."
    }

    n_channels = xml.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@SizeC != ''
        }
        .collect { node -> node.@SizeC.toInteger() }

    if (n_channels.toSet().size() != 1) {
        error "Found inconsistent number of channels in images."
    }

    size_units = xml.'**'.findAll {
            node -> node.name() == 'Pixels' && node.@PhysicalSizeXUnit != '' && node.@PhysicalSizeYUnit != ''
        }
        .collect { node -> [node.@PhysicalSizeXUnit.toString(), node.@PhysicalSizeYUnit.toString()] }

    if (size_units.toSet().size() != 1 || size_units.flatten().toSet().size() != 1) {
        error "Inconsistent pixels size unit in images."
    }

    // transform pixel size to microns
    s_units = size_units[0][0]

    if (size_units == 'mm'){
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

    pixel_datatype = xml.'**'.findAll {
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

    exposure_time = xml.'**'.findAll {
            node -> node.name() == 'Plane'
        }
        .collect {
            node ->
            if (node.@ExposureTime == '' || node.@ExposureTimeUnit == '') {
                return [
                    'cycle_number': meta.cycle_number,
                    //'channel_number': (meta.cycle_number - 1)*n_channels[0] + node.@TheC.toInteger() + 1, // channels on samplesheet start at 1
                    'channel_number': node.@TheC.toInteger() + 1,
                    'exposure_time': null,
                    'exposure_time_unit': null
                ]
            }
            else {
                return [
                    'cycle_number': meta.cycle_number,
                    //'channel_number': (meta.cycle_number - 1)*n_channels[0] + node.@TheC.toInteger() + 1, // channels on samplesheet start at 1
                    'channel_number': node.@TheC.toInteger() + 1,
                    'exposure_time': node.@ExposureTime.toDouble(),
                    'exposure_time_unit': node.@ExposureTimeUnit.toString()
                ]
            }
        }.toSet()

    if (exposure_time.size() != n_channels[0]) {
    //only needed inter cycle
        //println exposure_time
        error "Inconsistent number of exposure time entries, found " + exposure_time.size() + " expected " + n_channels[0]
    }

    sample_meta = [
                    'pixel_size': pixels,
                    'channel_count':n_channels[0],
                    'pixel_size_unit':'µm',
                    'pixel_datatype':pixel_datatype,
                    'tile_count': tile_count,
                    'tile_size_x': tile_size[0],
                    'tile_size_y': tile_size[1]
                ]

    marker_meta = exposure_time
}
