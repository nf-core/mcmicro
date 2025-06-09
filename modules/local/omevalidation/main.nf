process OMEVALIDATION {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(xmlPath)

    output:
    tuple val(meta), val(sample_meta), val(marker_meta)

    script:
    xml = new XmlSlurper(xmlPath)
    pixels = xml.'**'.findAll { node -> node.name() == 'Pixels' && node.@PhysicalSizeX != '' && node.@PhysicalSizeY != ''}.collect { node -> [node.@PhysicalSizeX.toDouble(), node.@PhysicalSizeY.toDouble()] }
    if (pixels.toSet().size() != 1 || (pixels[0][0]).round(3) != (pixels[0][1]).round(3)) {
       error "Found non consistent pixels sizes in images."
    }

    n_channels = xml.'**'.findAll { node -> node.name() == 'Pixels' && node.@SizeC != '' }.collect { node -> node.@SizeC.toInteger() }
    if (n_channels.toSet().size() != 1) {
       error "Found inconsistent number of channels in images."
    }

    size_units = xml.'**'.findAll { node -> node.name() == 'Pixels' && node.@PhysicalSizeXUnit != '' && node.@PhysicalSizeYUnit != ''}.collect { node -> [node.@PhysicalSizeXUnit.toString(), node.@PhysicalSizeYUnit.toString()] }
    if (size_units.toSet().size() != 1 || size_units.flatten().toSet().size() != 1) {
       error "Inconsistent pixels size unit in images."
    }
    // TODO transform pixel size to microns
    s_units = size_units[0][0]

    if (size_units == 'mm'){
      pixels = pixels[0][0] / 1000
    }
    else if (s_units == 'cm'){
      pixels = pixels[0][0] / 10000
    }
    else if (s_units == 'um' || s_units == 'µm'){
      pixels = pixels[0][0]
    }
    else{
      error "Invalid pixel size unit found."
    }

    pixel = pixel.round(3)

    pixel_datatype = xml.'**'.findAll { node -> node.name() == 'Pixels' && node.@Type}.collect { node -> node.@Type.toString() }
    if (pixels_datatype.toSet().size() != 1) {
       error "Inconsistent pixels datatype in images."
    }

    exposure_time = xml.'**'.findAll { node -> node.name() == 'Plane' && node.@ExposureTime != ''}.collect { node -> [node.@ExposureTime.toDouble(), node.@ExposureTimeUnit.toString()] }
    //only needed inter cycle
    //if (exposure_time.toSet().size() != 1) {
    //   error "Inconsistent exposure time"
    //}

    sample_meta = ['pixelsSize': pixel, 'nChannels':n_channels[0][0], 'pixelSizeUnit':'um', 'pixelDatatype':pixel_datatype[0][0]]
    marker_meta = ['exposureTime':exposure_time[0], 'exposureTimeUnits':exposure_time[1]]

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omevalidation: \$(omevalidation --version)
    END_VERSIONS
}
