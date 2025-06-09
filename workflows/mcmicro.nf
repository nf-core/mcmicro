/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

import groovy.io.FileType
import groovy.xml.XmlSlurper
import nextflow.Nextflow

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_mcmicro_pipeline'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { BASICPY                } from '../modules/nf-core/basicpy/main'
include { ASHLAR                 } from '../modules/nf-core/ashlar/main'
include { BACKSUB                } from '../modules/nf-core/backsub/main'
include { CELLPOSE               } from '../modules/nf-core/cellpose/main'
include { COREOGRAPH             } from '../modules/nf-core/coreograph/main'
include { DEEPCELL_MESMER        } from '../modules/nf-core/deepcell/mesmer/main'
include { SCIMAP_MCMICRO         } from '../modules/nf-core/scimap/mcmicro/main'
include { MCQUANT                } from '../modules/nf-core/mcquant/main'
include { OMEXTRACTOR            } from '../modules/nf-core/omextractor/main'
include { OMEVALIDATION          } from '../modules/local/omevalidation/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


process validateOmeXMLData {
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

    sample_meta = ['pixelsSize': pixels[0][0], 'nChannels':n_channels[0][0], 'pixelSizeUnit':'um', 'pixelDatatype':pixel_datatype[0][0]]
    marker_meta = ['exposureTime':exposure_time[0], 'exposureTimeUnits':exposure_time[1]]
}


workflow MCMICRO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input_cycle or --input_sample
    ch_markersheet // channel: markersheet read in from --marker_sheet

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // ch_samplesheet.multimap{meta, image_tiles, dfp, ffp -> meta: meta, image: image_tiles} | OMEXTRACTOR

    // val_data = validateOmeXMLData(OMEXTRACTOR.out.xml)  // TODO add inter meta checks and add values to samplesheet and markersheet

    ch_samplesheet.multimap{meta, image_tiles, dfp, ffp -> meta: meta, image: image_tiles} | OMEXTRACTOR
    ch_versions = ch_versions.mix(OMEXTRACTOR.out.versions)

    val_data = OMEXTRACTOR.out.xml | OMEVALIDATION

    // ch_samplesheet.join(val_data).map { original_meta, image_tiles, dfp, ffp, xml_data -> [xml_data + original_meta, image_tiles, dfp, ffp]}.dump(tag="ch_samplesheet (meta)").set { ch_samplesheet }
    ch_samplesheet.join(val_data)
        // .map { original_meta, image_tiles, dfp, ffp, xml_meta, xml_data, marker_data -> // i think if xml_data = originalmeta; then only original meta is used
         .map { original_meta, image_tiles, dfp, ffp, xml_meta, marker_data ->
                [xml_meta + original_meta, image_tiles, dfp, ffp]
        }
        .dump(tag="ch_samplesheet_meta")
        .set { ch_samplesheet }

    ch_markersheet.join(val_data)
        .map {
            'channel_number,cycle_number,marker_name,exposure,background,remove' ->
        }

    //
    // MODULE: BASICPY
    //
    if (params.illumination == 'basicpy') {
        ch_samplesheet
            .map{ meta, image_tiles, dfp, ffp ->
                [meta.subMap('id', 'cycle_number'), image_tiles]
            }
            .dump(tag: 'BASICPY in')
            | BASICPY
        ch_versions = ch_versions.mix(BASICPY.out.versions)
        ch_samplesheet = ch_samplesheet
            .map{ meta, image_tiles, dfp, ffp ->
                [meta.subMap('id', 'cycle_number'), image_tiles]
            }
            .join(BASICPY.out.profiles)
            .dump(tag: 'ch_samplesheet (after BASICPY)')
    }

    ch_samplesheet
        .map{ meta, image_tiles, dfp, ffp ->
            [[id: meta.id], [meta.cycle_number, image_tiles, dfp, ffp]]
        }
        // FIXME: pass groupTuple size: from samplesheet cycle count
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map{ meta, cycles -> [meta, *cycles.collect{ it[1..-1] }.transpose()]}
        .dump(tag: 'ASHLAR in')
        // flatten() handles list of empty-lists, turning it into a single empty list.
        .multiMap{ meta, images, dfps, ffps ->
            images: [meta, images]
            dfps: dfps.flatten()
            ffps: ffps.flatten()
        }
        | ASHLAR
    ch_versions = ch_versions.mix(ASHLAR.out.versions)

    // Run Background Correction
    if (params.backsub) {
        ch_backsub_markers = ch_markersheet
            .map { ['channel_number,cycle_number,marker_name,exposure,background,remove',
                it.collect{ it.channel_number + "," + it.cycle_number + "," + it.marker_name + "," + it.exposure + "," + it.background + "," + it.remove}] }
            .flatten()
            .map { it.replaceAll('(?<=,|^)null(?=,|$)', '') }
            .collectFile(name: 'markers_backsub.csv', sort: false, newLine: true)

        ASHLAR.out.tif
            .combine(ch_backsub_markers)
            .dump(tag: 'BACKSUB IN')
            .multiMap{ meta, image, marker ->
                image: [meta, image]
                markers: [meta, marker]
            }
            | BACKSUB

        post_registration = BACKSUB.out.backsub_tif
        ch_versions = ch_versions.mix(BACKSUB.out.versions)
    } else {
        post_registration = ASHLAR.out.tif
    }

    // Run Coreograph
    if (params.tma_dearray) {
        COREOGRAPH(post_registration)
        COREOGRAPH.out.cores
            .transpose()
            .map { meta, img -> [[id: meta.id + '_' + img.fileName.toString().tokenize('.')[0]], img]}
            .set { ch_segmentation_input }
    } else {
        ch_segmentation_input = post_registration
    }

    // Run Segmentation

    ch_masks = Channel.empty()

    ch_segmentation_input
        .multiMap{ meta, image ->
            img: [meta + [segmenter: 'mesmer'], image]
            membrane_img: [[:], []]
        }
        | DEEPCELL_MESMER
    ch_masks = ch_masks.mix(DEEPCELL_MESMER.out.mask)
    ch_versions = ch_versions.mix(DEEPCELL_MESMER.out.versions)

    ch_segmentation_input
        .multiMap{ meta, image ->
            image: [meta + [segmenter: 'cellpose'], image]
            model: params.cellpose_model
        }
        | CELLPOSE
    ch_masks = ch_masks.mix(CELLPOSE.out.mask)
    ch_versions = ch_versions.mix(CELLPOSE.out.versions)

    // Run Quantification

    // Generate markers.csv for mcquant with just the marker_name column, and
    // omitting rows removed by backsub.
    ch_mcquant_markers = Channel.of('marker_name')
        .concat(
            ch_markersheet
                .flatten()
                .filter{ row -> !(params.backsub && row.remove) }
                .map{ row -> '"' + row.marker_name + '"' }
        )
        .dump(tag: "MARKERS")
        .collectFile(name: 'markers.csv', sort: false, newLine: true)

    ch_segmentation_input
        .cross(ch_masks) { it[0]['id'] }
        .map{ t_ashlar, t_mask -> [t_mask[0], t_ashlar[1], t_mask[1]] }
        .combine(ch_mcquant_markers)
        .dump(tag: 'MCQUANT IN')
        .multiMap{ meta, image, mask, marker ->
            image: [meta, image]
            mask: [meta, mask]
            markers: [meta, marker]
        }
        | MCQUANT

    ch_versions = ch_versions.mix(MCQUANT.out.versions)

    /*
    // // Run Reporting
    SCIMAP_MCMICRO(MCQUANT.out.csv)
    ch_versions = ch_versions.mix(SCIMAP_MCMICRO.out.versions)
    */

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'mcmicro_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
