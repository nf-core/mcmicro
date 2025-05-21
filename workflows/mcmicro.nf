/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

import groovy.io.FileType
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
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


process validateOmeXMLData {
    input:
    path xmlPath
    
    output:
    val data

    script:
    xml = new XmlSlurper(xmlPath)
    pixels = xml.'**'.findAll { node -> node.name() == 'Pixels' && node.@PhysicalSizeX != '' && node.@PhysicalSizeY != ''}.collect { node -> [node.@PhysicalSizeX, node.@PhysicalSizeY] }
    if (pixels.toSet().size() != 1 || pixels[0][0].round(3) != pixels[0][1].round(3)) {
       error "Found non consistent pixels sizes in images."
    }

    n_channels = xml.'**'.findAll { node -> node.name() == 'Pixels' && node.@SizeC != '' }.collect { node -> node.@SizeC }
    if (n_channels.toSet().size() != 1) {
       error "Found inconsistent number of channels in images."
    }

    size_units = xml.'**'.findAll { node -> node.name() == 'Pixels' && node.@PhysicalSizeXUnit != '' && node.@PhysicalSizeYUnit != ''}.collect { node -> [node.@PhysicalSizeXUnit, node.@PhysicalSizeYUnit] }
    if (size_units.toSet().size() != 1 || size_units.flatten().toSet().size() != 1) {
       error "Inconsistent pixels size unit in images."
    }
    // TODO transform pixel size to microns

    switch(size_units[0][0]){
       case 'mm':
          pixels = pixels[0][0] / 1000
          break
       case 'cm':
          pixels = pixels[0][0] / 10000
          break
       case 'um':
          pixels = pixels[0][0]
          break
       case 'µm':
          pixels = pixels[0][0]
          break
    }

    pixel = pixel.round(3)

    pixel_datatype = xml.'**'.findAll { node -> node.name() == 'Pixels' && node.@Type}.collect { node -> node.@Type }
    if (pixels_datatype.toSet().size() != 1) {
       error "Inconsistent pixels datatype in images."
    }

    exposure_time = xml.'**'.findAll { node -> node.name() == 'Plane' && node.@ExposureTime != ''}.collect { node -> [node.@ExposureTime, node.@ExposureTimeUnit] }
    //only needed inter cycle
    //if (exposure_time.toSet().size() != 1) {
    //   error "Inconsistent exposure time"
    //}

    data = ['pixelsSize': pixels[0][0], 'nChannels':n_channels[0][0], 'pixelSizeUnit':size_units[0][0], 'pixelDatatype':pixel_datatype[0][0], 'exposureTime':exposure_time]


}


workflow MCMICRO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input_cycle or --input_sample
    ch_markersheet // channel: markersheet read in from --marker_sheet

    main:

    meta = ch_samplesheet.map{meta, image_tiles, dfp, ffp -> image_tiles} | OMEXTRACTOR

    val_data = validateOmeXMLData(meta)  // TODO add inter meta checks and add values to samplesheet and markersheet

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
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
                it.collect{ channel_number, cycle_number, marker_name, _1, _2, _3, exposure, background, remove ->
                    channel_number + "," + cycle_number + "," + marker_name + "," + exposure + "," + background + "," + remove}] }
            .flatten()
            .map { it.replace('[]', '') }
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

    // Generate markers.csv for mcquant with just the marker_name column.
    ch_mcquant_markers = ch_markersheet
        .flatMap{
            ['marker_name'] +
            it.collect{ _1, _2, marker_name, _4, _5, _6, _7, _8, _9 -> '"' + marker_name + '"' }
        }
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
