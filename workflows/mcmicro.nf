/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

import groovy.io.FileType
import nextflow.Nextflow

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_mcmicro_pipeline'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { BASICPY                } from '../modules/nf-core/basicpy/main'
include { ASHLAR                 } from '../modules/nf-core/ashlar/main'
include { BACKSUB                } from '../modules/nf-core/backsub/main'
include { CELLPOSE               } from '../modules/nf-core/cellpose/main'
include { DEEPCELL_MESMER        } from '../modules/nf-core/deepcell/mesmer/main'
include { MCQUANT                } from '../modules/nf-core/mcquant/main'
include { SCIMAP_MCMICRO         } from '../modules/nf-core/scimap/mcmicro/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MCMICRO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input_cycle or --input_sample
    ch_markersheet // channel: markersheet read in from --marker_sheet

    main:

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

    // // Run Background Correction
    // BACKSUB(ASHLAR.out.tif, ch_markers)
    //BACKSUB(ASHLAR.out.tif, [[id: "backsub"], params.marker_sheet])
    /*
    if (params.backsub) {
        BACKSUB(ASHLAR.out.tif, [[id:"$ASHLAR.out.tif[0]['id']"], params.marker_sheet])
        ch_segmentation_input = BACKSUB.out.backsub_tif
        ch_versions = ch_versions.mix(BACKSUB.out.versions)
    } else {
    */
        ch_segmentation_input = ASHLAR.out.tif
    /*
    }
    */

    //ch_versions = ch_versions.mix(BACKSUB.out.versions)

    // Run Segmentation

    params.segmentation_list = params.segmentation.split(',') as List
    mesmer_out_mask = channel.empty()
    if ("mesmer" in params.segmentation_list) {
        DEEPCELL_MESMER(ch_segmentation_input, [[:],[]])
        ch_versions = ch_versions.mix(DEEPCELL_MESMER.out.versions)
        ch_segmentation_input
            .join(DEEPCELL_MESMER.out.mask)
            .dump(tag: 'MCQUANT in mesmer')
            .multiMap { it ->
                image: [it[0], it[1]]
                mask: [[id: it[0]['id'], segmenter: 'mesmer'], it[2]]
            }
            .set { ch_mesmer_out }
        mesmer_out_mask = ch_mesmer_out.mask
        ch_versions = ch_versions.mix(DEEPCELL_MESMER.out.versions)
    }

    cellpose_out_mask = channel.empty()
    if ("cellpose" in params.segmentation_list) {
        CELLPOSE( ch_segmentation_input, [] )
        ch_versions = ch_versions.mix(CELLPOSE.out.versions)
        ch_segmentation_input
            .join(CELLPOSE.out.mask)
            .dump(tag: 'MCQUANT in cellpose')
            .multiMap { it ->
                image: [it[0], it[1]]
                mask: [[id: it[0]['id'], segmenter: 'cellpose'], it[2]]
            }
            .set { ch_cellpose_out }
        cellpose_out_mask = ch_cellpose_out.mask
            .dump(tag: 'MASK')
        ch_versions = ch_versions.mix(CELLPOSE.out.versions)
    }

    // Run Quantification

    // Generate markers.csv for mcquant with just the marker_name column.
    ch_mcquant_markers = ch_markersheet
        .flatMap{
            ['marker_name'] +
            it.collect{ _1, _2, marker_name, _4, _5, _6 -> '"' + marker_name + '"' }
        }
        .collectFile(name: 'markers.csv', sort: false, newLine: true)

    ch_mcquant_masks_markers_key = channel.empty()
        .mix(mesmer_out_mask, cellpose_out_mask) // add new seg_out here when seg module added above
        .combine(ch_mcquant_markers)
        .map { meta, mask, marker -> [meta.id, meta, mask, marker ] }
        .dump(tag: 'MCQUANT MM KEY')

    ch_mcquant_ashlar_key = ASHLAR.out.tif.map { meta, img -> [meta.id, meta, img]}
        .dump(tag: 'MCQUANT IMG KEY')

    ch_mcquant_ashlar_key
        .combine(ch_mcquant_masks_markers_key, by: 0)
        .dump(tag: 'MCQUANT IN')
        .multiMap { it ->
            image: [it[1], it[2]]
            mask: [it[3], it[4]]
            markers: [it[3], it[5]]
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
            name: 'nf_core_pipeline_software_mqc_versions.yml',
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

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
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
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
