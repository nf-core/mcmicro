/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { UPDATE_FROM_OME        } from '../subworkflows/local/update_from_ome'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_mcmicro_pipeline'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { BASICPY                } from '../modules/nf-core/basicpy/main'
include { ASHLAR                 } from '../modules/nf-core/ashlar/main'
include { BACKSUB                } from '../modules/nf-core/backsub/main'
include { CELLPOSE               } from '../modules/nf-core/cellpose/main'
include { MCCELLPOSE             } from '../modules/local/mccellpose/main'
include { COREOGRAPH             } from '../modules/nf-core/coreograph/main'
include { DEEPCELL_MESMER        } from '../modules/nf-core/deepcell/mesmer/main'
include { SCIMAP_MCMICRO         } from '../modules/nf-core/scimap/mcmicro/main'
include { MCQUANT                } from '../modules/nf-core/mcquant/main'
include { BFTOOLS_SHOWINF        } from '../modules/nf-core/bftools/showinf/main'
include { PRELUDE                } from '../subworkflows/local/prelude/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MCMICRO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input_cycle or --input_sample
    ch_markersheet // channel: markersheet read in from --marker_sheet
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_samplesheet.map{meta, image_tiles, dfp, ffp -> [meta, image_tiles]} | BFTOOLS_SHOWINF
    ch_versions = ch_versions.mix(BFTOOLS_SHOWINF.out.versions)

    metadata    = UPDATE_FROM_OME(ch_samplesheet, ch_markersheet, BFTOOLS_SHOWINF.out.xml)

    ch_samplesheet = metadata.samplesheet
    ch_markersheet = metadata.markersheet

    ch_samplesheet.dump(tag: "ch_samplesheet")
    ch_markersheet.dump(tag: "ch_markersheet")

    summaries = PRELUDE(ch_markersheet, ch_samplesheet, BFTOOLS_SHOWINF.out.xml)

    ch_multiqc_files = ch_multiqc_files.mix(summaries.output_file_mixed_matrix_summary)

    ch_has_no_errors = summaries.output_file_error_merged.map{
            meta, errors -> return tuple([meta['id'], meta['cycle_number']], errors.readLines().size() == 1)
        }.dump(tag:"METAHASNOERRORS")

    ch_samplesheet.dump(tag:"SampleBeforeFiltering")
    .combine(ch_has_no_errors.map{meta, has_no_error -> has_no_error})
    .filter{ meta, image_tiles, dfp, ffp, has_no_error -> has_no_error && !params.prelude}
    .map{meta, image_tiles, dfp, ffp, error -> [meta, image_tiles, dfp, ffp]}.dump(tag:"SampleAfterFiltering")
    .set{ch_samplesheet}

    //
    // MODULE: BASICPY
    //
    if (params.illumination == 'basicpy') {
        ch_samplesheet
            .map{ meta, image_tiles, dfp, ffp -> [meta, image_tiles] }
            .dump(tag: 'BASICPY in')
            | BASICPY
        ch_versions = ch_versions.mix(BASICPY.out.versions)
        ch_samplesheet = ch_samplesheet
            .map{ meta, image_tiles, dfp, ffp -> [meta, image_tiles] }
            .join(BASICPY.out.profiles)
            .dump(tag: 'ch_samplesheet (after BASICPY)')
    }

    ch_samplesheet
        .map{ meta, image_tiles, dfp, ffp ->
            [meta.subMap('id', 'pixel_size'), [meta.cycle_number, image_tiles, dfp, ffp]]
        }
        // FIXME: pass groupTuple size: from samplesheet cycle count
        .groupTuple(sort: { a, b -> a[0] <=> b[0] } )
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
            .map { meta, img -> [meta + [id: meta.id + '_' + img.fileName.toString().tokenize('.')[0]], img]}
            .set { ch_segmentation_input }
    } else {
        ch_segmentation_input = post_registration
    }

    // Run Segmentation

    ch_masks = channel.empty()

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

    ch_segmentation_input
        .multiMap{ meta, image ->
            image: [meta + [segmenter: 'mccellpose'], image]
        }
        | MCCELLPOSE
    ch_masks = ch_masks.mix(MCCELLPOSE.out.mask)
    ch_versions = ch_versions.mix(MCCELLPOSE.out.versions)

    // Run Quantification

    // Generate markers.csv for mcquant with just the marker_name column, and
    // omitting rows removed by backsub.
    ch_mcquant_markers = channel.of('marker_name')
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



    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'mcmicro_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'mcmicro'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    ch_has_no_errors.map{ if (!it[1]) error "QC Error found" }

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
