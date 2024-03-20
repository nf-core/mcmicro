/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_mcmicro_pipeline'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'

include { BASICPY } from '../modules/nf-core/basicpy/main'
include { ASHLAR } from '../modules/nf-core/ashlar/main'
include { BACKSUB } from '../modules/nf-core/backsub/main'
include { CELLPOSE } from '../modules/nf-core/cellpose/main'
include { DEEPCELL_MESMER } from '../modules/nf-core/deepcell/mesmer/main'
include { MCQUANT } from '../modules/nf-core/mcquant/main'
include { SCIMAP_MCMICRO } from '../modules/nf-core/scimap/mcmicro/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Manually define inputs here
//image_tuple = tuple([ id:'image' ], '/home/florian/Documents/tmp_data_folder/cycif_tonsil_registered.ome.tif')
//marker_tuple = tuple([ id:'marker'], '/home/florian/Documents/tmp_data_folder/markers.csv')

workflow MCMICRO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet
        .multiMap {
            meta, data, marker, tissue ->
            registration: tuple(meta, data)
            marker: tuple(meta, marker)
        }
        .set { input }

    // Format input for BASICPY
    // data_path = ch_from_samplesheet
    //     .map(it->"${it[1]}/*.ome.tif")
    // raw_cycles = Channel.of([[id:"exemplar-001"],"/workspace/data/exemplar-001/raw/exemplar-001-cycle-06.ome.tiff"])

    //
    // MODULE: BASICPY
    //

    // BASICPY(raw_cycles)
    //ch_versions = ch_versions.mix(BASICPY.out.versions)

    // /*
    // if ( params.illumination ) {
    //     BASICPY(ch_images)
    //     ch_tif = BASICPY.out.fields
    //     ch_versions = ch_versions.mix(BASICPY.out.versions)

    //     ch_dfp = ch_tif.filter { file -> file.name.endsWith('.dfp.tiff') }
    //     ch_ffp = ch_tif.filter { file -> file.name.endsWith('.ffp.tiff') }
    // }
    // */

    ASHLAR(input.registration, [], [])
    ch_versions = ch_versions.mix(ASHLAR.out.versions)

    // // Run Background Correction
    // BACKSUB(ASHLAR.out.tif, ch_markers)
    // ch_versions = ch_versions.mix(BACKSUB.out.versions)

    // // Run Segmentation
    DEEPCELL_MESMER(ASHLAR.out.tif, [[:],[]])
    ch_versions = ch_versions.mix(DEEPCELL_MESMER.out.versions)

    // // Run Quantification
    MCQUANT(ASHLAR.out.tif,
            DEEPCELL_MESMER.out.mask,
            input.marker)
    ch_versions = ch_versions.mix(MCQUANT.out.versions)

    // // Run Reporting
    // SCIMAP_MCMICRO(MCQUANT.out.csv)
    // ch_versions = ch_versions.mix(SCIMAP_MCMICRO.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

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
