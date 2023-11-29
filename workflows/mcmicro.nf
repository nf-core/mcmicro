/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

//WorkflowMcmicro.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { BASICPY                     } from '../modules/nf-core/basicpy/main'
include { ASHLAR                      } from '../modules/nf-core/ashlar/main'
include { BACKSUB                     } from '../modules/nf-core/backsub/main'
include { CELLPOSE                    } from '../modules/nf-core/cellpose/main'
include { DEEPCELL_MESMER             } from '../modules/nf-core/deepcell/mesmer/main'
include { MCQUANT                     } from '../modules/nf-core/mcquant/main'
include { SCIMAP_MCMICRO              } from '../modules/nf-core/scimap/mcmicro/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

include { ILASTIK_PIXELCLASSIFICATION } from '../modules/nf-core/ilastik/pixelclassification/main'
include { ILASTIK_MULTICUT            } from '../modules/nf-core/ilastik/multicut/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Manually define inputs here
//image_tuple = tuple([ id:'image' ], '/home/florian/Documents/tmp_data_folder/cycif_tonsil_registered.ome.tif')
//marker_tuple = tuple([ id:'marker'], '/home/florian/Documents/tmp_data_folder/markers.csv')

// Info required for completion email and summary
def multiqc_report = []

workflow MCMICRO {

    // maybe better to just use input_check for validation and then parse 
    //   the channel created below from the sample sheet

    INPUT_CHECK(params.input)

//    input_check_channel
//        .multiMap 
//            { it -> 
//                ashlar: it[1][2]
//                foo: it[1][0]
//            }

    //input_check_channels.ashlar.view { "ashlar $it" }
    //input_check_channels.foo.view { "foo $it" }

    ch_versions = Channel.empty()

    ch_from_samplesheet = Channel.fromSamplesheet("input")
        .view { "all $it" } 
        .multiMap 
            { it -> 
                ashlar: [ [id:it[0]], [file(it[3])] ]
                foo: it[0]
            }

    ch_from_samplesheet.ashlar.view { "ashlar $it" }
    ch_from_samplesheet.foo.view { "foo $it" }

    // markerFile = [[id:"test_all" ], file("/workspace/data/cycif-tonsil-channels.csv")]
    marker_sheet = [[id:"test_all" ], file("/Users/robertyoung/DATA/exemplar/exemplar-001/markers.csv")]

    //raw_images = [ [ id:'test_all' ],[
    //            file("/workspace/data/ashlar/cycif-tonsil-cycle1.ome.tif"),
    //            file("/workspace/data/ashlar/cycif-tonsil-cycle2.ome.tif"),
    //            file("/workspace/data/ashlar/cycif-tonsil-cycle3.ome.tif")]
    //            ]

    //raw_images = [ [ id:'test_all' ],[
    //            file("/Users/robertyoung/DATA/cycif/tonsil/cycif-tonsil-cycle1.ome.tif")]
    //            ]

    //raw_images = [ [ id:'test_all' ],[
    //            file("/Users/robertyoung/DATA/exemplar/exemplar-001/raw/exemplar-001-cycle-06.ome.tiff")]
    //            ]

    raw_images = [ [ id:'test_all' ],[
                file("/Users/robertyoung/DATA/cycif/tonsil/cycif-tonsil-cycle1.ome.tif"),
                file("/Users/robertyoung/DATA/cycif/tonsil/cycif-tonsil-cycle2.ome.tif"),
                file("/Users/robertyoung/DATA/cycif/tonsil/cycif-tonsil-cycle3.ome.tif")]
                ]

    //dfp =file("/workspace/data/cycif-tonsil-dfp.ome.tif")
    //ffp =file("/workspace/data/cycif-tonsil-ffp.ome.tif")

    dfp =file("/Users/robertyoung/DATA/exemplar/exemplar-001/illumination/exemplar-001-cycle-06-dfp.tif")
    ffp =file("/Users/robertyoung/DATA/exemplar/exemplar-001/illumination/exemplar-001-cycle-06-ffp.tif")

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

    // ASHLAR(raw_images, dfp, ffp)
    //ASHLAR(raw_images, [], [])
    ASHLAR(ch_from_samplesheet.ashlar, [], [])
    ch_versions = ch_versions.mix(ASHLAR.out.versions)

    // // Run Background Correction
    // BACKSUB(ASHLAR.out.tif, ch_markers)
    // ch_versions = ch_versions.mix(BACKSUB.out.versions)

    // // Run Segmentation
    /*
    DEEPCELL_MESMER(ASHLAR.out.tif, [[:],[]])
    ch_versions = ch_versions.mix(DEEPCELL_MESMER.out.versions)
    */

    // fails
    // ILASTIK_PIXELCLASSIFICATION(ASHLAR.out.tif, 1)

    // // Run Quantification
    /*
    MCQUANT(ASHLAR.out.tif,
            DEEPCELL_MESMER.out.mask,
            markerFile)
    ch_versions = ch_versions.mix(MCQUANT.out.versions)
    */
    // // Run Reporting
    // SCIMAP_MCMICRO(MCQUANT.out.csv)
    // ch_versions = ch_versions.mix(SCIMAP_MCMICRO.out.versions)

    /*
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    */

    //
    // MODULE: MultiQC
    //
    /*
    workflow_summary    = WorkflowMcmicro.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMcmicro.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
    */
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
