/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

import groovy.io.FileType

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
// include { INPUT_CHECK } from '../subworkflows/local/input_check'

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

    ch_versions = Channel.empty()

    ch_from_samplesheet = Channel.fromSamplesheet("input")
        .view { "all $it" }
        .multiMap
            { it ->
                ashlar: make_ashlar_input(it)
            }

    ch_from_samplesheet.ashlar.view { "ashlar $it" }

    ch_from_marker_sheet = Channel.fromSamplesheet("marker_sheet")
        .map { validate_marker_sheet(it) }

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

def make_ashlar_input(ArrayList samplesheet_row) {
    files = []
    def image_dir = new File(samplesheet_row[1])
    image_dir.eachFileRecurse (FileType.FILES) {
        if(it.toString().endsWith(".ome.tif")){
            files << file(it)
        }
    }


    ashlar_input = [[id:samplesheet_row[0]], files]

    return ashlar_input
}

marker_name_list = []
cycle_channel_tuple_list = []

def validate_marker_sheet(ArrayList sample_sheet_row) {
    // check marker name uniqueness
    if(marker_name_list.contains(sample_sheet_row[2])){
        throw new Exception("Error: duplicate marker name in marker sheet! Marker names must be unique.")
        System.exit(1)
    } else {
        marker_name_list.add(sample_sheet_row[2])
    }

    // check that (cycle, channel) tuple is unique
    curr_tuple = new Tuple(sample_sheet_row[0], sample_sheet_row[1])
    if(cycle_channel_tuple_list.contains(curr_tuple)){
        throw new Exception("Error: duplicate cycle_number & channel_number pair! cycle_number & channel_number pairs must be unique.")
        System.exit(1)
    } else if(sample_sheet_row[0] < 1 || sample_sheet_row[1] < 1) {
        throw new Exception("Error: cycle_number and channel_number are 1-based.  Values less than 1 are not allowed.")
        System.exit(1)
    } else if((cycle_channel_tuple_list.size() > 0) &&
              ( ((sample_sheet_row[0] != cycle_channel_tuple_list.last()[0]) &&
                 (sample_sheet_row[0] != cycle_channel_tuple_list.last()[0]+1)) ||
                ((sample_sheet_row[1] != cycle_channel_tuple_list.last()[1]) &&
                 (sample_sheet_row[1] != cycle_channel_tuple_list.last()[1]+1)) )) {
        throw new Exception("Error: cycle_number and channel_number must be sequential with no gaps")
        System.exit(1)
    } else {
        cycle_channel_tuple_list.add(curr_tuple)
    }
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
