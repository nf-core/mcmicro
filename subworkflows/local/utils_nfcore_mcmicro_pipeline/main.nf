//
// Subworkflow with functionality specific to the nf-core/mcmicro pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

import groovy.io.FileType
import groovy.json.JsonSlurper

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input_cycle       //  string: Path to input_cycle samplesheet
    input_sample      //  string: Path to input_sample samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )
    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    def input_type
    if (input_sample) {
        input_type = "sample"
        ch_samplesheet = Channel.fromSamplesheet(
            "input_sample",
            skip_duplicate_check: false
        )
        .tap { ch_raw_samplesheet }
        .map { validateInputSamplesheetRow(it, input_type) }
        .map { make_ashlar_input_sample(it) }

    } else if (input_cycle) {
        input_type = "cycle"
        ch_samplesheet = Channel.fromSamplesheet(
            "input_cycle",
            skip_duplicate_check: false
        )
        .tap { ch_raw_samplesheet }
        .map { validateInputSamplesheetRow(it, input_type) }
        .map { [[id:it[0]], it[3]] }
        .groupTuple()

    } else {
        error "Either input_sample or input_cycle is required."
    }

    Channel.fromSamplesheet(
        "marker_sheet",
        skip_duplicate_check: false
        )
        .tap { markersheet_data }
        .toList()
        .map{ validateInputMarkersheet(it) }

    ch_raw_samplesheet
        .toList()
        .concat(markersheet_data)
        .toList()
        .map { validateInputSamplesheetMarkersheet(it, input_type) }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {

    if (params.input_sample && params.input_cycle) {
        error "You must specify EITHER input_sample OR input_cycle, but not both."
    } else if(!params.input_sample && !params.input_cycle) {
        error "You must specify either input_sample or input_cycle."
    }
}

//
// Validate channels from input samplesheet
//

def validateInputMarkersheet( sheet_data ) {

    def idx_marker_name = input_sheet_index("marker", "marker_name")
    def idx_channel_number = input_sheet_index("marker", "channel_number")
    def idx_cycle_number = input_sheet_index("marker", "cycle_number")
    def marker_name_list = []
    def channel_number_list = []
    def cycle_number_list = []

    sheet_data.each { curr_list ->
        def idx = 0
        curr_list.each { curr_val ->
            curr_pair = [-1, -1]
            if (idx == idx_marker_name) {
                if (marker_name_list.contains(curr_val)) {
                    error("Error: duplicate marker name found in marker sheet!")
                } else {
                    marker_name_list.add(curr_val)
                }
            } else if (idx == idx_channel_number) {
                if (curr_val <= 0) {
                    error("Error: channel_number must be >= 1")
                } else if (channel_number_list && (curr_val != channel_number_list[-1] && curr_val != channel_number_list[-1].toInteger() + 1)) {
                    error("Error: channel_number cannot skip values and must be in order!")
                } else {
                    channel_number_list.add(curr_val)
                }
            } else if (idx == idx_cycle_number) {
                if (curr_val <= 0) {
                    error("Error: cycle_number must be >= 1")
                } else if (cycle_number_list && (curr_val != cycle_number_list[-1] && curr_val != cycle_number_list[-1].toInteger() + 1)) {
                    error("Error: cycle_number cannot skip values and must be in order!")
                } else {
                    cycle_number_list.add(curr_val)
                }
            }
            idx++
        }
    }

    // uniqueness of (channel, cycle) tuple in marker sheet
    def test_tuples = [channel_number_list, cycle_number_list].transpose()
    if ( test_tuples.size() != test_tuples.unique( false ).size() ) {
        error("Error: duplicate (channel,cycle) pair")
    }

    return sheet_data
}

// function that returns the index of a given column from a given sheet
//   as defined in the schema file.
//   (will need to be updated when the schema files change)
def input_sheet_index( sheet_type, column_name ) {

    def index_map = [:]
    if (sheet_type == "sample") {
        index_map = [sample: 0, image_directory: 1, cycle_images: 2, dfp: 3, ffp: 4]
    } else if (sheet_type == "cycle") {
        index_map = [sample: 0, cycle_number: 1, channel_count: 2, image_tiles: 3, dfp: 4, ffp: 5]
    } else if (sheet_type == "marker") {
        index_map = [channel_number: 0, cycle_number: 1, marker_name: 2, filter: 3,
                        excitation_wavelength: 4, emission_wavelength: 5]
    } else {
        error("Error: bad sheet type: $sheet_type")
    }

    return index_map[column_name]
}

def validateInputSamplesheetRow ( row, mode ) {
    // TODO: Add sample sheet validation for cycle inputs

    if (mode == "sample") {
        // check for the existence of all files under cycle_image column in the given image_directory
        if (row.size() >= 3 && row[2] != []) {
            def file_list = row[2].split(" ")
            file_list.each { curr_file ->
                def curr_path = new File(row[1].toString() + "/" + curr_file)
                if (!curr_path.exists()) {
                    error("Error: file in samplesheet not found: $curr_path")
                }
            }
        }
    }

    return row
}

def validateInputSamplesheetMarkersheet ( sheet_data, mode ) {
    if (mode == 'cycle' ) {
        def ctr = 0
        def sample_cycle_list = []
        def marker_cycle_list = []
        def idx_sample_cycle = input_sheet_index("cycle", "cycle_number")
        def idx_marker_cycle = input_sheet_index("marker", "cycle_number")

        sheet_data.each { curr_list ->
            if (ctr == 0) {
                curr_list.each { curr_sublist ->
                    sample_cycle_list.add(curr_sublist[idx_sample_cycle])
                }
            } else {
                marker_cycle_list.add(curr_list[idx_marker_cycle])
            }
            ctr++
        }
        if (marker_cycle_list.unique() != sample_cycle_list.unique() ) {
            error("Error: cycle_number in sample and marker sheets must match 1:1!")
        }
    } else if ( mode == 'sample' ) {
        // TODO: add validation for 1 row per sample samplesheet & markersheet correspondence
    } else {
        error("Error: bad mode $mode in validateInputSamplesheetMarkersheet()")
    }
}

def make_ashlar_input_sample( samplesheet_row ) {

    def cycle_images = []
    def index_sample_cycle_images = input_sheet_index("sample", "cycle_images")
    def index_sample_image_directory = input_sheet_index("sample","image_directory")
    def index_sample_sample = input_sheet_index("sample", "sample")

    if (samplesheet_row[index_sample_cycle_images]) {
        def tmp_path = samplesheet_row[index_sample_image_directory]
        if (tmp_path[-1] != "/") {
            tmp_path = "${tmp_path}/"
        }
        cycle_images = samplesheet_row[index_sample_cycle_images].split(' ').collect{ "${tmp_path}${it}" }
        cycle_images.each{ file_path ->
            def file_test = new File(file_path)
            if (!file_test.exists()) {
                Nextflow.error("Error: ${file_path} does not exist!")
            }
        }
    } else {
        // TODO: remove this option or allow it to grab all files when no column in the samplesheet?
        def image_dir = samplesheet_row[index_sample_image_directory]
        image_dir.eachFileRecurse (FileType.FILES) {
            if(it.toString().endsWith(".ome.tif")){
                cycle_images << file(it)
            }
        }
    }

    return [[id:samplesheet_row[index_sample_sample]], cycle_images]
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText( mqc_methods_yaml ) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
