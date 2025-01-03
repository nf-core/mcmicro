//
// Subworkflow with functionality specific to the nf-core/mcmicro pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
import groovy.io.FileType

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input_cycle       //  string: Path to input_cycle samplesheet
    input_sample      //  string: Path to input_sample samplesheet
    marker_sheet      //  string: Path to marker_sheet

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
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
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
    // Create channel from input file provided through params.input_cycle or .input_sample
    //
    if (input_cycle) {
        // TODO: Validate that cycle_number is 1..N, in order, for all samples.
        ch_samplesheet = Channel.fromList(samplesheetToList(params.input_cycle, "${projectDir}/assets/schema_input_cycle.json"))
            .map{
                sample, cycle_number, channel_count, image_tiles, dfp, ffp ->
                [
                    [id: sample, cycle_number: cycle_number, channel_count: channel_count],
                    image_tiles,
                    dfp,
                    ffp
                ]
            }
            .dump(tag: 'ch_samplesheet (cycle)')
    } else if (input_sample) {
        ch_samplesheet = Channel.fromList(samplesheetToList(params.input_sample, "${projectDir}/assets/schema_input_sample.json"))
            .flatMap { expandSampleRow(it) }
            .dump(tag: 'ch_samplesheet (sample)')
    }

    ch_markersheet = Channel.fromList(samplesheetToList(params.marker_sheet, "${projectDir}/assets/schema_marker.json"))
        .toList()
        .map{ validateInputMarkersheet(it, "${projectDir}/assets/schema_marker.json", params) }
        .dump(tag: 'ch_markersheet')

    ch_samplesheet.toList()
        .concat(ch_markersheet)
        .toList()
        .map{ samples, markers -> validateInputSamplesheetMarkersheet(samples, markers) }

    emit:
    samplesheet = ch_samplesheet
    markersheet = ch_markersheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_report.toList()
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

    if (params.cellpose_model && !segmentation_list.contains('cellpose')) {
        error "You can only provide a cellpose model if you have selected cellpose as one of your segmentation methods"
    }

    if (params.tma_dearray && !params.pixel_size) {
        error "You must also provide the pixel_size parameter (image pixel width in microns) when enabling tma_dearray."
    }
}

//
// Validate channels from input samplesheet
//
def validateInputMarkersheet(markersheet_data, schema_file, params) {
    def schema = new groovy.json.JsonSlurper().parse(new File(schema_file))

    def required_columns = schema.items.required ?: []
    def properties = schema.items.properties

    // Create a map of column names to their indices based on the order in the schema
    def column_indices = properties.keySet().withIndex().collectEntries { column, index ->
        [(column): index]
    }

    def marker_name_list = []
    def channel_number_list = []
    def cycle_number_list = []

    // Collect all marker data
    def marker_data = markersheet_data.collect { row ->
        [
            marker_name: row[column_indices['marker_name']],
            channel_number: row[column_indices['channel_number']] as int,
            cycle_number: row[column_indices['cycle_number']] as int,
            background: row[column_indices['background']] != [] ? row[column_indices['background']] : null,
            exposure: row[column_indices['exposure']] != [] ? row[column_indices['exposure']] : null,
            remove: column_indices.containsKey('remove') ? (row[column_indices['remove']] != [] ? row[column_indices['remove']] : null) : null
        ]
    }

    // Validate basic requirements and collect lists
    marker_data.each { row ->
        // Check if all required columns are present and non-empty
        required_columns.each { column ->
            if (row[column] == null || row[column].toString().trim().isEmpty()) {
                error("Missing required value in column '${column}' for marker '${row.marker_name}'.")
            }
        }

        if (marker_name_list.contains(row.marker_name)) {
            error("Duplicate marker name found: '${row.marker_name}'")
        } else {
            marker_name_list.add(row.marker_name)
        }

        if (channel_number_list && (row.channel_number != channel_number_list[-1] && row.channel_number != channel_number_list[-1] + 1)) {
            error("Channel_number cannot skip values and must be in order! Error at marker '${row.marker_name}'")
        } else {
            channel_number_list.add(row.channel_number)
        }

        if (cycle_number_list && (row.cycle_number != cycle_number_list[-1] && row.cycle_number != cycle_number_list[-1] + 1)) {
            error("Cycle_number cannot skip values and must be in order! Error at marker '${row.marker_name}'")
        } else {
            cycle_number_list.add(row.cycle_number)
        }
    }

    // uniqueness of (channel, cycle) tuple in marker sheet
    def test_tuples = [channel_number_list, cycle_number_list].transpose()
    def dups = test_tuples.countBy{ it }.findAll{ _, count -> count > 1 }*.key
    if (dups) {
        error("Duplicate [channel, cycle] pairs: ${dups}")
    }

    // Validate backsub data
    // Check if backsub columns are present
    def backsub_columns = ['exposure', 'background', 'remove']
    def has_backsub_columns = markersheet_data.any { row ->
        backsub_columns.any { column ->
            column_indices.containsKey(column) && row[column_indices[column]] && row[column_indices[column]] != []
        }
    }

    // Throw error if backsub = false but backsub columns are present
    if (!params.backsub && has_backsub_columns) {
        error("Error: exposure, background, or remove columns are present in the marker sheet, but params.backsub is set to false. Either remove these columns or set params.backsub to true.")
    }

    def has_any_background = marker_data.any { it.background != null }

    if (params.backsub && !has_any_background) {
        error("Backsub is enabled, but all values in the background column are empty. No subtraction will occur. Either set params.backsub=false or specify how the channel subtraction should be performed.")
    } else if (has_any_background) {
        // Create a set of all markers used as background
        def markers_used_as_background = marker_data.findAll { it.background != null }.collect { it.background }.toSet()

        marker_data.each { row ->
            if (row.background != null) {
                if (row.exposure == null) {
                    error("Missing exposure value for marker '${row.marker_name}' with background '${row.background}'.")
                }
                if (!marker_name_list.contains(row.background)) {
                    error("Background value '${row.background}' specified for marker '${row.marker_name}' does not exist in the marker_name column.")
                }
            }

            // Check if this marker is used as a background and ensure it has an exposure value
            if (markers_used_as_background.contains(row.marker_name) && row.exposure == null) {
                error("Marker '${row.marker_name}' is used as a background for another marker but does not have an exposure value.")
            }
        }
    }

    return markersheet_data
}

def validateInputSamplesheetMarkersheet ( samples, markers ) {
    def sample_cycles = samples.collect{ meta, image_tiles, dfp, ffp -> meta.cycle_number }
    def marker_cycles = markers.collect{ channel_number, cycle_number, marker_name, _1, _2, _3, _4, _5, _6 -> cycle_number }

    if (marker_cycles.unique(false) != sample_cycles.unique(false) ) {
        error("cycle_number values must match between sample and marker sheets")
    }
}

def expandSampleRow( row ) {
    def (sample, image_directory, dfp, ffp) = row
    def files = []

    file(image_directory).eachFileRecurse (FileType.FILES) {
        if(it.toString().endsWith(".ome.tif")){
            files << file(it)
        }
    }

    return files.withIndex(1).collect{ f, i ->
        [[id: sample, cycle_number: i], f, dfp, ffp]
    }
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
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

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

