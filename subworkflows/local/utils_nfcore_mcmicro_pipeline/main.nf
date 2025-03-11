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
        .map{ validateInputMarkersheet(it) }
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
    def multiqc_reports = multiqc_report.toList()

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
                multiqc_reports.getVal(),
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

def validateInputMarkersheet( markersheet_data ) {

    def marker_name_list = []
    def channel_number_list = []
    def cycle_number_list = []

    markersheet_data.each { row ->
        def (channel_number, cycle_number, marker_name) = row

        if (marker_name_list.contains(marker_name)) {
            error("Duplicate marker name found in marker sheet!")
        } else {
            marker_name_list.add(marker_name)
        }

        if (channel_number_list && (channel_number != channel_number_list[-1] && channel_number != channel_number_list[-1] + 1)) {
            error("Channel_number cannot skip values and must be in order!")
        } else {
            channel_number_list.add(channel_number)
        }

        if (cycle_number_list && (cycle_number != cycle_number_list[-1] && cycle_number != cycle_number_list[-1] + 1)) {
            error("Cycle_number cannot skip values and must be in order!")
        } else {
            cycle_number_list.add(cycle_number)
        }
    }


    if (params.segmentation_recyze) {
            // Extract segmentation_channel values
            def segmentation_channel_list = markersheet_data.findResults{ _1, _2, _3, _4, _5, _6, _7, _8, _9, segmentation_channels -> segmentation_channels ?: null }
            // Check if the 'segmentation_channel' column exists
            if (segmentation_channel_list.isEmpty()) {
                error "The 'segmentation_channel' column is missing from the markersheet or all values are null. This column is required when params.segmentation_recyze is given."
            }
            // Check if at least one value in the 'segmentation_channel' column is TRUE
            def hasSegmentationChannel = segmentation_channel_list.any { it?.toString().toLowerCase() == 'true' }
            if (!hasSegmentationChannel) {
                error "The 'segmentation_channel' column must have at least one TRUE value when params.segmentation_recyze is given."
            }
        }
    // uniqueness of (channel, cycle) tuple in marker sheet
    def test_tuples = [channel_number_list, cycle_number_list].transpose()
    def dups = test_tuples.countBy{ it }.findAll{ _, count -> count > 1 }*.key
    if (dups) {
        error("Duplicate [channel, cycle] pairs: ${dups}")
    }

    // validate backsub columns if present
    def exposure_list = markersheet_data.findResults{ _1, _2, _3, _4, _5, _6, exposure, _8, _9, _10 -> exposure ?: null }
    def background_list = markersheet_data.findResults{ _1, _2, _3, _4, _5, _6, _7, background, _9, _10 -> background ?: null }
    def remove_list = markersheet_data.findResults{ _1, _2, _3, _4, _5, _6, _7, _8, remove, _10 -> remove ?: null }

    if (!background_list && (exposure_list || remove_list)) {
        error("No values in background column, but values in either exposure or remove columns.  Must have background column values to perform background subtraction.")
    } else if (background_list) {
        inter_list = marker_name_list.intersect(background_list)
        if (inter_list.size() != background_list.size()) {
            outliers_list = background_list - inter_list
            error('background column values must exist in the marker_name column. The following background column values do not exist in the marker_name column: ' + outliers_list)
        }

        if (!exposure_list) {
            error('You must have at least one value in the exposure column to perform background subtraction')
        }

        if (!remove_list) {
            error ('You must have at least one value in the remove column to perform background subtraction')
        }
    }

    return markersheet_data
}

def validateInputSamplesheetMarkersheet ( samples, markers ) {
    def sample_cycles = samples.collect{ meta, image_tiles, dfp, ffp -> meta.cycle_number }
    def marker_cycles = markers.collect{ channel_number, cycle_number, marker_name, _1, _2, _3, _4, _5, _6, _7 -> cycle_number }

    if (marker_cycles.unique(false) != sample_cycles.unique(false) ) {
        error("cycle_number values must match between sample and marker sheets")
    }

    // TODO: should the following test be in a separate validateInputSamplesheet() function?

    def channel_cycle_map = samples.collect{ meta, image_tiles, dfp, ffp -> [meta.id,meta.cycle_number] }.groupBy{ it[0] }
    channel_cycle_map.each { entry ->
        last_val = -1
        entry.value.collect{ it[1] }.each{ curr_val ->
            if (last_val != -1 && (curr_val > (last_val + 1) || curr_val <= last_val)) {
                error("cycle_number values must be increasing with no gaps")
            }
            last_val = curr_val
        }
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
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            params["illumination"] ? "Basicpy (Peng et al. 2017)," : "",
            "Ashlar (Muhlich et al. 2022),",
            params["segmentation"].contains("cellpose") ? "Cellpose (Stringer et al. 2021)," : "",
            params["segmentation"].contains("mesmer") ? "Mesmer (Van Valen et al. 2016)," : "",
            "MCQuant (Schapiro et al. 2022),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            params["illumination"] ? "<li>Peng, T., Thorn, K., Schroeder, T., Wang, L., Theis, F.J., Marr*, C., Navab*, N. (2017). A BaSiC Tool for Background and Shading Correction of Optical Microscopy Images Nature Communication 8(14836). doi: 10.1038/ncomms14836</li>" : "",
            "<li>Muhlich, J.L., Chen, Y., Yapp, C., Russell, D., Santagata, S., Sorger, P.K. (2022) Stitching and registering highly multiplexed whole-slide images of tissues and tumors using ASHLAR, Bioinformatics 38(19), 4613–4621. doi: 10.1093/bioinformatics/btac544</li>",
            params["segmentation"].contains("cellpose") ? "<li>Stringer, C., Wang, T., Michaelos, M., & Pachitariu, M. (2021). Cellpose: a generalist algorithm for cellular segmentation. Nature methods, 18(1), 100-106.</li>" : "",
            params["segmentation"].contains("mesmer") ? "<li>Van Valen, D.A., Kudo, T., Lane, K.M., Macklin, D.N., Quach, N.T., DeFelice, M.M., Maayan, I., Tanouchi, Y., Ashley, E.A., Covert, M.W. (2016). Deep Learning Automates the Quantitative Analysis of Individual Cells in Live-Cell Imaging Experiments. PLOS Computational Biology 12(11), doi: 10.1371/journal.pcbi.1005177.</li>" : "",
            "<li>Schapiro, D., Sokolov, A., Yapp, C. et al. MCMICRO: a scalable, modular image-processing pipeline for multiplexed tissue imaging. Nat Methods 19, 311–315 (2022). doi: 10.1038/s41592-021-01308-y</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
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

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
