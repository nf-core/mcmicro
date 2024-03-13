//
// Check input samplesheet and get read channels
//
/* TODO: commented out SAMPLESHEET_CHECK because it doesn't check for anything more than
         fromSamplesheet does.
         Leaving commented code in just place in case we want to do additional validation later
include { SAMPLESHEET_CHECK  } from '../../modules/local/samplesheet_check' */
include { MARKER_SHEET_CHECK  } from '../../modules/local/marker_sheet_check'
include { SAMPLE_AND_MARKER_SHEET_CHECK } from '../../modules/local/sample_and_marker_sheet_check'

parameters_schema = "assets/nextflow_schema.json"

workflow INPUT_CHECK {
    take:
    input_type          // either 'sample' or 'cycle'
    samplesheet_sample  // file: /path/to/input_sample.csv
    samplesheet_cycle   // file: /path/to/input_cycle.csv
    marker_sheet        // file: /path/to/marker_sheet.csv

    main:

    if ( input_type == "sample" ) {
        /* commented out redundant checks
        SAMPLESHEET_CHECK ( input_type, samplesheet_sample )
        input = Channel.fromSamplesheet(
            "input_sample",
            parameters_schema: parameters_schema,
            skip_duplicate_check: false)
        */
        SAMPLE_AND_MARKER_SHEET_CHECK ( params.input_sample, params.marker_sheet )
    } else if ( input_type == "cycle" ) {
        /* commented out redundant checks
        SAMPLESHEET_CHECK ( input_type, samplesheet_cycle )
        input = Channel.fromSamplesheet(
            "input_cycle",
            parameters_schema: parameters_schema,
            skip_duplicate_check: false)
        */
        SAMPLE_AND_MARKER_SHEET_CHECK ( params.input_cycle, params.marker_sheet )
    }

    MARKER_SHEET_CHECK ( params.marker_sheet )
    marker = Channel.fromSamplesheet(
        "marker_sheet",
        parameters_schema: parameters_schema,
        skip_duplicate_check: false
        )

    emit:
    csv = MARKER_SHEET_CHECK.out.csv            // channel: [ marker_sheet.valid.csv ]
    versions = MARKER_SHEET_CHECK.out.versions  // channel: [ versions.yml ]
}
