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

    marker_map = make_sheet_map(marker_sheet)
    check_marker_sheet(marker_map)
    make_version_file()
    if (input_type == 'cycle') {
        check_sample_and_marker_sheet(marker_map, params.input_cycle)
    }

    emit:
    // csv = MARKER_SHEET_CHECK.out.csv             // channel: [ marker_sheet.valid.csv ]
    csv = marker_sheet                              // channel: [ marker_sheet.valid.csv ]
    versions = MARKER_SHEET_CHECK.out.versions      // channel: [ versions.yml ]

}

// Functions

def make_sheet_map(path_sheet) {
    marker_map = [:]
    ctr = 0
    new File(path_sheet).splitEachLine(",") { fields ->
        if ( ctr == 0 ) {
            keys = fields.unique( false )
            if ( keys.size() != fields.size() ) {
                throw new Exception("Error: duplicate header name found in marker sheet!")
            }
            keys.each { value ->
                marker_map[value] = []
            }
        } else {
            idx = 0
            fields.each { value ->
                marker_map[keys[idx]].add(value)
                idx++
            }
        }
        ctr++
    }
    return marker_map
}

def check_sample_and_marker_sheet(marker_map, path_sample_sheet) {

    // cycle_number in marker and sample sheet must match 1:1

    sample_map = make_sheet_map(path_sample_sheet)

    if (marker_map['cycle_number'].unique(false) != sample_map['cycle_number'].unique(false)) {
        throw new Exception("Error: cycle_number in marker and sample sheets must match 1:1!")
    }
}

def make_version_file() {
}

def check_marker_sheet(marker_map) {

    /*
    marker_map = [:]
    ctr = 0
    new File(path_marker_sheet).splitEachLine(",") { fields ->
        if ( ctr == 0 ) {
            keys = fields.unique( false )
            if ( keys.size() != fields.size() ) {
                throw new Exception("Error: duplicate header name found in marker sheet!")
            }
            keys.each { value ->
                marker_map[value] = []
            }
        } else {
            idx = 0
            fields.each { value ->
                marker_map[keys[idx]].add(value)
                idx++
            }
        }
        ctr++
    }
    */

    // uniqueness of marker name in marker sheet
    if ( marker_map['marker_name'].size() != marker_map['marker_name'].unique( false ).size() ) {
        throw new Exception("Error: duplicate marker name found in marker sheet!")
    }

    // uniqueness of (channel, cycle) tuple in marker sheet
    test_tuples = [marker_map['channel_number'], marker_map['cycle_number']].transpose()
    if ( test_tuples.size() != test_tuples.unique( false ).size() ) {
        throw new Exception("Error: duplicate (channel,cycle) pair")
    }

    // cycle and channel are 1-based so 0 should throw an exception
    if (marker_map['channel_number'].stream().anyMatch { it.toInteger() < 1 }) {
        throw new Exception("Error: channel_number must be >= 1")
    }

    // cycle and channel cannot skip values and must be in order
    prev_cycle = marker_map['cycle_number'][0]
    marker_map['cycle_number'].each { curr_cycle ->
        if ( (curr_cycle.toInteger() != prev_cycle.toInteger() ) && (curr_cycle.toInteger() != prev_cycle.toInteger() + 1) ) {
            throw new Exception("Error: cycle_number cannot skip values and must be in order!")
        }
        prev_cycle = curr_cycle
    }
    prev_channel = marker_map['channel_number'][0]
    marker_map['channel_number'].each { curr_channel ->
        if ( (curr_channel.toInteger() != prev_channel.toInteger() ) && (curr_channel.toInteger() != (prev_channel.toInteger() + 1)) ) {
            throw new Exception("Error: channel_number cannot skip values and must be in order!")
        }
        prev_channel = curr_channel
    }
}
