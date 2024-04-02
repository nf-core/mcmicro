//
// Check input samplesheet and get read channels
//
/* TODO: commented out SAMPLESHEET_CHECK because it doesn't check for anything more than
         fromSamplesheet does.
         Leaving commented code in just place in case we want to do additional validation later
include { SAMPLESHEET_CHECK  } from '../../modules/local/samplesheet_check' */
include { MARKER_SHEET_CHECK  } from '../../modules/local/marker_sheet_check'
include { SAMPLE_AND_MARKER_SHEET_CHECK } from '../../modules/local/sample_and_marker_sheet_check'

import groovy.json.JsonSlurper

parameters_schema = "${projectDir}/assets/nextflow_schema.json"
markersheet_schema = "${projectDir}/assets/schema_marker.json"
input_cycle_schema = "${projectDir}/assets/schema_input_cycle.json"

workflow INPUT_CHECK {
    take:
    input_type          // either 'sample' or 'cycle'
    samplesheet_sample  // file: /path/to/input_sample.csv
    samplesheet_cycle   // file: /path/to/input_cycle.csv
    marker_sheet        // file: /path/to/marker_sheet.csv

    main:

    /* need needed with the new groovy validation
    if ( input_type == "sample" ) {
        SAMPLE_AND_MARKER_SHEET_CHECK ( params.input_sample, params.marker_sheet )
    } else if ( input_type == "cycle" ) {
        SAMPLE_AND_MARKER_SHEET_CHECK ( params.input_cycle, params.marker_sheet )
    }
    */

    // MARKER_SHEET_CHECK ( params.marker_sheet )

    marker = Channel.fromSamplesheet(
        "marker_sheet",
        parameters_schema: parameters_schema,
        skip_duplicate_check: false
        )
        .set { marker_sheet_data }

    // NOTE: fromSamplesheet data is in column order defined in the schema not the csv
    marker_header = sheet_keys(markersheet_schema)
    Channel.from(marker_header)
        .toList()
        .concat(marker_sheet_data)
        .toList()
        .map{ validate_marker_sheet(it) }

    marker_map_new = make_sheet_map_new(markersheet_schema, marker)
    marker_map = make_sheet_map(marker_sheet)
    //check_marker_sheet(marker_map)

    if (input_type == 'cycle') {

        sample = Channel.fromSamplesheet(
            "input_cycle",
            parameters_schema: parameters_schema,
            skip_duplicate_check: false
        )
        .set { sample_cycle_data }

        Channel.from(marker_header)
            .toList()
            .concat(marker_sheet_data)
            .toList()
            .set { marker_data }

        Channel.of(["divider"])
            .set { divider }

        sample_cycle_header = sheet_keys(input_cycle_schema)
        Channel.from(sample_cycle_header)
            .toList()
            .concat(sample_cycle_data)
            .concat(divider)
            .concat(marker_data)
            .toList()
            .map { ch_validate_sample_and_marker_sheets(it) }

        check_sample_and_marker_sheet(marker_map, params.input_cycle)
    }

    emit:
    // csv = MARKER_SHEET_CHECK.out.csv             // channel: [ marker_sheet.valid.csv ]
    csv = marker_sheet                              // channel: [ marker_sheet.valid.csv ]
    //versions = MARKER_SHEET_CHECK.out.versions      // channel: [ versions.yml ]
    versions = "\"NFCORE_MCMICRO:MCMICRO:INPUT_CHECK:MARKER_SHEET_CHECK\":\n\tpython: 3.8.3"
}

// Functions

def ch_validate_sample_and_marker_sheets ( sheet_data ) {
    // we're going to validate cycle number in both sheets are 1:1 match
    // TODO: there must be a less complicated way to code this!
    sample_cycle_list = []
    marker_cycle_list = []
    marker_mode = false
    ctr_list = 0
    idx_cycle_number = -1      // store index of cycle_number for each sheet
    sheet_data.each { curr_list_cvsams ->
        if ( marker_mode) {
            curr_list_cvsams.each { curr_sublist_cvsams ->
                if(ctr_list == 0){
                    idx_cycle_number = curr_list_cvsams.indexOf('cycle_number')
                    ctr_list++
                } else {
                    marker_cycle_list.add(curr_sublist_cvsams[1])
                }
            }
        } else {
            if (curr_list_cvsams == ["divider"]) {
                marker_mode = true
                ctr_list = 0
            } else {
                if (ctr_list == 0) {
                    idx_cycle_number = curr_list_cvsams.indexOf('cycle_number')
                    ctr_list++
                } else {
                    sample_cycle_list.add(curr_list_cvsams[idx_cycle_number])
                }
            }
        }
    }

    if (marker_cycle_list.unique() != sample_cycle_list.unique() ) {
        print("Error: cycle_number in sample and marker sheets must match 1:1!")
    }

    /*
    sample_cycle_map = [:]
    ctr_vscs = 0
    sheet_data.each { curr_list_vscs->
        if ( ctr_vscs == 0 ) {
            keys_vscs = curr_list_vscs.unique( false )
            keys_vscs.each { curr_val_vscs ->
                sample_cycle_map[curr_val_vscs] = []
            }
        } else {
            idx_vscs = 0
            curr_list_vscs.each { curr_val_vscs ->
                sample_cycle_map[keys_vscs[idx_vscs]].add(curr_val_vscs)
                idx_vscs++
            }
        }
        ctr_vscs++
    }
    */
}

def validate_marker_sheet( sheet_data ) {

    marker_map = [:]
    ctr = 0
    sheet_data.each { curr_list ->
        if ( ctr == 0 ) {
            keys = curr_list.unique( false )
            keys.each { curr_val ->
                marker_map[curr_val] = []
            }
        } else {
            idx = 0
            curr_list.each { curr_val ->
                marker_map[keys[idx]].add(curr_val)
                idx++
            }
        }
        ctr++
    }

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
    if (marker_map['cycle_number'].stream().anyMatch { it.toInteger() < 1 }) {
        throw new Exception("Error: cycle_number must be >= 1")
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

def sheet_keys(path_marker_schema) {
    def inputFile = new File(path_marker_schema)
    def InputJSON = new JsonSlurper().parseText(inputFile.text)
    //InputJSON.each{ println it }
    return InputJSON['items']['properties'].keySet()
}

def make_sheet_map_new(path_schema, marker_data) {

    def inputFile = new File("${projectDir}/assets/schema_marker.json")
    def InputJSON = new JsonSlurper().parseText(inputFile.text)
    //InputJSON.each{ println it }
    //print(InputJSON['items']['properties'].keySet())

}

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
