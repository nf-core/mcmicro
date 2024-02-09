//
// Check input samplesheet and get read channels
//
include { SAMPLESHEET_CHECK  } from '../../modules/local/samplesheet_check'
include { MARKER_SHEET_CHECK  } from '../../modules/local/marker_sheet_check'
include { SAMPLE_AND_MARKER_SHEET_CHECK } from '../../modules/local/sample_and_marker_sheet_check'

// TODO: fails nf-test with relative path; I don't like having this full path buried in a subworkflow; can this be improved?
parameters_schema = '/home/pollen/github/mcmicro-nf-core/nextflow_schema.json'

workflow INPUT_CHECK {
    take:
    input_type          // either 'sample' or 'cycle'
    samplesheet_sample  // file: /path/to/input_sample.csv
    samplesheet_cycle   // file: /path/to/input_cycle.csv
    marker_sheet        // file: /path/to/marker_sheet.csv

    main:

    // not running this check because fromSamplesheet checks columns and data format
    //   we can add it back if we'd like to do more in depth validation
    // currently just running to output version.yml
    SAMPLESHEET_CHECK ( input_type )

    if ( input_type == "sample" ) {
        input = Channel.fromSamplesheet(
            "input_sample",
            parameters_schema: parameters_schema,
            skip_duplicate_check: false)
        SAMPLE_AND_MARKER_SHEET_CHECK ( params.input_sample, params.marker_sheet )
    } else if ( input_type == "cycle" ) {
        input = Channel.fromSamplesheet(
            "input_cycle",
            parameters_schema: parameters_schema,
            skip_duplicate_check: false)
        SAMPLE_AND_MARKER_SHEET_CHECK ( params.input_cycle, params.marker_sheet )
    }

    MARKER_SHEET_CHECK ( params.marker_sheet )
    marker = Channel.fromSamplesheet(
        "marker_sheet",
        parameters_schema: parameters_schema,
        skip_duplicate_check: false
        )

    /*
    if( params.illumination){
        SAMPLESHEET_CHECK ( samplesheet )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_fastq_channel(it) }
            .groupTuple(by: [0])
            .map { meta, reads -> [ meta, reads.flatten() ] }
            .set { images_merged }
    }
    */

    emit:
    input                                     // channel: [ val(meta), [ image ], [ marker ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ cycle_number, channel_count, image_tiles ] ]
def create_fastq_channel(LinkedHashMap row) {
    print("*** create_fastq_channel: entering... ***")

    // create meta map
    def meta = [:]
    meta.id         = row.sample
    // meta.tiff = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    //def fastq_meta = []

    '''
    if (!file(row.image).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> image file does not exist!\n${row.image}"
    }
    if (!file(row.marker).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> image file does not exist!\n${row.marker}"
    }
    '''

    // if (meta.single_end) {
    //     fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    // } else {
    //     if (!file(row.fastq_2).exists()) {
    //         exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    //     }
    //}
    //image_meta = [ meta, [ file(row.image) ], [file(row.marker)] ]
    image_meta = [ meta, [ row.cycle_number, row.channel_count, file(row.image_tiles) ] ]

    return image_meta
}
