//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:

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
    
    SAMPLESHEET_CHECK ( samplesheet )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_fastq_channel(it) }
            .set { input }

    // need to emit all sample sheet contents here
    //    you can split them up in mcmicro.nf

    // input.view { "all $it" }
    //input.ashlar.view { "ashlar $it" }
    //input.foo.view { "foo $it" }

    emit:
    input                                    // channel: [ val(meta), [ image ], [ marker ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ cycle_number, channel_count, image_tiles ] ]
def create_fastq_channel(LinkedHashMap row) {
    print("*** create_fastq_channel: entering... ***")

    // create meta map
    def meta = [:]
    meta.id         = row.sample + '_cycle_' + row.cycle_number
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
