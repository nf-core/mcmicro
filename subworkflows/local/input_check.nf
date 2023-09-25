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

    emit:
    input                                    // channel: [ val(meta), [ image ], [ marker ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    // meta.tiff = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    //def fastq_meta = []
    if (!file(row.image).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> image file does not exist!\n${row.image}"
    }
    if (!file(row.marker).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> image file does not exist!\n${row.marker}"
    }

    // if (meta.single_end) {
    //     fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    // } else {
    //     if (!file(row.fastq_2).exists()) {
    //         exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    //     }
    //}
    image_meta = [ meta, [ file(row.image) ], [file(row.marker)] ]

    return image_meta
}
