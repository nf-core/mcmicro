//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { images }

    emit:
    images                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.tiff = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.tiff).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> tiff file does not exist!\n${row.tiff}"
    }
    // if (meta.single_end) {
    //     fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    // } else {
    //     if (!file(row.fastq_2).exists()) {
    //         exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    //     }
    //}
    fastq_meta = [ meta, [ file(row.tiff) ] ]

    return fastq_meta
}
