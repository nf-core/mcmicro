include { SUMMARY_XML           } from '../../../modules/local/prelude/main'
include { SUMMARY_MARKERSHEET_LITERAL   } from '../../../modules/local/prelude/main'
include { SUMMARY_SAMPLESHEET   } from '../../../modules/local/prelude/main'

workflow PRELUDE {
    take:
    markersheet
    samplesheet
    xml

    emit:
    output_file_xml
    output_file_samplesheet
    output_file_markersheet

    main:
    samplesheet.map{meta, image_tiles, dfp, ffp -> meta}.set{meta}

    output_file_xml         = SUMMARY_XML(meta, xml).output
    output_file_markersheet = SUMMARY_MARKERSHEET_LITERAL(meta, markersheet).output
    output_file_samplesheet = SUMMARY_SAMPLESHEET(meta, samplesheet).output
}
