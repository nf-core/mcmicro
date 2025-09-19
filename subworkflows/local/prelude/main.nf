include { SUMMARY_XML                   } from '../../../modules/local/prelude/main'
include { SUMMARY_MARKERSHEET_LITERAL   } from '../../../modules/local/prelude/main'
include { SUMMARY_SAMPLESHEET           } from '../../../modules/local/prelude/main'

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
    output_file_xml         = ( xml.map{meta, xmlpath -> [meta, xmlpath]} | SUMMARY_XML ).output
    output_file_markersheet = SUMMARY_MARKERSHEET_LITERAL( markersheet.map{[["id": "markers"], it]} ).output
    output_file_samplesheet = (samplesheet.map{meta, image_tiles, dfp, ffp -> [meta, samplesheet] } | SUMMARY_SAMPLESHEET).output
}
