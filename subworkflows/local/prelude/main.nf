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
    output_file_xml         = ( xml.map {
                                    meta, xmlpath -> [meta, xmlpath]
                                    } | SUMMARY_XML ).output

    output_file_markersheet = ( markersheet.map {
                                    [["id": "markers"], it]
                                    } | SUMMARY_MARKERSHEET_LITERAL ).output

    output_file_samplesheet = ( samplesheet.map {
                                    meta, image_tiles, dfp, ffp ->
                                    [meta, [image_tiles, dfp, ffp]]
                                    } | SUMMARY_SAMPLESHEET).output
}
