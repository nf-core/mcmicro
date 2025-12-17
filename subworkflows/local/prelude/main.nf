include { SUMMARY_XML                   } from '../../../modules/local/prelude/main'
include { SUMMARY_MARKERSHEET_LITERAL   } from '../../../modules/local/prelude/main'
include { SUMMARY_SAMPLESHEET           } from '../../../modules/local/prelude/main'

workflow PRELUDE {
    take:
    markersheet
    samplesheet
    xml

    main:
    ch_output_xml         = ( xml.map {
                                    meta, xmlpath -> [meta, xmlpath]
                                    } | SUMMARY_XML ).output

    ch_output_markersheet = ( markersheet.map {
                                    [["id": "markers"], it]
                                    } | SUMMARY_MARKERSHEET_LITERAL ).output

    ch_output_samplesheet = ( samplesheet.map {
                                    meta, image_tiles, dfp, ffp ->
                                    [meta, [image_tiles, dfp, ffp]]
                                    } | SUMMARY_SAMPLESHEET).output

    emit:
    output_file_xml         = ch_output_xml
    output_file_samplesheet = ch_output_samplesheet
    output_file_markersheet = ch_output_markersheet
}
