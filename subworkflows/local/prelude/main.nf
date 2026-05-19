include { PRELUDE_SUMMARY_XML          } from '../../../modules/local/summary_xml/main'
include { PRELUDE_SUMMARY_MARKERSHEET  } from '../../../modules/local/summary_markersheet/main'
include { PRELUDE_SUMMARY_SAMPLESHEET  } from '../../../modules/local/summary_samplesheet/main'
include { PRELUDE_MULTI_SUMMARY        } from '../../../modules/local/multi_xml_summary/main'
include { PRELUDE_MULTI_MATRIX_SUMMARY } from '../../../modules/local/multi_matrix_summary/main'

workflow PRELUDE {
    take:
    markersheet
    samplesheet
    xml

    main:
    ch_output_xml         = xml.map {
                                    meta, xmlpath -> [meta, xmlpath]
                                    } | PRELUDE_SUMMARY_XML

    ch_output_markersheet = ( markersheet.map {
                                    it ->
                                    def f = file("$it".md5() + ".json")
                                    f.write(new groovy.json.JsonBuilder(it).toString())
                                    return [["id": "markers"], f]
                                    } | PRELUDE_SUMMARY_MARKERSHEET ).output

    ch_output_samplesheet = ( samplesheet.map {
                                    meta, image_tiles, dfp, ffp ->
                                    def f = file("$meta".md5() + ".json")
                                    f.write(new groovy.json.JsonBuilder(meta).toString())
                                    return [meta, f]
                                    } | PRELUDE_SUMMARY_SAMPLESHEET ).output

    ch_output_merged_xml  = PRELUDE_MULTI_SUMMARY([], ch_output_xml.variables.map{
                                                        meta, file -> file
                                                        }
                                                        .collect()
                                                        .dump(tag:"ch_output_xml.variables"))

    ch_output_mixed_matrix_summary  = PRELUDE_MULTI_MATRIX_SUMMARY(
                                ch_output_samplesheet.map{
                                    meta, files -> files
                                    }.collect().dump(tag: "prelude_samplesheet"),
                                ch_output_xml.output.map{
                                    meta, files -> files
                                    }.collect().dump(tag: "prelude_xml"),
                                ch_output_merged_xml.output_mqc.collect().dump(tag: "prelude_multi_xml")
                                )


    emit:
    output_file_xml                  = ch_output_xml.output.map{ meta, files -> files}
    output_file_samplesheet          = ch_output_samplesheet.map{ meta, files -> files}
    output_file_markersheet          = ch_output_markersheet.map{ meta, files -> files}
    output_file_merged_xml           = ch_output_merged_xml.output_mqc
    output_file_error_merged         = ch_output_merged_xml.output_errors
    output_file_mixed_matrix_summary = ch_output_mixed_matrix_summary.output
}
