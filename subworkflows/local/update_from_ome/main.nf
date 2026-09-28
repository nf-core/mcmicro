include { OMEVALIDATION   } from '../../../modules/local/omevalidation/main'

workflow UPDATE_FROM_OME {
    take:
    samplesheet
    markersheet
    xml

    main:

    if (workflow.stubRun) {
        return
    }

    def val_data = xml.dump(tag: "XMLpreOME").map{meta, x -> [meta, x]} | OMEVALIDATION

    samplesheet.join(val_data)
         .map { original_meta, image_tiles, dfp, ffp, xml_meta, marker_data ->
                [xml_meta + original_meta, image_tiles, dfp, ffp]
        }
        .dump(tag: "ch_samplesheet_meta")
        .set { samplesheet_meta }

    def c_sum = 0
    def agg   = samplesheet_meta
                    .map {
                        meta, image_tiles, dfp, ffp -> [meta.cycle_number, meta.channel_count]
                    }
                    .unique()
                    .toSortedList()
                    .flatMap()
                    .map {
                        cn, cc ->
                        def temp = [cn, c_sum]
                        c_sum += cc
                        return temp
                    }
                    .dump(tag: "CHANNEL_DELTA_INDEX")

    // update val_data channel_number so it matches with samplesheet
    val_data.map {
            meta, xml_meta, marker_meta -> marker_meta.collect{ it -> meta.subMap('id') + it }
        }
        .flatten()
        .map{
            meta -> [meta.cycle_number, meta]
        }
        .combine(agg, by: 0).dump(tag:"COMBINE")
        .map {
            key, meta, counter ->
            def temp = meta + [channel_number: meta.channel_number + counter]
            return temp
            //meta.channel_number += counter
            //return meta
        }
        .dump(tag:'val_data_markers')
        .set{ val_data_markers }

    def markersheet_template =
    [
        'channel_number':null, 'cycle_number':null,
        'marker_name':null, 'filter':null,
        'excitation_wavelength':null, 'emission_wavelength':null,
        'exposure':null, 'background':null,
        'remove':null
        ]

    markersheet
        .flatten()
        .combine(samplesheet_meta.map{ it[0].subMap('id') }.unique())
        .map{ it -> it[0] + it[1] }
        .dump(tag: "ch_markersheet_premeta")
        .map{ e -> [e.subMap('id', 'channel_number', 'cycle_number'), e.findAll{ it.value != null }] }
        .join(
            val_data_markers.map{ x -> [x.subMap('id', 'channel_number', 'cycle_number'), x] },
            remainder: true
        )
        .map{ it -> it.drop(1) }
        .dump(tag:'ch_markersheet_mismatch_check')
        .filter{ e -> !e.any{it == null}} // ignore errors, will be caught in PRELUDE
        .map{ orig, validated ->  markersheet_template + validated + orig }
        .map{ meta -> meta - meta.subMap('id') }
        .unique()
        .toSortedList { a, b -> a.channel_number <=> b.channel_number }
        .dump(tag: "ch_markersheet_meta")
        .set { markersheet_meta }



    if (params.backsub) {
        markersheet_meta
            .map { entry ->
                if (entry.exposure_time == null){
                    error "Exposure time cannot be NULL if doing backsub"
                }
            }
    }

    emit:
    samplesheet = samplesheet_meta
    markersheet = markersheet_meta
}
