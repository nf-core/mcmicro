include { BFTOOLS_SHOWINF } from '../../../modules/nf-core/bftools/showinf/main'
include { OMEVALIDATION   } from '../../../modules/local/omevalidation/main'

workflow UPDATE_FROM_OME {
    take:
    samplesheet
    markersheet

    main:
    samplesheet.map{meta, image_tiles, dfp, ffp -> [meta, image_tiles]} | BFTOOLS_SHOWINF

    val_data = BFTOOLS_SHOWINF.out.xml | OMEVALIDATION

    samplesheet.join(val_data)
         .map { original_meta, image_tiles, dfp, ffp, xml_meta, marker_data ->
                [xml_meta + original_meta, image_tiles, dfp, ffp]
        }
        .dump(tag: "ch_samplesheet_meta")
        .set { samplesheet_meta }

    c_sum = 0
    agg   = channel.empty()

    samplesheet_meta
        .map {
            meta, image_tiles, dfp, ffp -> [meta.cycle_number, meta.channel_count]
        }
        .unique()
        .toSortedList()
        .flatMap()
        .map {
            cn, cc ->
            temp = [cn, c_sum]
            c_sum += cc
            return temp
        }
        .dump(tag: "CHANNEL_DELTA_INDEX")
        .set { agg }

    // update val_data channel_number so it matches with samplesheet
    val_data.map {
            meta, xml_meta, marker_meta -> marker_meta
        }
        .flatten()
        .map{
            meta -> [meta.cycle_number, meta]
        }
        .combine(agg, by: 0)
        .map {
            key, meta, counter ->
            meta.channel_number += counter
            return meta
        }
        //.dump(tag:"vd")
        .set{ val_data }


    markersheet_template = //markersheet.flatten().first().keySet().collectEntries {key -> [key, null]}.toList().first().dump(tag: "template")
    [
        'channel_number':null, 'cycle_number':null,
        'marker_name':null, 'filter':null,
        'excitation_wavelength':null, 'emission_wavelength':null,
        'exposure':null, 'background':null,
        'remove':null
        ]

    markersheet
        .flatten()
        .dump(tag: "ch_markersheet_premeta")

        .map {
            entry ->
            [entry.subMap('channel_number', 'cycle_number'), entry.findAll {it.value != null}]
        }
        .join(
            val_data
            .map {
                x-> [x.subMap('channel_number', 'cycle_number'), x]
            }
        )
        .map {
            channel, orig, validated ->
            markersheet_template + validated + orig
        }
        .toSortedList {
            a, b -> a.channel_number <=> b.channel_number
        }
        .dump(tag: "ch_markersheet_meta")
        .set { markersheet_meta }


    markersheet_meta  // Inter sample marker setup check
        .map {
            entry -> [[entry.cycle_number, entry.channel_number], [entry.exposure_time, entry.exposure_time_unit]]
        }
        .groupTuple()
        .map {
            key, values ->
                if (values.unique().size() != 1)
                    error "Inconsistent marker exposure data across samples."
        }

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
    versions    = BFTOOLS_SHOWINF.out.versions
}
