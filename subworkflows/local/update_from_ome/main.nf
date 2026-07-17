include { OMEVALIDATION   } from '../../../modules/local/omevalidation/main'

workflow UPDATE_FROM_OME {
    take:
    samplesheet
    markersheet
    xml

    main:
    val_data = xml.map{meta, x -> [meta, x]} | OMEVALIDATION

    samplesheet.join(val_data)
         .map { original_meta, image_tiles, dfp, ffp, xml_meta, marker_data ->
                [xml_meta + original_meta, image_tiles, dfp, ffp]
        }
        .dump(tag: "ch_samplesheet_meta")
        .set { samplesheet_meta }

    agg = channel.empty()

    samplesheet_meta
        .map {
            meta, image_tiles, dfp, ffp -> [meta.cycle_number, meta.channel_count]
        }
        .unique()
        .toSortedList()
        .map { pairs ->
            // Samples must agree on the channel count of a given cycle, otherwise a single
            // shared markers.csv cannot describe every image and MCQUANT will later fail with
            // an opaque "number of channels doesn't match" error. Catch it here instead.
            pairs.groupBy { it[0] }.each { cycle, entries ->
                def counts = entries.collect { it[1] }.unique()
                if (counts.size() != 1) {
                    error "Samples disagree on the channel count for cycle ${cycle} (found ${counts}). " +
                          "A shared markers.csv cannot be built across images with different channel layouts."
                }
            }
            // Cumulative per-cycle offset so each cycle's channels are numbered contiguously
            // across the stacked image. Computed inside a single map so it does not rely on
            // downstream operators preserving order.
            def c_sum = 0
            pairs.collect { cn, cc -> def temp = [cn, c_sum]; c_sum += cc; temp }
        }
        .flatMap()
        .dump(tag: "CHANNEL_DELTA_INDEX")
        .set { agg }

    // update val_data channel_number so it matches with samplesheet
    val_data.map {
            meta, xml_meta, marker_meta -> marker_meta.collect{ meta.subMap('id') + it }
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
        .dump(tag:'val_data_markers')
        .set{ val_data_markers }

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
        .combine(samplesheet_meta.map{ it[0].subMap('id') }.unique())
        .map{ it[0] + it[1] }
        .dump(tag: "ch_markersheet_premeta")
        .map{ e -> [e.subMap('id', 'channel_number', 'cycle_number'), e.findAll{ it.value != null }] }
        .join(
            val_data_markers.map{ x -> [x.subMap('id', 'channel_number', 'cycle_number'), x] },
            remainder: true
        )
        .map{ it.drop(1) }
        .dump(tag:'ch_markersheet_mismatch_check')
        .map{ e ->
            if (e.any{ it == null }) {
                error('Markersheet cycle/channel numbering does not match image file metadata')
            }
            e
        }
        .map{ orig, validated ->  markersheet_template + validated + orig }
        .map{ meta -> meta - meta.subMap('id') }
        // Collapse the per-sample rows down to one shared row per channel. Deduplicating on
        // the whole map (the previous `.unique()`) leaves duplicates whenever samples differ
        // in any per-sample field (e.g. exposure time), inflating markers.csv beyond the image
        // channel count (nf-core/mcmicro#165). Group on the channel identity instead.
        .map{ meta -> [[meta.cycle_number, meta.channel_number], meta] }
        .groupTuple()
        .map{ key, rows ->
            def (cycle_number, channel_number) = key
            // Per-sample-varying fields are allowed to differ; ignore them when checking that
            // samples agree on the actual marker definition for this channel.
            def per_sample_fields = ['exposure_time', 'exposure_time_unit']
            def core = rows.collect{ it - it.subMap(per_sample_fields) }.unique()
            if (core.size() != 1) {
                error "Samples disagree on the marker definition for cycle ${cycle_number}, " +
                      "channel ${channel_number}: ${core}."
            }
            if (params.backsub) {
                def exposures = rows.collect{ [it.exposure_time, it.exposure_time_unit] }.unique()
                if (exposures.size() != 1) {
                    error "Samples report inconsistent exposure for cycle ${cycle_number}, " +
                          "channel ${channel_number} (${exposures}); required for background subtraction."
                }
                if (rows.any{ it.exposure_time == null }) {
                    error "Exposure time cannot be null for cycle ${cycle_number}, " +
                          "channel ${channel_number} when doing backsub."
                }
            }
            rows.first()
        }
        .toSortedList { a, b -> a.channel_number <=> b.channel_number }
        .map { rows ->
            // Defense in depth: the shared markers.csv must have exactly one contiguously
            // numbered row per image channel, or MCQUANT will reject it downstream.
            def channels = rows.collect { it.channel_number }
            if (channels != (1..rows.size()).toList()) {
                error "Constructed markers.csv has ${rows.size()} rows numbered ${channels}; " +
                      "expected a contiguous 1..N with one row per image channel."
            }
            rows
        }
        .dump(tag: "ch_markersheet_meta")
        .set { markersheet_meta }

    emit:
    samplesheet = samplesheet_meta
    markersheet = markersheet_meta
}
