process ALIGNMULTIPLECYCLES {
    tag "$meta.id"
    label 'process_single'
    container "docker.io/voigtgesa/palom:v2024.4.1_cli2"

    input:
    tuple val(meta), path(images, stageAs: 'image*/*')

    output:
    tuple val(meta), path("*.ome.tif"), emit: tif
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args           ?: ''
    def prefix        = task.ext.prefix         ?: "${meta.id}"

    """
    python /palom/align_multiple_cycles.py \\
        --img_list $images \\
        --out_name ${prefix}.ome.tif \\
        --out_dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        palom: \$(palom-svs --version | sed 's/palom v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ome.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        palom: \$(palom-svs --version | sed 's/palom v//')
    END_VERSIONS
    """
}
