process MCCELLPOSE {
    tag "$meta.id"
    label 'process_low'
    label 'process_gpu'

    container "docker.io/labsyspharm/mccellpose:1.0.3"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*_mask.ome.tif"), emit: mask
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--channel 1'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gpu_args = task.ext.use_gpu ? "--use-gpu --jobs ${task.cpus}" : ''
    """
    export HOME=\$PWD
    export NUMBA_CACHE_DIR=\$PWD

    mccellpose \
        --input $image \
        --output-cell ${prefix}_mask.ome.tif \
        --channel 1 \
        --expand-size 2 \
        $gpu_args \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mccellpose: \$(mccellpose --version | awk '{print \$2}')
        cellpose: \$(cellpose --version | awk 'NR==2 {print \$3}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mask.ome.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mccellpose: \$(mccellpose --version | awk '{print \$2}')
        cellpose: \$(cellpose --version | awk 'NR==2 {print \$3}')
    END_VERSIONS
    """
}
