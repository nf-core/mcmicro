process ASHLAR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ashlar:1.17.0--pyh5e36f6f_0' :
        'biocontainers/ashlar:1.17.0--pyh5e36f6f_0' }"

    input:
    // Stage everything in its own subdirectory to avoid collisions.
    // dfp and ffp may be empty lists.
    tuple val(meta), path(images, stageAs: 'image*/*'), path(dfp, stageAs: 'dfp*/*'), path(ffp, stageAs: 'ffp*/*')

    output:
    tuple val(meta), path("*.ome.tif"), emit: tif
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def n_images = images instanceof List ? images.size() : 1
    def n_dfp = dfp instanceof List ? dfp.size() : 1
    def n_ffp = ffp instanceof List ? ffp.size() : 1
    if ( dfp && n_dfp != n_images ) { error 'You must provide the same number of dfp paths as image paths' }
    if ( ffp && n_ffp != n_images ) { error 'You must provide the same number of ffp paths as image paths' }

    def args          = task.ext.args           ?: ''
    def prefix        = task.ext.prefix         ?: "${meta.id}"
    def dfp_arg       = dfp                     ? "--dfp $dfp" : ""
    def ffp_arg       = ffp                     ? "--ffp $ffp" : ""
    """

    ashlar \\
        -o ${prefix}.ome.tif \\
        $images \\
        $args \\
        $dfp_arg \\
        $ffp_arg

    sed -i -E 's/UUID="urn:uuid:[[:xdigit:]]{8}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{12}"/                                                    /g' ${prefix}.ome.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ashlar: \$(ashlar --version | sed 's/^.*ashlar //' )
    END_VERSIONS
    """
}
