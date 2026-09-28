process PRELUDE_SUMMARY_MARKERSHEET {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(markersheet)

    output:
    tuple val(meta), path("*.tsv")                                                    , emit: output
    tuple val("${task.process}"), val("summary_markersheet"), eval("python --version"), emit: versions_summarymarker, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"

    """
    summary_markersheet.py ${markersheet} ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_markersheet_mqc.tsv
    """
}
