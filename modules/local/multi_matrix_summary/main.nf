process PRELUDE_MULTI_MATRIX_SUMMARY {
    tag "matrix_summary"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    path(sample)
    path(xml)
    path(multi)

    output:
    path("*_mqc.tsv")                                                         , emit: output
    tuple val("${task.process}"), val("summary_xml"), eval("python --version"), topic: versions, emit: versions_summarymatrix

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "matrix_summary"
    """
    multi_matrix_summary.py --xml ${xml} --samplesheet ${sample} --merged ${multi} --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "matrix_summary"
    """
    touch ${prefix}_mqc.tsv
    """
}
