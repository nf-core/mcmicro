process PRELUDE_SUMMARY_SAMPLESHEET {
    tag "${meta.id}_${meta.cycle_number}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(samplesheet)

    output:
    tuple val(meta), path("*.tsv"), emit: output
    tuple val("${task.process}"), val("summary_samplesheet"), eval("python --version"), topic: versions, emit: versions_summarysample

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.cycle_number}"
    """
    summary_samplesheet.py ${samplesheet} ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.cycle_number}"
    """
    touch ${prefix}_samplesheet_mqc.tsv
    """
}
