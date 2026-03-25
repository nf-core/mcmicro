process PRELUDE_MULTI_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(xmls)

    output:
    path "*_mqc.tsv"                                                            , emit: output_mqc
    tuple val(meta), path("*_errors.tsv")                                                         , emit: output_errors
    tuple val("${task.process}"), val("multi_summary"), eval("python --version"), topic: versions, emit: versions_multisummary

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    multi_summary.py ${xmls}
    """

    stub:
    """
    touch multi_summary_mqc.tsv
    touch multi_summary_errors.tsv
    """
}
