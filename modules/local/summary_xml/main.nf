process PRELUDE_SUMMARY_XML {
    tag "${meta.id}_${meta.cycle_number}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(xml)

    output:
    tuple val(meta), path("*_xml_mqc.tsv")                                    , emit: output
    tuple val(meta), path("*_variables.tsv")                                  , emit: variables
    tuple val("${task.process}"), val("summary_xml"), eval("python --version"), emit: versions_summaryxml, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.cycle_number}"
    """
    summary_xml.py ${xml} ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.cycle_number}"
    """
    touch ${prefix}_xml_mqc.tsv
    touch ${prefix}_variables.tsv
    """
}
