process METHYLDACKEL_EXTRACT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7' :
        'biocontainers/methyldackel:0.6.1--he4a0461_7' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)

    output:
    tuple val(meta), path("*.bedGraph")             , optional: true, emit: bedgraph
    tuple val(meta), path("*.methylKit")            , optional: true, emit: methylkit
    tuple val(meta), path("*.cytosine_report.txt")  , optional: true, emit: cytosine_report
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def base_args = args.replaceAll('--methylKit', '').replaceAll('--cytosine_report', '').replaceAll('\\s+', ' ').trim()
    def base_args_no_merge = base_args.replaceAll('--mergeContext', '').replaceAll('\\s+', ' ').trim()
    def run_methylkit       = task.ext.run_methylkit       == null ? true : task.ext.run_methylkit
    def run_cytosine_report = task.ext.run_cytosine_report == null ? true : task.ext.run_cytosine_report
    """
    # Run 1: bedGraph (default)
    MethylDackel extract \\
        ${base_args} \\
        $fasta \\
        $bam

    # Run 2: methylKit format (incompatible with --mergeContext)
    if [[ "${run_methylkit}" == "true" ]]; then
        MethylDackel extract \\
            --methylKit \\
            ${base_args_no_merge} \\
            $fasta \\
            $bam
    fi

    # Run 3: cytosine_report (incompatible with --mergeContext)
    if [[ "${run_cytosine_report}" == "true" ]]; then
        MethylDackel extract \\
            --cytosine_report \\
            ${base_args_no_merge} \\
            $fasta \\
            $bam
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methyldackel: \$(MethylDackel --version 2>&1 | cut -f1 -d" ")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch ${bam.baseName}_CpG.bedGraph
    touch ${bam.baseName}_CpG.methylKit
    touch ${bam.baseName}.cytosine_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methyldackel: \$(MethylDackel --version 2>&1 | cut -f1 -d" ")
    END_VERSIONS
    """
}
