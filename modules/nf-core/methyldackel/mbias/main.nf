process METHYLDACKEL_MBIAS {
    tag "$meta.id"
    label 'process_low'

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
    tuple val(meta), path("*.mbias.txt"), emit: txt
    tuple val(meta), path("*.svg")      , optional: true, emit: svg
    path  "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def run_all_contexts = task.ext.all_contexts ? 'true' : 'false'
    """
    MethylDackel mbias \\
        $args \\
        $fasta \\
        $bam \\
        ${prefix}_CpG \\
        --txt \\
        > ${prefix}_CpG.mbias.txt

    if [[ "${run_all_contexts}" == "true" ]]; then
        MethylDackel mbias \\
            $args \\
            --CHG --noCpG \\
            $fasta \\
            $bam \\
            ${prefix}_CHG \\
            --txt \\
            > ${prefix}_CHG.mbias.txt

        MethylDackel mbias \\
            $args \\
            --CHH --noCpG \\
            $fasta \\
            $bam \\
            ${prefix}_CHH \\
            --txt \\
            > ${prefix}_CHH.mbias.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methyldackel: \$(MethylDackel --version 2>&1 | cut -f1 -d" ")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def run_all_contexts = task.ext.all_contexts ? 'true' : 'false'
    """
    touch ${prefix}_CpG.mbias.txt
    touch ${prefix}_CpG_OT.svg
    touch ${prefix}_CpG_OB.svg

    if [[ "${run_all_contexts}" == "true" ]]; then
        touch ${prefix}_CHG.mbias.txt
        touch ${prefix}_CHG_OT.svg
        touch ${prefix}_CHG_OB.svg
        touch ${prefix}_CHH.mbias.txt
        touch ${prefix}_CHH_OT.svg
        touch ${prefix}_CHH_OB.svg
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methyldackel: \$(MethylDackel --version 2>&1 | cut -f1 -d" ")
    END_VERSIONS
    """
}
