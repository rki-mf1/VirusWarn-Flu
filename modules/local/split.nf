process SPLIT {
    label 'split'

    publishDir "${params.output}/${params.split_dir}", mode: params.publish_dir_mode

    input:
        path input_fasta

    output:
        path "h1n1_ha.fasta",   emit: h1n1_ha,      optional: true
        path "h3n2_ha.fasta",   emit: h3n2_ha,      optional: true
        path "vic_ha.fasta",   emit: vic_ha,        optional: true

    script:
    """
    echo "Preprocessing Step: Splitting fasta file into seperate files for subtypes / segments"

    splitting.py \
        -i ${input_fasta} \
        -m ${params.split} \
        --h1n1 "h1n1_ha.fasta" \
        --h3n2 "h3n2_ha.fasta" \
        --vic "vic_ha.fasta"
    """

    stub:
    """
    touch h1n1_ha.fasta
    touch h3n2_ha.fasta
    touch vic_ha.fasta
    """
}