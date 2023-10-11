process VOCAL {
    label 'vocal'

    publishDir "${params.output}/${params.vocal_dir}", mode: params.publish_dir_mode

    input:
        val mode
        path input_fasta
        path output_psl

    output:
        path "variant_table.tsv",   emit: variant_table

    script:

    def filter = opt.name != 'NO_FILE' ? "--filter $output_psl" : ''

    """
    echo "Step 1: Annotate mutations in the proteins"

    if [ ${mode} = 'PSL' ]
    then
        python vocal/vocal.py \
            -i ${input_fasta} \
            --PSL ${filter} \
            -o "variant_table.tsv"
    else
        python vocal/vocal.py \
            -i ${input_fasta} \
            -o "variant_table.tsv"
    fi
    """

    stub:
    """
    touch variant_table.tsv
    """
}