process VOCAL {
    label 'vocal'

    publishDir "${params.output}/${params.vocal_dir}", mode: params.publish_dir_mode

    input:
        val mode
        path input_fasta
        path ref_nt
        path ref_aa
        path output_psl

    output:
        path "variant_table.tsv",   emit: variant_table

    script:
    
    // to make output_psl optional
    def filter = opt.name != 'NO_FILE' ? "--filter $output_psl" : ''

    """
    echo "Step 1: Annotate mutations in the proteins"

    if [ ${mode} = 'PSL' ]
    then
        python vocal/vocal.py \
            -i ${input_fasta} \
            --ref_nt ${ref_nt} \
            --ref_aa ${ref_aa} \
            --PSL ${filter} \
            -o "variant_table.tsv"
    else
        python vocal/vocal.py \
            -i ${input_fasta} \
            --ref_nt ${ref_nt} \
            --ref_aa ${ref_aa} \
            -o "variant_table.tsv"
    fi
    """

    stub:
    """
    touch variant_table.tsv
    """
}