process NEXTCLADE {
    label 'nextclade'

    publishDir "${params.output}/${params.nextclade_dir}", mode: params.publish_dir_mode

    input:
        path input_fasta
        val ref

    output:
        path "variant_table.tsv",   emit: variant_table
        path "sigpep_table.tsv",    emit: sigpep_table

    script:
    """
    echo "Step 1: Annotate mutations in the sequences with Nextclade"

    if [[ (-z "${ref}") && ("${params.subtype}" == "h1n1") ]]
    then
        nextclade dataset get \
        --name 'nextstrain/flu/h1n1pdm/ha/MW626062' \
        --output-dir 'data/h1n1pdm_ha'

        nextclade run \
            --input-dataset data/h1n1pdm_ha \
            --output-all=output/ \
            ${input_fasta}
    elif [[ (! -z "${ref}") && ("${params.subtype}" == "h1n1") ]]
    then
        nextclade dataset get \
        --name 'nextstrain/flu/h1n1pdm/ha/MW626062' \
        --output-dir 'data/h1n1pdm_ha'

        nextclade run \
            --input-dataset data/h1n1pdm_ha \
            --input-ref=${ref} \
            --output-all=output/ \
            ${input_fasta}
    elif [[ (-z "${ref}") && ("${params.subtype}" == "h3n2") ]]
    then
        nextclade dataset get \
        --name 'nextstrain/flu/h3n2/ha/EPI1857216' \
        --output-dir 'data/h3n2_ha'

        nextclade run \
            --input-dataset data/h3n2_ha \
            --output-all=output/ \
            ${input_fasta}
    elif [[ (! -z "${ref}") && ("${params.subtype}" == "h3n2") ]]
    then
        nextclade dataset get \
        --name 'nextstrain/flu/h3n2/ha/EPI1857216' \
        --output-dir 'data/h3n2_ha'

        nextclade run \
            --input-dataset data/h3n2_ha \
            --input-ref=${ref} \
            --output-all=output/ \
            ${input_fasta}
    elif [[ (-z "${ref}") && ("${params.subtype}" == "vic") ]]
    then
        nextclade dataset get \
        --name 'nextstrain/flu/vic/ha/KX058884' \
        --output-dir 'data/vic_ha'

        nextclade run \
            --input-dataset data/vic_ha \
            --output-all=output/ \
            ${input_fasta} 
    elif [[ (! -z "${ref}") && ("${params.subtype}" == "vic") ]]
    then
        nextclade dataset get \
        --name 'nextstrain/flu/vic/ha/KX058884' \
        --output-dir 'data/vic_ha'

        nextclade run \
            --input-dataset data/vic_ha \
            --input-ref=${ref} \
            --output-all=output/ \
            ${input_fasta} 
    fi
    
    nextclade.py \
        -i "output/nextclade.csv" \
        -o "variant_table.tsv" \
        --sigpep "sigpep_table.tsv"
    """

    stub:
    """
    touch variant_table.tsv
    touch sigpep_table.tsv
    """
}