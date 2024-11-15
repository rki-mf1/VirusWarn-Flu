process NEXTCLADE {
    label 'nextclade'

    publishDir "${params.output}/${params.nextclade_dir}", mode: params.publish_dir_mode

    input:
        path input_fasta
        val subtype

    output:
        path "variant_table.tsv",       emit: variant_table
        path "output/nextclade.csv",    emit: nextclade
        path "sigpep_table.tsv",        emit: sigpep_table

    script:
    """
    echo "Step 1: Annotate mutations in the sequences with Nextclade"

    if [[ (-z "${params.ref}") && ("${subtype}" == "h1n1") ]]
    then
        echo "Reference used for Influenza A H1N1: A/Wisconsin/588/2019 (MW626065)"

        nextclade dataset get \
            --name 'nextstrain/flu/h1n1pdm/ha/MW626062' \
            --output-dir 'data/h1n1pdm_ha'

        nextclade run \
            --input-dataset data/h1n1pdm_ha \
            --output-all=output/ \
            ${input_fasta}
    elif [[ ("${params.ref}" == "old") && ("${subtype}" == "h1n1") ]]
    then
        echo "Reference used for Influenza A H1N1: A/California/7/2009 (CY121680)"

        nextclade dataset get \
            --name 'nextstrain/flu/h1n1pdm/ha/CY121680' \
            --output-dir 'data/h1n1pdm_ha'

        nextclade run \
            --input-dataset data/h1n1pdm_ha \
            --output-all=output/ \
            ${input_fasta}
    elif [[ (-z "${params.ref}") && ("${subtype}" == "h3n2") ]]
    then
        echo "Reference used for Influenza A H3N2: A/Darwin/6/2021 (EPI1857216)"

        nextclade dataset get \
            --name 'nextstrain/flu/h3n2/ha/EPI1857216' \
            --output-dir 'data/h3n2_ha'

        nextclade run \
            --input-dataset data/h3n2_ha \
            --output-all=output/ \
            ${input_fasta}
    elif [[ ("${params.ref}" == "old") && ("${subtype}" == "h3n2") ]]
    then
        echo "Reference used for Influenza A H3N2: A/Wisconsin/67/2005 (CY163680)"

        nextclade dataset get \
            --name 'nextstrain/flu/h3n2/ha/CY163680' \
            --output-dir 'data/h3n2_ha'

        nextclade run \
            --input-dataset data/h3n2_ha \
            --output-all=output/ \
            ${input_fasta}
    elif [[ (-z "${params.ref}") && ("${subtype}" == "vic") ]]
    then
        echo "Reference used for Influenza B Victoria: B/Brisbane/60/2008 (KX058884)"

        nextclade dataset get \
            --name 'nextstrain/flu/vic/ha/KX058884' \
            --output-dir 'data/vic_ha'

        nextclade run \
            --input-dataset data/vic_ha \
            --output-all=output/ \
            ${input_fasta} 
    fi
    
    nextclade.py \
        -i "output/nextclade.csv" \
        --subtype ${subtype} \
        -o "variant_table.tsv" \
        --sigpep "sigpep_table.tsv"
    """

    stub:
    """
    touch variant_table.tsv
    touch sigpep_table.tsv
    """
}