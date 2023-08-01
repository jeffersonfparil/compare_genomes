/////////////////////////////////////////////////
// INFER GENE FAMILY CONTRACTION AND EXPANSION //
/////////////////////////////////////////////////
// Assumed extension names:
//  - genome: *.fna
//  - annotation: *.gff
//  - coding DNA: *.cds
//  - proteome: *.faa

process CAFE5_GENE_FAMILY_CONTRACTION_EXPANSION {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
        val cafe5_n_gamma_cats
        val cafe5_pvalue
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    
    echo "Define the location of the results of OrthoFinder run, i.e. the most recent output folder."
    DIR_ORTHOFINDER_OUT=$(ls -tr ORTHOGROUPS/OrthoFinder/ | tail -n1)
    DIR_ORTHOGROUPS=$(pwd)/ORTHOGROUPS/OrthoFinder/${DIR_ORTHOFINDER_OUT}

    ORTHOUT=$(pwd)/ORTHOGROUPS/orthogroups_gene_counts_families_go.out
    rev ${ORTHOUT} | cut -f5- | rev > col2_to_coln.tmp
    awk -F'\t' '{print $(NF-1)}' ${ORTHOUT} > col1.tmp
    paste -d'\t' col1.tmp col2_to_coln.tmp > counts.tmp
    # TREE=${DIR_ORTHOGROUPS}/Species_Tree/SpeciesTree_rooted.txt ### May not be generated successfully due to branch conflicts/inconsistencies
    TREE=ORTHOGROUPS_SINGLE_GENE.NT.treefile ### NOTE: Make sure there are no redundant species or else the following analysis will not be successful.

    echo "Remove the top 1% of orthogroups with the highest between species differences in an attempt to reduce too much bias."
    echo 'using Statistics
    file = open("counts.tmp", "r")
    vec_name = []
    vec_id = []
    vec_count = []
    vec_std = []
    vec_with_zero = []
    header = readline(file)
    while !eof(file)
        line = split(readline(file), "\t")
        push!(vec_name, line[1])
        push!(vec_id, line[2])
        push!(vec_count, parse.(Int64, line[3:end]))
        push!(vec_std, std(vec_count[end]))
    end
    close(file)
    n = length(vec_std)
    temp_idx = sortperm(vec_std, rev=true)
    temp_std = vec_std[temp_idx]
    maximum_threshold = temp_std[Int64(ceil(n*0.01))]
    vec_idx = vec_std .< maximum_threshold
    # sum(vec_idx)
    # using UnicodePlots
    # UnicodePlots.histogram(Float64.(vec_std[vec_idx]))
    # Write-out the new counts file
    file = open("counts_filtered.tmp", "a")
    write(file, string(header, "\n"))
    for i in 1:length(vec_std)
        # i = 50
        if vec_idx[i]
            line = join([join([vec_name[i], vec_id[i]], "\t"), join(vec_count[i], "\t")], "\t")
            write(file, string(line, "\n"))
        end
    end
    close(file)
    ' > remove_max_var_orthogroup.jl
    julia remove_max_var_orthogroup.jl

    echo "Gene family contraction/expansion analysis."
    if [ !{cafe5_n_gamma_cats} -gt 1 ]
    then
        cafe5 \
            --infile counts_filtered.tmp \
            --tree ${TREE} \
            --n_gamma_cats !{cafe5_n_gamma_cats} \
            --cores !{task.cpus} \
            --pvalue !{cafe5_pvalue} \
            --output_prefix CAFE_results
        echo -e "Species\tExpansion\tContraction" > CONTRACTION_EXPANSION.txt
        grep -v "^#" CAFE_results/Gamma_clade_results.txt | \
            grep -v "^<" | \
            sed 's/<..>//g' | \
            sed 's/<.>//g' >> CONTRACTION_EXPANSION.txt
        ### If something goes wrong, then revert to the base model.      
        if [ $(cat CONTRACTION_EXPANSION.txt | wc -l) -eq 1 ]
        then
            cafe5 \
                --infile counts_filtered.tmp \
                --tree ${TREE} \
                --cores !{task.cpus} \
                --pvalue !{cafe5_pvalue} \
                --output_prefix CAFE_results
            echo -e "Species\tExpansion\tContraction" > CONTRACTION_EXPANSION.txt
            grep -v "^#" CAFE_results/Base_clade_results.txt | \
                grep -v "^<" | \
                sed 's/<..>//g' | \
                sed 's/<.>//g' >> CONTRACTION_EXPANSION.txt
        fi
    else
        cafe5 \
            --infile counts_filtered.tmp \
            --tree ${TREE} \
            --cores !{task.cpus} \
            --pvalue !{cafe5_pvalue} \
            --output_prefix CAFE_results
        echo -e "Species\tExpansion\tContraction" > CONTRACTION_EXPANSION.txt
        grep -v "^#" CAFE_results/Base_clade_results.txt | \
            grep -v "^<" | \
            sed 's/<..>//g' | \
            sed 's/<.>//g' >> CONTRACTION_EXPANSION.txt
    fi

    echo "Cleanup"
    rm col1.tmp counts.tmp counts_filtered.tmp

    echo "Output:"
    echo "  (1/1) ORTHOGROUPS_SINGLE_GENE.NT.4DTv"
    '''
}

workflow {
    CAFE5_GENE_FAMILY_CONTRACTION_EXPANSION(params.dir,
                                            params.cafe5_n_gamma_cats,
                                            params.cafe5_pvalue)
}
