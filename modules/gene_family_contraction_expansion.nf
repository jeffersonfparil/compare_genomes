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
    DIR_ORTHOFINDER_OUT=$(ls -tr PROTEOMES/OrthoFinder/ | tail -n1)
    DIR_ORTHOGROUPS=$(pwd)/PROTEOMES/OrthoFinder/${DIR_ORTHOFINDER_OUT}

    ORTHOUT=$(pwd)/PROTEOMES/orthogroups_gene_counts_families_go.out
    rev ${ORTHOUT} | cut -f5- | rev > col2_to_coln.tmp
    awk -F'\t' '{print $(NF-1)}' ${ORTHOUT} > col1.tmp
    paste -d'\t' col1.tmp col2_to_coln.tmp > counts.tmp
    # TREE=${DIR_ORTHOGROUPS}/Species_Tree/SpeciesTree_rooted.txt
    TREE=ORTHOGROUPS_SINGLE_GENE.NT.timetree.nwk
    cafe5 \
        --infile counts.tmp \
        --tree ${TREE} \
        --n_gamma_cats !{cafe5_n_gamma_cats} \
        --cores tasks.cpus \
        --pvalue !{cafe5_pvalue} \
        --output_prefix CAFE_results
    echo -e "Species\tExpansion\tContraction" > CONTRACTION_EXPANSION.txt
    grep -v "^#" CAFE_results/Gamma_clade_results.txt | \
        grep -v "^<" | \
        sed 's/<..>//g' | \
        sed 's/<.>//g' >> CONTRACTION_EXPANSION.txt

    echo "Cleanup"
    rm col1.tmp counts.tmp

    echo "Output:"
    echo "  (1/1) ORTHOGROUPS_SINGLE_GENE.NT.4DTv"
    '''
}

workflow {
    CAFE5_GENE_FAMILY_CONTRACTION_EXPANSION(params.dir,
                                            params.cafe5_n_gamma_cats,
                                            params.cafe5_pvalue)
}
