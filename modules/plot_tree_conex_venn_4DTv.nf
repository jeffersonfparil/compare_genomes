////////////////////////////
// ASSESS SPECIFIC GENES  //
////////////////////////////

process PLOT {
    label "LOW_MEM_LOW_CPU"
    input:
        val dir
        val comparisons_4DTv
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    TREE=ORTHOGROUPS_SINGLE_GENE.NT.treefile
    Rscript !{projectDir}/../scripts/plot_tree_singleGeneConex_venn_4DTv.R \
        ORTHOGROUPS_SINGLE_GENE.NT.timetree.nex \
        CONTRACTION_EXPANSION.txt \
        PROTEOMES/orthogroups_summarised_gene_counts.csv \
        PROTEOMES/orthogroups_gene_counts_families_go.out \
        !{dir} \
        .4DTv \
        ORTHOGROUPS_SINGLE_GENE.NT.4DTv \
        !{comparisons_4DTv} \
        !{params.species_of_interest}_comparative_genomics.svg
    '''
}

workflow {
    PLOT(params.dir, params.comparisons_4DTv)
}
