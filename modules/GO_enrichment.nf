/////////////////////////////////////////////////////////////////////
// GENE ONTOLOGY TERM ENRICHMENT ANALSIS OF EXPANDED GENE FAMILIES //
/////////////////////////////////////////////////////////////////////
// Assumption:
//  - !{species_of_interest} matches of the the species in "config/urls.txt"

process GO_TERM_ENRICHMENT {
    label "LOW_MEM_LOW_CPU"
    input:
        val dir
        val species_of_interest
        val species_of_interest_panther_HMM_for_gene_names_url
        val go_term_enrich_genome_id
        val go_term_enrich_annotation_id
        val go_term_enrich_test
        val go_term_enrich_correction
        val go_term_enrich_ngenes_per_test
        val go_term_enrich_ntests
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    ORTHOUT=PROTEOMES/orthogroups_gene_counts_families_go.out
    n=$(head -n1 CAFE_results/Gamma_change.tab | sed -z "s/\\t/\\n/g" | grep -n "!{species_of_interest}" | cut -d":" -f1)
    cut -f1,${n} CAFE_results/Gamma_change.tab | grep -v "+0" | grep "+" | cut -f1 > expanded_orthogroups_for_grep.tmp
    cut -f1,${n} CAFE_results/Gamma_change.tab | grep -v "+" | cut -f1 > contracted_orthogroups_for_grep.tmp

    grep -wf expanded_orthogroups_for_grep.tmp $ORTHOUT | cut -f$(head -n1 $ORTHOUT | awk '{printf NF-2}') | grep -v "^$" > expanded_orthogroups.pthr.tmp
    grep -wf contracted_orthogroups_for_grep.tmp $ORTHOUT | cut -f$(head -n1 $ORTHOUT | awk '{printf NF-2}') | grep -v "^$" > contracted_orthogroups.pthr.tmp
    wget !{species_of_interest_panther_HMM_for_gene_names_url} -O PTHR_unzipped
    grep -wf expanded_orthogroups.pthr.tmp PTHR_unzipped | cut -f3 | sed "s/^LOC_//g"> expanded_orthogroups.forgo
    grep -wf contracted_orthogroups.pthr.tmp PTHR_unzipped | cut -f3 > contracted_orthogroups.forgo

    echo "{
        \\"organism\\": \\"!{go_term_enrich_genome_id}\\",
        \\"refOrganism\\": \\"!{go_term_enrich_genome_id}\\",
        \\"annotDataSet\\": \\"!{go_term_enrich_annotation_id}\\",
        \\"enrichmentTestType\\": \\"!{go_term_enrich_test}\\",
        \\"correction\\": \\"!{go_term_enrich_correction}\\"
    }" > go_term_enrich.json

    ### Shuffle lines
    for i in $(seq 1 !{go_term_enrich_ntests})
    do
        shuf \
            -n !{go_term_enrich_ngenes_per_test} \
            expanded_orthogroups.forgo > \
            expanded_orthogroups.forgo.${i}.shuffled
        python3 pantherapi-pyclient/pthr_go_annots.py \
            --service enrich \
            --params_file go_term_enrich.json \
            --seq_id_file expanded_orthogroups.forgo.${i}.shuffled > \
            expanded_orthogroups.${i}.goout
    done
    rm PTHR_unzipped *.tmp
    '''
}

workflow {
    GO_TERM_ENRICHMENT(params.dir,
                       params.species_of_interest,
                       params.species_of_interest_panther_HMM_for_gene_names_url,
                       params.go_term_enrich_genome_id,
                       params.go_term_enrich_annotation_id,
                       params.go_term_enrich_test,
                       params.go_term_enrich_correction,
                       params.go_term_enrich_ngenes_per_test,
                       params.go_term_enrich_ntests)
}

