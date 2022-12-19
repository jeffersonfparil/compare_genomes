#!/usr/bin/env bash
echo "========================================================="
echo "Running the entire comparative genomics workflow."
echo "========================================================="
echo "SETUP"
time nextflow run modules/setup.nf                               -c config/params.config
echo "DETERMINE ORTHOGROUPS"
time nextflow run modules/orthofinder.nf                         -c config/params.config
echo "BUILD A PHYLOGENETIC TREE USING SINGLE-COPY GENE FAMILIES"
time nextflow run modules/single_gene_orthogroups_tree.nf        -c config/params.config
echo "COMPUTE THE TRANSVERSION RATES AMONNG 4-FOLD DEGENERATE SITES (4DTv)"
time nextflow run modules/single_gene_orthogroups_4DTv.nf        -c config/params.config
echo "ASSESS SIGNIFICANT GENE FAMILY CONTRACTION AND EXPANSION"
time nextflow run modules/gene_family_contraction_expansion.nf   -c config/params.config ### This may take days depending on your machine, number of species, and proteome sizes
echo "GENE ONTOLOGY TERM ENRICHMENT ANALYSIS OF SIGNIFICANTLY EXPANDED GENE FAMILIES"
time nextflow run modules/GO_enrichment.nf                       -c config/params.config
echo "ASSESS WHOLE GENOME DUPLICATION EVENTS"
time nextflow run modules/assess_WGD.nf                          -c config/params.config
echo "PLOT THE PHYLOGENETIC TREE, CONTRACTION/EXPANSION, GENE SETS VENN DIAGRAM AND 4DTv"
time nextflow run modules/plot_tree_conex_venn_4DTv.nf           -c config/params.config
echo "ASSESS CONTRACTION/EXPANSION AND NON-SYNONYMOUS TO SYNONYMOUS NUCLEOTIDE SUBSTITION RATIOS"
time nextflow run modules/assess_specific_genes.nf               -c config/params.config
