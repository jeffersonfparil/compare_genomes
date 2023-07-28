#!/usr/bin/env bash
echo "========================================================="
echo "Running the entire comparative genomics workflow."
echo "========================================================="
nextflow run modules/setup.nf                               -c config/params.config ### SETUP
nextflow run modules/orthofinder.nf                         -c config/params.config ### DETERMINE ORTHOGROUPS
nextflow run modules/single_gene_orthogroups_tree.nf        -c config/params.config ### BUILD A PHYLOGENETIC TREE USING SINGLE-COPY GENE FAMILIES
nextflow run modules/single_gene_orthogroups_4DTv.nf        -c config/params.config ### COMPUTE THE TRANSVERSION RATES AMONNG 4-FOLD DEGENERATE SITES (4DTv)
nextflow run modules/gene_family_contraction_expansion.nf   -c config/params.config ### ASSESS SIGNIFICANT GENE FAMILY CONTRACTION AND EXPANSION (This may take days depending on your machine, number of species, and proteome sizes)
nextflow run modules/GO_enrichment.nf                       -c config/params.config ### GENE ONTOLOGY TERM ENRICHMENT ANALYSIS OF SIGNIFICANTLY EXPANDED GENE FAMILIES
nextflow run modules/assess_WGD.nf                          -c config/params.config ### ASSESS WHOLE GENOME DUPLICATION EVENTS
nextflow run modules/plot_tree_conex_venn_4DTv.nf           -c config/params.config ### PLOT THE PHYLOGENETIC TREE, CONTRACTION/EXPANSION, GENE SETS VENN DIAGRAM AND 4DTv
nextflow run modules/assess_specific_genes.nf               -c config/params.config ### ASSESS CONTRACTION/EXPANSION AND NON-SYNONYMOUS TO SYNONYMOUS NUCLEOTIDE SUBSTITION RATIOS
