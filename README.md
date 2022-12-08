# compare_genomes: a comparative genomics workflow

- Comparative genomics papers seldomly fully disclose their comparative genomics workflow
- This hinders replicability and transferability of novel methods
- Here we use conda and nextflow to improve the transferability of our specific comparative genomics workflow which have been used in the comparative genomics of new [*Lolium rigidum*](https://doi.org/10.3389/fgene.2022.1012694) and *Salvia hispanica* genome assemblies.

## Installation
1. Install conda
```shell
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh ./Miniconda3-latest-Linux-x86_64.sh
echo "To finish the conda setup - log-out then log back in."
```
2. Configure conda
```shell
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

3. Create a new conda environment
```shell
conda create --name "compare_genomes"
conda env list
conda activate compare_genomes
```

4. Install software
```shell
conda install -y parallel
conda install -y wget r-base julia
conda install -y -c bioconda orthofinder
conda install -y -c bioconda hmmer
conda install -y -c bioconda cafe
conda install -y -c bioconda macse
conda install -y -c bioconda iqtree
conda install -y -c bioconda kakscalculator2
conda install -y -c conda-forge r-ape
conda install -y -c conda-forge r-VennDiagram
conda install -y -c conda-forge r-png
```

5. Export environment and create a new environment based on the exported settings
```shell
conda env export -n compare_genomes > compare_genomes.yml
```

6. Import conda environment
```shell
conda env create -n compare_genomes_import --file compare_genomes.yml
```

7. Install nextflow but first install java
```shell
sudo apt install -y default-jre
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/bin
```

8. Dowload this comparative_genomics repository
```shell
git clone https://github.com/jeffersonfparil/compare_genomes.git
```

## Setup parameters

Setup the parameters for the comparative genomics analysis (Note that for gene family contraction and expansion analysis to work, you need to include at least 3 species).

1. [`urls.txt`](config/urls.txt): links to the genome sequence, genome annotation, coding DNA sequence, and amino acid sequences for at least 3 species you wish to include in the analyses.
    - Formatted as headerless, two-columned, comma-separated file
    - Column 1: filename of the genome sequence, genome annotation, coding DNA sequence, and amino acid sequences (**Note**: the species names and extension names should be the consistent across these files)
    - Column 2: URL (uniform resource locator) of the gunzipped files (*.gz) for download

2. [`dates.txt`](config/dates.txt): pairwise divergence times between the species you wish to include in the analyses (e.g. look up dvergence times, e.g. from [http://timetree.org/](http://timetree.org/))
    - Formatted as headerless, two-columned, tab-delimited file
    - Column 1: two species separated by a comma with the same names used in the [`urls.txt`](config/urls.txt)
    - Column 2: time in million years, e.g. -160 for 160 million years ago

3. [`genes.txt`](config/genes.txt): links to the gene sequences you wish to look at individually across the species included
    - Formatted as headerless, three-columned, comma-separated file
    - Column 1: phenotype name or some identifying name
    - Column 2: species name which should be the same as in [`urls.txt`](config/urls.txt) and [`dates.txt`](config/dates.txt)
    - Column 3: URL of the genes for download

4. [`params.config`](config/params.config): configuration file listing the variables specific to the analyses you wish to perform
    - **dir**: output directory.
    - **species_of_interest**: the focal species of interest which should the same as in [`urls.txt`](config/urls.txt), [`dates.txt`](config/dates.txt), and [`genes.txt`](config/genes.txt).
    - **species_of_interest_panther_HMM_for_gene_names_url**: URL to the specific Panther HMM database to extract gene names from, preferrably from the same species which will be used for gene ontology (GO) term enrichment analysis. See the current release list [here](http://data.pantherdb.org/ftp/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/).
    - **urls**: location of [`urls.txt`](config/urls.txt).
    - **dates** = location of [`dates.txt`](config/dates.txt), e.g. look up dvergence times from [http://timetree.org/](http://timetree.org/).
    - **genes**: location of [`genes.txt`](config/genes.txt).
    - **genomes**: extension name of the genome sequences (e.g. consistently '*.fna' for all species).
    - **gff**: extension name of the genome annotations (e.g. consistently '*.gff' for all species).
    - **cds**: extension name of the coding DNA sequences (e.g. consistently '*.cds' for all species).
    - **faa**: extension name of the protein sequences (e.g. consistently '*.faa' for all species).
    - **cafe5_n_gamma_cats**: number of the gamma values (parameter of the substittion model) to use for the assessment of significant gene family expansion and contraction using CAFE5.
    - **cafe5_pvalue**: signifcance threshold of gene family expansion and contraction.
    - **go_term_enrich_genome_id**: genome ID for the species specified in **species_of_interest_panther_HMM_for_gene_names_url**. Find the appropriate taxon ID from [here](http://pantherdb.org/services/oai/pantherdb/supportedgenomes).
    - **go_term_enrich_annotation_id**: code for the gene ontology level you with to use, e.g. `GO:0008150` for "Biological Process". See the list of GO codes [here](http://pantherdb.org/services/oai/pantherdb/supportedannotdatasets).
    - **go_term_enrich_test**: GO term enrichment test to perform which can be either `FISHER` (Fisher's Exact Test) or `BINOMIAL`" (binomial distribution test).
    - **go_term_enrich_correction**: multiple testing correction which can be `NONE`, `FDR` (False discovery rate), or `BONFERRONI` (Bonferroni correction).
    - **go_term_enrich_ngenes_per_test**: number of genes to include in each GO term enrichment analysis.
    - **go_term_enrich_ntests**: number GO term enrichment analyses to perform.

5. [`process.config`](config/process.config): the second and last cofiguration file listing the computing resource allocation availble to you. Assign the number of **cpus** and **memory** capacity to use for low and high resources intensive tasks:
    - **LOW_MEM_LOW_CPU**
    - **HIGH_MEM_HIGH_CPU**

## Run

1. Run the entire workflow
```shell
cd compare_genomes
chmod +x run.sh
./run.sh
```

2. Run each module individually but this needs to be in order as each module depends on the output of the previous one. Run these if you wish to troubleshoot and/or modify the workflow for analyses not implemented here, e.g. add a GO term enrichment analysis for significantly contracted gene families.
```shell
cd compare_genomes
nextflow run modules/setup.nf                               -c config/params.config ### SETUP
nextflow run modules/orthofinder.nf                         -c config/params.config ### DETERMINE ORTHOGROUPS
nextflow run modules/gene_family_contraction_expansion.nf   -c config/params.config ### ASSESS SIGNIFICANT GENE FAMILY CONTRACTION AND EXPANSION (This may take days depending on your machine, number of species, and proteome sizes)
nextflow run modules/GO_enrichment.nf                       -c config/params.config ### GENE ONTOLOGY TERM ENRICHMENT ANALYSIS OF SIGNIFICANTLY EXPANDED GENE FAMILIES
nextflow run modules/single_gene_orthogroups_tree.nf        -c config/params.config ### BUILD A PHYLOGENETIC TREE USING SINGLE-COPY GENE FAMILIES
nextflow run modules/single_gene_orthogroups_4DTv.nf        -c config/params.config ### COMPUTE THE TRANSVERSION RATES AMONNG 4-FOLD DEGENERATE SITES (4DTv)
nextflow run modules/assess_WGD.nf                          -c config/params.config ### ASSESS WHOLE GENOME DUPLICATION EVENTS
nextflow run modules/plot_tree_conex_venn_4DTv.nf           -c config/params.config ### PLOT THE PHYLOGENETIC TREE, CONTRACTION/EXPANSION, GENE SETS VENN DIAGRAM AND 4DTv
nextflow run modules/assess_specific_genes.nf               -c config/params.config ### ASSESS CONTRACTION/EXPANSION AND NON-SYNONYMOUS TO SYNONYMOUS NUCLEOTIDE SUBSTITION RATIOS
```

## References
- OrthoFinder: Emms, David M., and Steven Kelly. “OrthoFinder: Phylogenetic Orthology Inference for Comparative Genomics.” Genome Biology 20, no. 1 (November 14, 2019): 238. https://doi.org/10.1186/s13059-019-1832-y.was 
- HMMER: Mistry, Jaina, Robert D. Finn, Sean R. Eddy, Alex Bateman, and Marco Punta. “Challenges in Homology Search: HMMER3 and Convergent Evolution of Coiled-Coil Regions.” Nucleic Acids Research 41, no. 12 (July 1, 2013): e121. https://doi.org/10.1093/nar/gkt263.
- PantherHMM gene family models: Mi, Huaiyu, Anushya Muruganujan, Dustin Ebert, Xiaosong Huang, and Paul D Thomas. “PANTHER Version 14: More Genomes, a New PANTHER GO-Slim and Improvements in Enrichment Analysis Tools.” Nucleic Acids Research 47, no. D1 (January 8, 2019): D419–26. https://doi.org/10.1093/nar/gky1038.
- CAFE v5: De Bie, Tijl, Nello Cristianini, Jeffery P. Demuth, and Matthew W. Hahn. “CAFE: A Computational Tool for the Study of Gene Family Evolution.” Bioinformatics 22, no. 10 (May 15, 2006): 1269–71. https://doi.org/10.1093/bioinformatics/btl097.
- Gene ontology (GO) enrichment analysis: The UniProt Consortium. “UniProt: A Worldwide Hub of Protein Knowledge.” Nucleic Acids Research 47, no. D1 (January 8, 2019): D506–15. https://doi.org/10.1093/nar/gky1049.
- MACSE: Ranwez, Vincent, Sébastien Harispe, Frédéric Delsuc, and Emmanuel J. P. Douzery. “MACSE: Multiple Alignment of Coding SEquences Accounting for Frameshifts and Stop Codons.” PLOS ONE 6, no. 9 (September 16, 2011): e22594.
- IQ-TREE: Minh, Bui Quang, Heiko A Schmidt, Olga Chernomor, Dominik Schrempf, Michael D Woodhams, Arndt von Haeseler, and Robert Lanfear. “IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era.” Molecular Biology and Evolution 37, no. 5 (May 1, 2020): 1530–34. https://doi.org/10.1093/molbev/msaa015.
- ModelFinder (module of IQ-TREE): Kalyaanamoorthy, Subha, Minh, Bui Quang, Wong, Thomas KF, von Haeseler, Arndt, and Jermiin, Lars S. ”ModelFinder: Fast model selection for accurate phylogenetic estimates.” Nature Methods, 14 (2017):587–589. https://doi.org/10.1038/nmeth.4285
- KaKs_calculator v2: Wang, Da-Peng, Hao-Lei Wan, Song Zhang, and Jun Yu. “γ-MYN: A New Algorithm for Estimating Ka and Ks with Consideration of Variable Substitution Rates.” Biology Direct 4, no. 1 (June 16, 2009): 20. https://doi.org/10.1186/1745-6150-4-20.

