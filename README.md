# compare_genomes: a comparation genomics workflow


- Comparative genomics paper seldomly fully disclose their comparative genomics workflow
- This hinders replicability and transferability of novel methods
- Here we use conda and nextflow to improvie transferability of our specific comparative genomics workflow

## Configuration
1. Install conda
```shell
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh ./Miniconda3-latest-Linux-x86_64.sh
### log-out then log-in
```
2. Configure conda
```shell
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

3. Create a new conda environment
```shell
conda create --name "adder"
conda env list
conda activate adder
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
```

5. Export environment and create a new environment based on the exported settings
```shell
# Export
conda env export -n adder > adder.yml
```

6. Import conda environment
```shell
conda env create -n adder_import --file adder.yml
```

7. Install nextflow but first install java
```shell
sudo apt install -y default-jre
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/bin
```

8. Dowload repo and run
```shell
git clone https://github.com/jeffersonfparil/compare_genomes.git
cd compare_genomes
nextflow run modules/setup.nf                             -c config/params.config
nextflow run modules/orthofinder.nf                       -c config/params.config
nextflow run modules/gene_family_contraction_expansion.nf -c config/params.config ### Thismay take days depending on your machine, number of species, and proteome sizes
nextflow run modules/GO_enrichment.nf                     -c config/params.config
nextflow run modules/single_gene_orthogroups_tree.nf      -c config/params.config
```

