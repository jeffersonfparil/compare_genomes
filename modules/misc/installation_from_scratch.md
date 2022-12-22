# Installation from scratch

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
conda install -y nextflow
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

7. Dowload this comparative_genomics repository
```shell
git clone https://github.com/jeffersonfparil/compare_genomes.git
```