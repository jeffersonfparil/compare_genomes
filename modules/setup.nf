///////////////////////////////////
// SETUP THE WORKING ENVIRONMENT //
///////////////////////////////////
// Assumed extension names:
//  - genome: *.fna
//  - annotation: *.gff
//  - coding DNA: *.cds
//  - proteome: *.faa

process SETUP_DIRECTORIES {
    label "LOW_MEM_LOW_CPU"
    input:
        val dir
        val urls
    output:
        val dir
        val urls
    shell:
    """
    echo $dir
    mkdir ${dir}/GENOMES
    mkdir ${dir}/GFF
    mkdir ${dir}/CDS
    mkdir ${dir}/PROTEOMES
    """
}

// We asume the urls file has no header, where the first column is the name of the output file and the second and last column is the URL of the file
process DOWNLOAD_OMICS_DATA {
    label "LOW_MEM_LOW_CPU"
    input:
        val dir
        val urls
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    for line in $(cat !{urls})
    do
        spe=$(echo ${line} | cut -d, -f1)
        url=$(echo ${line} | cut -d, -f2)
        if [ $(echo ${spe} | rev | cut -d. -f1 | rev) = "fna" ]
        then
            DIR=$(echo !{dir}/GENOMES)
        elif [ $(echo ${spe} | rev | cut -d. -f1 | rev) = "gff" ]
        then
            DIR=$(echo !{dir}/GFF)
        elif [ $(echo ${spe} | rev | cut -d. -f1 | rev) = "cds" ]
        then
            DIR=$(echo !{dir}/CDS)
        else
            DIR=$(echo !{dir}/PROTEOMES)
        fi
        wget ${url} -O - | gunzip -c - > ${DIR}/${spe}
    done
    '''
}

process DOWNLOAD_PANTHER_DATABASE {
    label "LOW_MEM_LOW_CPU"
    input:
        val dir
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    wget 'http://data.pantherdb.org/ftp/panther_library/current_release/PANTHER17.0_hmmscoring.tgz'
    tar -xvzf PANTHER17.0_hmmscoring.tgz
    rm PANTHER17.0_hmmscoring.tgz
    mv target/ PantherHMM_17.0/
    cd PantherHMM_17.0/
    wget 'http://data.pantherdb.org/ftp/hmm_classifications/current_release/PANTHER17.0_HMM_classifications'
    grep -v ':SF' PANTHER17.0_HMM_classifications > Panther17.0_HMM_familyIDs.txt
    '''
}

process INSTALL_JULIA_PACKAGES {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    echo 'using Pkg; Pkg.add(["Plots", "DataFrames", "CSV", "ProgressMeter", "JLD2"])' > install_julia_packages.jl
    julia install_julia_packages.jl
    rm install_julia_packages.jl
    '''
}

process INSTALL_PANTHER_API_FORGO {
    label "LOW_MEM_LOW_CPU"
    input:
        val dir
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    git clone https://github.com/pantherdb/pantherapi-pyclient.git
    cd pantherapi-pyclient
    python3 -m venv env
    # . env/bin/activate (bash) or source env/bin/activate.csh (C-shell or tcsh)
    pip3 install -r requirements.txt
    # test: python3 pthr_go_annots.py --service enrich --params_file params/enrich.json --seq_id_file resources/test_ids.txt
    '''
}

workflow {
    SETUP_DIRECTORIES(params.dir, params.urls) | \
        DOWNLOAD_OMICS_DATA
    DOWNLOAD_PANTHER_DATABASE(params.dir)
    INSTALL_JULIA_PACKAGES(params.dir)
    INSTALL_PANTHER_API_FORGO(params.dir)
}
