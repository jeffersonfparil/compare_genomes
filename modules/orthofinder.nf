/////////////////////////////////
// FIND ORTHOLOGS AND PARALOGS //
/////////////////////////////////
// Assumed extension names:
//  - genome: *.fna
//  - annotation: *.gff
//  - coding DNA: *.cds
//  - proteome: *.faa

process FIND_ORTHOGROUPS {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    echo "Append species names into protein sequence names."
    for f in PROTEOMES/*.faa
    do
        # f=$(ls ORTHOGROUPS/*.faa | head -n1)
        fname=$(basename ${f})
        species=${fname%.faa*}
        sed -i "s/^>/>$species|/g" $f
    done
    echo "Run OrthoFinder."
    orthofinder \
        -f PROTEOMES/ \
        -t !{task.cpus}
    '''
}

process ASSIGN_GENE_FAMILIES_TO_ORTHOGROUPS {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
        val src_julia_1
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}

    echo "Define the location of the results of OrthoFinder run, i.e. the most recent output folder."
    DIR_ORTHOFINDER_OUT=$(ls -tr PROTEOMES/OrthoFinder/ | tail -n1)
    DIR_ORTHOGROUPS=$(pwd)/PROTEOMES/OrthoFinder/${DIR_ORTHOFINDER_OUT}
    
    echo "Define the location of the 15,619 protein family HMMs."
    DIR_PANTHER=$(pwd)/PantherHMM_17.0/famlib/rel/PANTHER17.0_altVersion/hmmscoring/PANTHER17.0/books
    GOT_PATHER=$(pwd)/PantherHMM_17.0/PANTHER17.0_HMM_classifications
    
    echo "Add Orthogroup name to each protein sequence name and merge so that we can be more efficient with hmmsearch."
    echo '#!/bin/bash
    f=$1
    ORTHOGROUP=$(basename ${f} | sed s/.fa$//g)
    sed "s/^>/>${ORTHOGROUP}:/g" ${f} > ${ORTHOGROUP}.tmp
    ' > rename_sequences_with_orthogroup_ID.sh
    chmod +x rename_sequences_with_orthogroup_ID.sh

    echo "Split the large list of orthogroup filenames (too long for bash commands)."
    find ${DIR_ORTHOGROUPS}/Orthogroup_Sequences/ -name '*.fa' > orthogroup_filenames.tmp
    split -l 10000 orthogroup_filenames.tmp orthogroup_filenames-SPLIT-

    echo "For each chunk of orthogroup filenames append the orthogroup names into the protein sequence names in parallel."
    time \
    for F in $(ls orthogroup_filenames-SPLIT-*)
    do
        parallel -j !{task.cpus}\
            ./rename_sequences_with_orthogroup_ID.sh \
            {} \
            ::: $(cat $F)
    done
    rm rename_sequences_with_orthogroup_ID.sh

    echo "Merge orthogroup sequences into a single fasta file."
    MERGED_ORTHOGROUPS=$(pwd)/PROTEOMES/orthogroups.faa
    touch $MERGED_ORTHOGROUPS
    time \
    for f in $(ls | grep "^OG" | grep "tmp$")
    do
        cat $f >> ${MERGED_ORTHOGROUPS}
        rm $f
    done
    echo "Cleanup"
    rm orthogroup_filenames*

    echo "Find PantherHMM protein families for each orthogroup."
    echo '#!/bin/bash
    PROTFA=$1
    DIR_PANTHER=$2
    d=$3
    HMMODL=${DIR_PANTHER}/${d}/hmmer.hmm
    OUTEMP=${PROTFA}-hhmer_gene_family_hits-${d}.tmp
    # hmmsearch -E 0.0001 --tblout ${OUTEMP} ${HMMODL} ${PROTFA}
    sed "/^#/d" ${OUTEMP} | awk @{print $1,$3,$5}@ > ${OUTEMP}.tmp
    if [ $(cat ${OUTEMP}.tmp | wc -l) -eq 0 ]
    then
        rm ${OUTEMP} ${OUTEMP}.tmp
    else
        mv ${OUTEMP}.tmp ${OUTEMP}
    fi
    ' | sed "s/@/'/g" > hmmsearch_for_parallel_execution.sh
    chmod +x hmmsearch_for_parallel_execution.sh
    parallel \
    ./hmmsearch_for_parallel_execution.sh \
        ${MERGED_ORTHOGROUPS} \
        ${DIR_PANTHER} \
        {} \
        ::: $(ls $DIR_PANTHER)
    rm hmmsearch_for_parallel_execution.sh

    echo "Concatenate hmmsearch output for each protein family into a single output file."
    PANTHER_ORTHOGROUPS=$(echo ${MERGED_ORTHOGROUPS} | sed s/.faa$//g).pthr
    cat ${MERGED_ORTHOGROUPS}-hhmer_gene_family_hits-* > ${PANTHER_ORTHOGROUPS}
    rm ${MERGED_ORTHOGROUPS}-hhmer_gene_family_hits-*
    
    echo "Find the best fitting gene family to each unique sequence per orthogroup. This means that each orthogroup can have multiple gene families. Next, add family name and GO terms to each gene family."
    grep "^>" ${MERGED_ORTHOGROUPS} | cut -d':' -f1 | sed 's/>//g' | sort | uniq > all_orthogroups.tmp
    julia !{src_julia_1} \
            PROTEOMES/orthogroups.pthr \
            PantherHMM_17.0/Panther17.0_HMM_familyIDs.txt \
            all_orthogroups.tmp \
            ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.GeneCount.tsv \
            ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups_UnassignedGenes.tsv \
            PROTEOMES/orthogroups_gene_counts_families_go.out
    '''
}

process ASSESS_ORTHOGROUPS_DISTRIBUTIONS {
    label "LOW_MEM_LOW_CPU"
    input:
        val dir
        val src_julia_2
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    echo "Preliminary assessment of the distribution of the genes, orthogroups and gene family classifications."
    julia !{src_julia_2} \
        PROTEOMES/orthogroups_gene_counts_families_go.out \
        PROTEOMES/orthogroups_summarised_gene_counts.csv
    '''
}

workflow {
    FIND_ORTHOGROUPS(params.dir)
    ASSIGN_GENE_FAMILIES_TO_ORTHOGROUPS(params.dir, params.src_julia_1)
    ASSESS_ORTHOGROUPS_DISTRIBUTIONS(params.dir, params.src_julia_2)
}
