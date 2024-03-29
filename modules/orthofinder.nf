/////////////////////////////////
// FIND ORTHOLOGS AND PARALOGS //
/////////////////////////////////
// Assumed extension names:
//  - genome: *.fna
//  - annotation: *.gff
//  - coding DNA: *.cds
//  - proteome: *.faa

//NOTE: 20221207 BUG WITH ORTHOFINDER MAKING TREES!!!!! UUUGGGGHHHHHHHHHHHHH

process FIND_ORTHOGROUPS {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
    output:
        val dir
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    echo "Append species names into protein sequence names."
    for f in $(ls PROTEOMES/*.faa)
    do
        fname=$(basename ${f})
        species=${fname%.faa*}
        sed "s/^>/>$species|/g" $f > ORTHOGROUPS/${fname}
    done
    
    echo "Remove gaps (i.e. '.' and '-') in amino acid sequences"
    # for f in ORTHOGROUPS/*.faa
    for f in ORTHOGROUPS/*.faa
    do
        sed -i -e '/^>/!s/[.]//g' $f
        sed -i -e '/^>/!s/-//g' $f
    done
    
    echo "Run OrthoFinder."
    orthofinder \
        -f ORTHOGROUPS/ \
        -t !{task.cpus}

    echo "Define the location of the results of OrthoFinder run, i.e. the most recent output folder."
    DIR_ORTHOFINDER_OUT=$(ls -tr ORTHOGROUPS/OrthoFinder/ | tail -n1)
    DIR_ORTHOGROUPS=$(pwd)/ORTHOGROUPS/OrthoFinder/${DIR_ORTHOFINDER_OUT}
    TREE=${DIR_ORTHOGROUPS}/Species_Tree/SpeciesTree_rooted.txt ### May not be generated successfully due to branch conflicts/inconsistencies

    echo "Output:"
    echo "  (1/1) ORTHOGROUPS/OrhoFinder/Results_{Mmmdd}"
    '''
}

process ASSIGN_GENE_FAMILIES_TO_ORTHOGROUPS {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
    output:
        val dir
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}

    echo "Define the location of the results of OrthoFinder run, i.e. the most recent output folder."
    DIR_ORTHOFINDER_OUT=$(ls -tr ORTHOGROUPS/OrthoFinder/ | tail -n1)
    DIR_ORTHOGROUPS=$(pwd)/ORTHOGROUPS/OrthoFinder/${DIR_ORTHOFINDER_OUT}
    
    echo "Define the location of the 15,619 protein family HMMs."
    DIR_PANTHER=$(pwd)/PantherHMM/famlib/rel/PANTHER*_altVersion/hmmscoring/PANTHER*/books
    GOT_PATHER=$(pwd)/PantherHMM/PANTHER_HMM_classifications
    
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
        parallel -j !{task.cpus} \
            ./rename_sequences_with_orthogroup_ID.sh \
            {} \
            ::: $(cat $F)
    done
    rm rename_sequences_with_orthogroup_ID.sh

    echo "Merge orthogroup sequences into a single fasta file."
    MERGED_ORTHOGROUPS=$(pwd)/ORTHOGROUPS/orthogroups.faa
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
    hmmsearch -E 0.0001 --tblout ${OUTEMP} ${HMMODL} ${PROTFA}
    sed "/^#/d" ${OUTEMP} | awk @{print $1,$3,$5}@ > ${OUTEMP}.tmp
    if [ $(cat ${OUTEMP}.tmp | wc -l) -eq 0 ]
    then
        rm ${OUTEMP} ${OUTEMP}.tmp
    else
        mv ${OUTEMP}.tmp ${OUTEMP}
    fi
    ' | sed "s/@/'/g" > hmmsearch_for_parallel_execution.sh
    chmod +x hmmsearch_for_parallel_execution.sh
    parallel -j !{task.cpus} \
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
    julia !{projectDir}/../scripts/orthogroup_classification_gene_family_GO_terms.jl \
            ORTHOGROUPS/orthogroups.pthr \
            PantherHMM/Panther_HMM_familyIDs.txt \
            all_orthogroups.tmp \
            ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.GeneCount.tsv \
            ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups_UnassignedGenes.tsv \
            ORTHOGROUPS/orthogroups_gene_counts_families_go.out

    echo "Output:"
    echo "  (1/1) ORTHOGROUPS/orthogroups_gene_counts_families_go.out"
    '''
}

process ASSESS_ORTHOGROUPS_DISTRIBUTIONS {
    label "LOW_MEM_LOW_CPU"
    input:
        val dir
    output:
        0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    echo "Preliminary assessment of the distribution of the genes, orthogroups and gene family classifications."
    julia !{projectDir}/../scripts/count_genes_per_ortholog_paralog_classes.jl \
        ORTHOGROUPS/orthogroups_gene_counts_families_go.out \
        ORTHOGROUPS/orthogroups_summarised_gene_counts.csv

    echo "Output:"
    echo "  (1/1) ORTHOGROUPS/orthogroups_summarised_gene_counts.csv"
    '''
}

workflow {
    FIND_ORTHOGROUPS(params.dir) | \
        ASSIGN_GENE_FAMILIES_TO_ORTHOGROUPS | \
        ASSESS_ORTHOGROUPS_DISTRIBUTIONS
}
