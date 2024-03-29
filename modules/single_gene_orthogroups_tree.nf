//////////////////////////////////
// SINGLE-GENE ORTHOGROUPS TREE //
//////////////////////////////////
// Assumed extension names:
//  - genome: *.fna
//  - annotation: *.gff
//  - coding DNA: *.cds
//  - proteome: *.faa

process IDENTIFY_SINGLE_GENE_ORTHOGROUPS {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
        val dates
    output:
        val dir
        val dates
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}

    echo "Define the location of the results of OrthoFinder run, i.e. the most recent output folder."
    DIR_ORTHOFINDER_OUT=$(ls -tr ORTHOGROUPS/OrthoFinder/ | tail -n1)
    DIR_ORTHOGROUPS=$(pwd)/ORTHOGROUPS/OrthoFinder/${DIR_ORTHOFINDER_OUT}

    echo "List single-gene orthogroups for all species"
    ORTHOUT=$(pwd)/ORTHOGROUPS/orthogroups_gene_counts_families_go.out
    NSPECIES=$(ls ORTHOGROUPS/*.faa | grep -v "orthogroups.faa" | wc -l)
    julia !{projectDir}/../scripts/extract_single_gene_orthogroups.jl \
        ${ORTHOUT} \
        ${NSPECIES} ### output: single_gene_list.grep
    grep -f single_gene_list.grep ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv > single_gene_list.geneNames

    echo "Find the single-gene orthogroups for all species"
    echo '#!/bin/bash
    i=$1
    line=$(head -n${i} single_gene_list.geneNames | tail -n1)
    ORTHONAME=$(echo $line | cut -d" " -f1)
    for name in $(echo $line | cut -d" " -f2-)
    do
        SPECIES=$(echo $name | cut -d"|" -f1)
        GENE_NAME=$(echo $name | cut -d"|" -f2)
        julia !{projectDir}/../scripts/extract_sequence_using_name_query.jl \
            CDS/${SPECIES}.cds \
            ${GENE_NAME} \
            ${ORTHONAME}-${SPECIES}.fasta \
            ${SPECIES} \
            false
    done
    cat ${ORTHONAME}-*.fasta > ${ORTHONAME}.fasta
    rm ${ORTHONAME}-*.fasta
    ' > parallel_extract_single_gene_orthogroups.sh
    chmod +x parallel_extract_single_gene_orthogroups.sh
    parallel -j !{task.cpus} \
        ./parallel_extract_single_gene_orthogroups.sh \
            {} ::: $(seq 1 $(cat single_gene_list.geneNames | wc -l))
    
    echo "Cleanup"
    rm single_gene_list.grep 
    rm parallel_extract_single_gene_orthogroups.sh
    
    echo "Output:"
    echo "  (1/2) {ORTHONAME}.fasta"
    echo "  (2/2) {ORTHONAME}-{SPECIES}.fasta"
    '''
}

process ALIGN_SINGLE_GENE_ORTHOGROUPS {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
        val dates
    output:
        val dir
        val dates
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    echo '#!/bin/bash
    f=$1
    ORTHOLOG=${f%.fasta*}
    # Align the CDS across species
    macse \
        -prog alignSequences \
        -seq ${f} \
        -out_NT ${ORTHOLOG}.aligned.unsorted.cds.tmp \
        -out_AA ${ORTHOLOG}.aligned.unsorted.prot.tmp
    # Convert stop codons and frameshifts as "---" for compatibility with downstream tools
    macse \
        -prog exportAlignment \
        -align ${ORTHOLOG}.aligned.unsorted.cds.tmp \
        -codonForFinalStop --- \
        -codonForInternalStop NNN \
        -codonForExternalFS --- \
        -codonForInternalFS --- \
        -out_NT ${ORTHOLOG}.NT.cds \
        -out_AA ${ORTHOLOG}.AA.prot
    # Clean-up
    rm ${ORTHOLOG}*.tmp ${ORTHOLOG}.AA.prot # we are not using the amino acid sequences
    ' > parallel_align_cds.sh
    chmod +x parallel_align_cds.sh
    time \
    parallel -j !{task.cpus} \
    ./parallel_align_cds.sh \
        {} \
        ::: $(ls OG*.fasta)
    
    echo "Cleanup"
    rm parallel_align_cds.sh
    
    echo "Output:"
    echo "  (1/1) {ORTHOLOG}.NT.cds"
    '''
}

process BUILD_TREE {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
        val dates
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    
    echo "Extract sequences per species (Outputs: {ORTHONAME}-{SPECIES}.fasta)"
    TYPE=NT.cds
    parallel -j !{task.cpus} \
    julia !{projectDir}/../scripts/extract_sequence_using_name_query.jl \
        {1}.${TYPE} \
        {2} \
        {1}-{2}.fasta \
        {1}-{2} \
        false \
    ::: $(ls *.${TYPE} | sed "s/.$TYPE//g") \
    ::: $(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g')

    echo "Remove alignments which do not match lengths across species (i.e. in cases where there are actually multiple alignmets of an orthogroup probably due to mislabelling in the input coding sequences - sequences with the same names)"
    for ORTHO_EXPECTED in $(ls OG*-*_*.fasta | cut -d- -f1 | sort | uniq)
    do
        # ORTHO_EXPECTED=$(ls OG*-*_*.fasta | cut -d- -f1 | sort | uniq | head -n1)
        SPECIES_1=$(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g' | head -n1)
        N_1=$(grep -v "^>" ${ORTHO_EXPECTED}-${SPECIES_1}.fasta | sed -z 's/\\n//g' | wc -c)
        for SPECIES_2 in $(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g' | tail -n+2)
        do
            # SPECIES_2=$(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g' | tail -n+2 | head -n1)
            N_2=$(grep -v "^>" ${ORTHO_EXPECTED}-${SPECIES_2}.fasta | sed -z 's/\\n//g' | wc -c)
            if [ $N_1 -ne $N_2 ]
            then
                echo ${ORTHO_EXPECTED}-${SPECIES_1}.fasta AND ${ORTHO_EXPECTED}-${SPECIES_2}.fasta
                rm ${ORTHO_EXPECTED}-*_*.fasta
                break
            fi
        done
    done

    echo "Concatenate alignments per species (Outputs: {SPECIES}.aln)"
    for SPECIES in $(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g')
    do
        echo $SPECIES
        # Concatenate sequences
        echo ">${SPECIES}" > ${SPECIES}.aln.tmp
        for f in $(ls *-${SPECIES}.fasta)
        do
            sed '/^>/d' $f | sed -z 's/\\n//g' >> ${SPECIES}.aln.tmp  ### Note that you need to remove one of the escape characters in the newline character, if you need to run this on bash manually
        done
        echo "" >> ${SPECIES}.aln.tmp
        julia !{projectDir}/../scripts/reformat_fasta_sequence.jl \
            ${SPECIES}.aln.tmp \
            50 \
            ${SPECIES}-${TYPE%.*}.aln
        rm ${SPECIES}.aln.tmp
    done

    echo "Extract sequence lengths to build the sequence partitioning nexus file (Output: alignment_parition.NT.nex)"
    SPECIES=$(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g' | head -n1)
    echo '#nexus
    begin sets;' > alignment_parition.${TYPE%.*}.nex
    N0=0
    for f in $(ls *-${SPECIES}.fasta)
    do
        # f=$(ls *-${SPECIES}.fasta | head -n1)
        NAME=$(head -n1 $f | sed 's/>//g' | cut -d'-' -f1)
        N1=$(cat $f | sed '/^>/d' | sed -z 's/\\n//g' | wc -c) ### Note that you need to remove one of the escape characters in the newline character, if you need to run this on bash manually
        START=$(echo "$N0 + 1" | bc)
        END=$(echo "$N0 + $N1" | bc)
        echo "charset $NAME = $START-$END;" >> alignment_parition.${TYPE%.*}.nex
        N0=$END
    done
    echo 'end;' >> alignment_parition.${TYPE%.*}.nex

    echo "Concatenate species alignments (Output: ORTHOGROUPS_SINGLE_GENE.NT.aln)"
    cat *-${TYPE%.*}.aln > ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln.tmp
    mv ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln.tmp ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln

    echo "Build tree"
    BOOTSTRAP_REPS=1000
    THREADS=!{task.cpus}
    TIP_DATE=0
    ### Use bootstrapping if you have more than 4 species
    n_species=$(grep "^>" ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln | wc -l)
    if [ ${n_species} -gt 4 ]
    then
        iqtree2 \
            -s ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln \
            -p alignment_parition.${TYPE%.*}.nex \
            -B ${BOOTSTRAP_REPS} \
            -T ${THREADS} \
            --date !{dates} \
            --date-tip ${TIP_DATE} \
            --prefix ORTHOGROUPS_SINGLE_GENE.${TYPE%.*} \
            --redo
    else    
        iqtree2 \
            -s ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln \
            -p alignment_parition.${TYPE%.*}.nex \
            -T ${THREADS} \
            --date !{dates} \
            --date-tip ${TIP_DATE} \
            --prefix ORTHOGROUPS_SINGLE_GENE.${TYPE%.*} \
            --redo
    fi

    if [ -f "ORTHOGROUPS_SINGLE_GENE.NT.timetree.nex" ]
    then
        echo "The time tree was successfully created!"
    else
        echo "Creating the time tree failed! Please check `dates.txt` for dating inconsistencies and/or remove divergence times at your discretion."
    fi

    echo "Cleanup"
    rm OG*.fasta
    rm single_gene_list.*
    rm *-${TYPE%.*}.aln
    
    echo "Output:"
    echo "  (01/12) ORTHOGROUPS_SINGLE_GENE.NT.best_scheme.nex"
    echo "  (02/12) ORTHOGROUPS_SINGLE_GENE.NT.best_scheme"
    echo "  (03/12) ORTHOGROUPS_SINGLE_GENE.NT.model.gz"
    echo "  (04/12) ORTHOGROUPS_SINGLE_GENE.NT.mldist"
    echo "  (05/12) ORTHOGROUPS_SINGLE_GENE.NT.bionj"
    echo "  (06/12) ORTHOGROUPS_SINGLE_GENE.NT.best_model.nex"
    echo "  (07/12) ORTHOGROUPS_SINGLE_GENE.NT.treefile"
    echo "  (08/12) ORTHOGROUPS_SINGLE_GENE.NT.iqtree"
    echo "  (09/12) ORTHOGROUPS_SINGLE_GENE.NT.timetree.nwk"
    echo "  (10/12) ORTHOGROUPS_SINGLE_GENE.NT.timetree.nex"
    echo "  (11/12) ORTHOGROUPS_SINGLE_GENE.NT.timetree.lsd"
    echo "  (12/12) ORTHOGROUPS_SINGLE_GENE.NT.ckp.gz"
    '''
}

workflow {
    IDENTIFY_SINGLE_GENE_ORTHOGROUPS(params.dir, params.dates) | \
        ALIGN_SINGLE_GENE_ORTHOGROUPS | \
        BUILD_TREE
}
