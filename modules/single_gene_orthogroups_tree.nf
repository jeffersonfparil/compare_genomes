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
        val src_julia_3
        val src_julia_4
        val src_julia_5
        val dates
    output:
        val dir
        val src_julia_4
        val src_julia_5
        val dates
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}

    echo "Define the location of the results of OrthoFinder run, i.e. the most recent output folder."
    DIR_ORTHOFINDER_OUT=$(ls -tr PROTEOMES/OrthoFinder/ | tail -n1)
    DIR_ORTHOGROUPS=$(pwd)/PROTEOMES/OrthoFinder/${DIR_ORTHOFINDER_OUT}

    echo "Find the single-gene orthogroups for all species"
    ORTHOUT=PROTEOMES/orthogroups_gene_counts_families_go.out
    NSPECIES=$(ls PROTEOMES/*.faa | grep -v "orthogroups.faa" | wc -l)
    julia !{src_julia_3} \
        ${ORTHOUT} \
        ${NSPECIES} ### output: single_gene_list.grep
    grep -f single_gene_list.grep ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv > single_gene_list.geneNames

    echo '#!/bin/bash
    i=$1
    line=$(head -n${i} single_gene_list.geneNames | tail -n1)
    ORTHONAME=$(echo $line | cut -d" " -f1)
    for name in $(echo $line | cut -d" " -f2-)
    do
        SPECIES=$(echo $name | cut -d"|" -f1)
        GENE_NAME=$(echo $name | cut -d"|" -f2)
        julia !{src_julia_4} \
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
    ./parallel_extract_single_gene_orthogroups.sh {} ::: $(seq 1 $(cat single_gene_list.geneNames | wc -l))
    ### Cleanup
    rm single_gene_list.grep 
    rm parallel_extract_single_gene_orthogroups.sh
    '''
}

process ALIGN_SINGLE_GENE_ORTHOGROUPS {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
        val src_julia_4
        val src_julia_5
        val dates
    output:
        val dir
        val src_julia_4
        val src_julia_5
        val dates
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    
    echo '#!/bin/bash
    f=$1
    MACSE=$2
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
    ./parallel_align_cds.sh {} ${MACSE} \
    ::: $(ls OG*.fasta)
    ### Cleanup
    rm parallel_align_cds.sh
    '''
}

process BUILD_TREE {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
        val src_julia_4
        val src_julia_5
        val dates
    output:
        val dir
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    
    TYPE=NT.cds
    ### Extract sequences per species (Outputs: ${ORTHONAME}-${SPECIES}.fasta)
    parallel -j !{task.cpus} \
    julia !{src_julia_4} \
        {1}.${TYPE} \
        {2} \
        {1}-{2}.fasta \
        {1}-{2} \
        false \
    ::: $(ls *.${TYPE} | sed "s/.$TYPE//g") \
    ::: $(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g')

    ### Concatenate alignments per species (Outputs: ${SPECIES}.aln)
    for SPECIES in $(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g')
    do
        echo $SPECIES
        # Concatenate sequences
        echo ">${SPECIES}" > ${SPECIES}.aln.tmp
        for f in $(ls *-${SPECIES}.fasta)
        do
            sed '/^>/d' $f | sed -z 's/\\n//g' >> ${SPECIES}.aln.tmp
        done
        echo "" >> ${SPECIES}.aln.tmp
        julia !{src_julia_5} \
            ${SPECIES}.aln.tmp \
            50 \
            ${SPECIES}-${TYPE%.*}.aln
        rm ${SPECIES}.aln.tmp
    done

    ### Extract sequence lengths to build the sequence partitioning nexus file (Output: alignment_parition.${TYPE%.*}.nex)
    SPECIES=$(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g' | head -n1)
    echo '#nexus
    begin sets;' > alignment_parition.${TYPE%.*}.nex
    N0=0
    for f in $(ls *-${SPECIES}.fasta)
    do
        # f=$(ls *-${SPECIES}.fasta | head -n1)
        NAME=$(head -n1 $f | sed 's/>//g' | cut -d'-' -f1)
        N1=$(cat $f | sed '/^>/d' | sed -z 's/\\n//g' | wc -c)
        START=$(echo "$N0 + 1" | bc)
        END=$(echo "$N0 + $N1" | bc)
        echo "charset $NAME = $START-$END;" >> alignment_parition.${TYPE%.*}.nex
        N0=$END
    done
    echo 'end;' >> alignment_parition.${TYPE%.*}.nex

    ### Concatenate species alignments (Output: ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln)
    cat *-${TYPE%.*}.aln > ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln.tmp
    mv ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln.tmp ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln

    ### Build tree
    BOOTSTRAP_REPS=1000
    THREADS=!{task.cpus}
    TIP_DATE=0
    iqtree2 \
        -s ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln \
        -p alignment_parition.${TYPE%.*}.nex \
        -B ${BOOTSTRAP_REPS} \
        -T ${THREADS} \
        --date !{dates} \
        --date-tip ${TIP_DATE} \
        --prefix ORTHOGROUPS_SINGLE_GENE.${TYPE%.*} \
        --redo
    '''
}

workflow {
    // IDENTIFY_SINGLE_GENE_ORTHOGROUPS(params.dir, params.src_julia_3, params.src_julia_4, params.src_julia_5, params.dates) | \
    //     ALIGN_SINGLE_GENE_ORTHOGROUPS | \
        BUILD_TREE(params.dir, params.src_julia_4, params.src_julia_5, params.dates)
}


// 2. Extract the CDS of these genes (Outputs: ${ORTHONAME}.fasta [includes sequences from each species]):
// ```shell

// ```

// 3. Align CDS (Outputs: ${ORTHOLOG}.NT.cds [nucleotide alignments] and ${ORTHOLOG}.AA.prot [amino acid alignments])
// ```shell

// ```

// 4. Build the tree (Output: ORTHOGROUPS_SINGLE_GENE.NT.timetree.nex)
// ```shell




// 5. **Additional**: Compute pairwise 4DTv (Output: ORTHOGROUPS_SINGLE_GENE.NT.4DTv)
// ```shell
// ### Compute the transversion rate among 4-fold degenerate sites (Output: ${ORTHOLOG}.NT.cds.4DTv.tmp)
// time \
// parallel \
// julia calculate_4DTv.jl {1} {1}.4DTv.tmp \
//     ::: $(ls *.NT.cds)
// ### Concatenate pairwise 4DTv among species across single-copy orthogroups
// echo -e "ORTHOGROUP\tSPECIES_1\tSPECIES_2\tn4D_sites\tnTv4D_sites\t4DTv" > ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.4DTv
// for f in $(ls *.cds.4DTv.tmp)
// do
//     # f=$(ls *.cds.4DTv.tmp | head -n10 | tail -n1)
//     n=$(cat $f | wc -l)
//     printf "${f%.${TYPE}*}\n%.0s" $(seq 1 $n) > col1.tmp
//     sed -z "s/ /\t/g" $f | sed -z "s/:/\t/g"  > col2n.tmp
//     paste col1.tmp col2n.tmp >> ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.4DTv
// done
// ```

// 6. Clean-up
// ```shell
// rm OG*.fasta
// rm OG*.NT.cds
// rm single_gene_list.*
// rm dates.txt
// rm *-${TYPE%.*}.aln
// rm *.tmp
// ```
