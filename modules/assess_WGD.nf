//////////////////////////////////////
// WHOLE GENOME DUPLICATION EVENTS  //
//////////////////////////////////////
// Assumed extension names:
//  - genome: *.fna
//  - annotation: *.gff
//  - coding DNA: *.cds
//  - proteome: *.faa

process ASSESS_WGD {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
        val src_shell_1
        val src_julia_4
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    echo "Define the location of the results of OrthoFinder run, i.e. the most recent output folder."
    cd !{dir}
    DIR_ORTHOFINDER_OUT=$(ls -tr PROTEOMES/OrthoFinder/ | tail -n1)
    DIR_ORTHOGROUPS=$(pwd)/PROTEOMES/OrthoFinder/${DIR_ORTHOFINDER_OUT}
    ORTHOUT=PROTEOMES/orthogroups_gene_counts_families_go.out
    echo "Identify multi-copy paralogs (2 to 5 copies) per species, align, and estimate 4DTv"
    head -n1 ${ORTHOUT} | rev | cut -f5- | rev | cut -f2- | sed -z "s/\\t/\\n/g" > species_names.tmp
    # for i in $(seq 1 $(cat species_names.tmp | wc -l))
    # do
       i=1
       ### Extract species name
       SPECIES=$(head -n${i} species_names.tmp | tail -n1)
       idx=$(echo $i + 1 | bc)
       ### Exract names of orthogroups with 2 copies in the current species
       awk -v col="${idx}" '($col >= 2) && ($col <= 5)' $ORTHOUT | cut -f1 > multi_gene_list.grep
       ### Extract names of the genes of these multi-copy orthogroups
       grep -f multi_gene_list.grep ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv | cut -f1,${idx} > multi_gene_list.geneNames
       ### Extract CDS, and align in parallel
       parallel -j !{task.cpus} \
       bash !{src_shell_1} \
           {} \
           !{src_julia_4} \
           ::: $(seq 1 $(cat multi_gene_list.geneNames | wc -l))
       ### Concatenate 4DTv estimates
       cat *.4DTv.tmp > ${SPECIES}.4DTv
       ### Clean-up
       # rm *.4DTv.tmp
       # rm multi_gene_list.grep
       # rm multi_gene_list.geneNames
    # done
    ### Clean-up
    # rm species_names.tmp
    '''
}

workflow {
    ASSESS_WGD(params.dir,
               params.src_shell_1,
               params.src_julia_4)
}
