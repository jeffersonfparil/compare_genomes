///////////////////////////////////////////////////////////////////////////////////////
// TRANSVERSION RATES AMONG 4-FOLD DEGENRATE SITES USING SINGLE-COPY GENE ORTHOLOGS  //
///////////////////////////////////////////////////////////////////////////////////////
// Assumed extension names:
//  - genome: *.fna
//  - annotation: *.gff
//  - coding DNA: *.cds
//  - proteome: *.faa

process CALCULATE_4DTV {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
        val src_julia_6
    output:
        val dir
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    TYPE=NT
    ### Compute the transversion rate among 4-fold degenerate sites (Output: ${ORTHOLOG}.${TYPE}.cds.4DTv.tmp)
    ### Calculate 4DTv: ratio of transversions in 4-fold degenerate sites
    ### 4-fold degenerate codons:
    ### (1) Ala - GCN
    ### (2) Arg - CGN (etc)
    ### (3) Gly - GGN
    ### (4) Leu - CTN (etc)
    ### (5) Pro - CCN
    ### (6) Ser - TCN (etc)
    ### (7) Thr - ACN
    ### (8) Val - GTN
    ## Transversions:
    ## A <-> C
    ## A <-> T
    ## G <-> C
    ## G <-> T
    parallel -j !{task.cpus} \
        julia !{src_julia_6} \
            {1} \
            {1}.4DTv.tmp \
    ::: $(ls *.${TYPE}.cds)
    ### Concatenate pairwise 4DTv among species across single-copy orthogroups
    echo -e "ORTHOGROUP\\tSPECIES_1\\tSPECIES_2\\tn4D_sites\\tnTv4D_sites\\t4DTv" > ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.4DTv
    for f in $(ls *.cds.4DTv.tmp)
    do
        # f=$(ls *.cds.4DTv.tmp | head -n10 | tail -n1)
        n=$(cat $f | wc -l)
        printf "${f%.${TYPE}*}\n%.0s" $(seq 1 ${n}) > col1.tmp
        sed -z "s/ /\\t/g" $f | sed -z "s/:/\\t/g"  > col2n.tmp
        paste col1.tmp col2n.tmp >> ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.4DTv
    done
    # rm *.${TYPE}.cds
    # rm *.tmp
    '''
}

workflow {
    CALCULATE_4DTV(params.dir, params.src_julia_6)
}
