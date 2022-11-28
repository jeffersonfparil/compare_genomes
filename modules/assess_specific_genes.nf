////////////////////////////
// ASSESS SPECIFIC GENES  //
////////////////////////////
// Assumed extension names:
//  - genome: *.fna
//  - annotation: *.gff
//  - coding DNA: *.cds
//  - proteome: *.faa

process DOWNLOAD_GENE_AASEQ {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
        val genes
        val cafe5_n_gamma_cats
        val cafe5_pvalue
    output:
        val dir
        val cafe5_n_gamma_cats
        val cafe5_pvalue
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}/SPECIFIC_GENES
    for g in $(cut -d, -f1 !{genes} | uniq)
    do
        grep "${g}" !{genes} | cut -d, -f3 > ${g}.tmp
        parallel -j !{task.cpus} \
            wget {} ::: $(cat ${g}.tmp)
        cat *.fasta > ${g}.faa
        rm *.fasta *.tmp
    done
    '''
}

process IDENTIFY_GENES_FROM_ORTHOGROUPS {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
        val cafe5_n_gamma_cats
        val cafe5_pvalue
    output:
        val dir
        val cafe5_n_gamma_cats
        val cafe5_pvalue
    shell:
    '''
    #!/usr/bin/env bash
    echo "Define the location of the results of OrthoFinder run, i.e. the most recent output folder."
    cd !{dir}/
    DIR_ORTHOFINDER_OUT=$(ls -tr PROTEOMES/OrthoFinder/ | tail -n1)
    DIR_ORTHOGROUPS=$(pwd)/PROTEOMES/OrthoFinder/${DIR_ORTHOFINDER_OUT}
    DIR_ORTHOGROUP_SEQS=${DIR_ORTHOGROUPS}/Orthogroup_Sequences

    echo "Generate BLAST database for each orthogroup (Outputs: {ORTHOGROUP}.*)"
    echo '#!/bin/bash
    f=$1
    # f=$(find $DIR_ORTHOGROUP_SEQS -name "*fa" | head -n1)
    makeblastdb \
        -in $f \
        -dbtype prot
    ' > prepare_blastdb_per_orthogroup.sh
    chmod +x prepare_blastdb_per_orthogroup.sh
    find $DIR_ORTHOGROUP_SEQS -name "*.fa" > orthogroup_faa_list.tmp
    split --additional-suffix=".tmp" --lines 10000 orthogroup_faa_list.tmp
    for f in $(ls x*.tmp)
    do
        parallel  -j !{task.cpus} \
            ./prepare_blastdb_per_orthogroup.sh \
                {} ::: $(cat ${f})
    done

    echo "Blastp (Outputs: {GENE}-{ORTHOGROUP}.blastout)"
    echo '#!/bin/bash
    query=$1
    dbase=$2
    ortname=$(basename $dbase)
    out=${query}-${ortname%.fa*}.blastout
    blastp \
        -db ${dbase} \
        -query ${query} \
        -out  ${out} \
        -evalue 1e-10 \
        -max_hsps 1 \
        -outfmt "6 qseqid sseqid slen qlen length evalue bitscore qcovs qcovhsp"
    if [ $(cat ${out} | wc -l) -eq 0 ]
    then
        rm ${out}
    fi
    ' > blastp_and_remove_no_hits.sh
    chmod +x blastp_and_remove_no_hits.sh
    for f in $(ls x*.tmp)
    do
        parallel  -j !{task.cpus} \
            ./blastp_and_remove_no_hits.sh \
            {1} \
            {2} \
            ::: $(ls SPECIFIC_GENES/*.faa) \
            ::: $(cat ${f})
    done

    echo "Cleanup"
    rm *.tmp
    rm prepare_blastdb_per_orthogroup.sh
    rm blastp_and_remove_no_hits.sh

    echo "Output:"
    echo "  (1/1) {query}-{ortname}.blastout"
    '''
}

process GENE_EXPANSION_CONTRACTION {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
        val cafe5_n_gamma_cats
        val cafe5_pvalue
    output:
        val dir
    shell:
    '''
    #!/usr/bin/env bash
    echo "Define the location of the results of OrthoFinder run, i.e. the most recent output folder."
    cd !{dir}
    DIR_ORTHOFINDER_OUT=$(ls -tr PROTEOMES/OrthoFinder/ | tail -n1)
    DIR_ORTHOGROUPS=$(pwd)/PROTEOMES/OrthoFinder/${DIR_ORTHOFINDER_OUT}
    ORTHOUT=$(pwd)/PROTEOMES/orthogroups_gene_counts_families_go.out
    TREE=${DIR_ORTHOGROUPS}/Species_Tree/SpeciesTree_rooted.txt

    echo "Extract orthologs per gene (Outputs: {GENE}.ortho)"
    cd SPECIFIC_GENES/
    echo 'args = commandArgs(trailingOnly=TRUE)
    # args = c("EPSPS")
    gene = args[1]
    fnames = list.files()
    fnames = fnames[grep(paste0(gene, "-"), fnames)]
    fnames = fnames[grep("blastout$", fnames)]
    orthoout = c()
    for (f in fnames){
        # f = fnames[1]
        orthogroup = unlist(strsplit(unlist(strsplit(f, "-"))[2], "[.]"))[1]
        df = read.delim(f, header=FALSE)
        if (mean(df[,ncol(df)]) >= 50){
            orthoout = c(orthoout, orthogroup)
        }
    }
    out = data.frame(gene=rep(gene, times=length(orthoout)), orthogroup=orthoout)
    write.table(out, file=paste0(gene, ".ortho"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\\t")
    ' > extract_orthogroups_per_gene.R
    parallel -j !{task.cpus} \
        Rscript extract_orthogroups_per_gene.R \
            {} ::: $(ls *.blastout | cut -d"-" -f1 | sort | uniq)
    echo "Move blast output into a seperate directory to clean-up a bit."
    mkdir BLASTOUT/
    mv *.blastout BLASTOUT/

    echo "Extract TSR and NTSR gene families (Output: {GENE}.orthocounts)"
    rev ${ORTHOUT} | cut -f5- | rev > coln.tmp
    awk -F'\\t' '{print $(NF-1)}' ${ORTHOUT} > col1.tmp
    paste -d'\\t' col1.tmp coln.tmp > counts.tmp
    echo "Infer gene family expansion and contraction (Outpus: {GENE}.conex)"
    for GENE in $(ls *.ortho | sed 's/.ortho//g')
    do
        # GENE=$(ls *.ortho | sed 's/.ortho//g' | head -n1)
        ### Extract gene family counts
        cat ${GENE}.ortho | cut -f2 > ${GENE}.ortho.tmp
        head -n1 counts.tmp > ${GENE}.orthocounts
        grep -f ${GENE}.ortho.tmp counts.tmp >> ${GENE}.orthocounts
        ### Run with lambda_i ~ Gamma(alpha), for each of the i_th gene family category (Output: ${GENE}_CAFE_Gamma100_results/)
        cafe5 \
            --infile ${GENE}.orthocounts \
            --tree ${TREE} \
            --n_gamma_cats !{cafe5_n_gamma_cats} \
            --cores !{task.cpus} \
            --pvalue !{cafe5_pvalue} \
            --output_prefix ${GENE}_CAFE_results
        ### Output(Output: ${GENE}.conex)
        echo -e "Species\\tExpansion\\tContraction" > ${GENE}.conex
        grep -v "^#" ${GENE}_CAFE_results/Gamma_clade_results.txt | \
            grep -v "^<" | \
            sed 's/<..>//g' | \
            sed 's/<.>//g' >> ${GENE}.conex
    done

    echo "Clean-up."
    rm *.tmp

    echo "Output:"
    echo "  (1/2) {GENE}.ortho"
    echo "  (2/2) {GENE}.conex"
    '''
}

// process EXTRACT_GENE_CDS {
//     label "HIGH_MEM_HIGH_CPU"
//     input:
//         val dir
//     output:
//         val 0
//     shell:
//     '''
//     #!/usr/bin/env bash
//     echo "Define the location of the results of OrthoFinder run, i.e. the most recent output folder."
//     cd !{dir}
//     DIR_ORTHOFINDER_OUT=$(ls -tr PROTEOMES/OrthoFinder/ | tail -n1)
//     DIR_ORTHOGROUPS=$(pwd)/PROTEOMES/OrthoFinder/${DIR_ORTHOFINDER_OUT}
//     ORTHOUT=$(pwd)/PROTEOMES/orthogroups_gene_counts_families_go.out

    
//     1. Extract EPSPS CDS (i.e. all orthologs and paralogs within blast-hit orthologs) (Outputs: ${species}-${gene}-${ortho}.cds)
//     ### Extract species names and number of species
//     head -n1 ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv | cut -f2- > species_names.tmp
//     NSPECIES=$(awk -F"\t" "{print NF}" species_names.tmp)

//     ### Extract gene names
//     echo '#!/bin/bash
//     f=$1
//     DIR_ORTHOGROUPS=$2
//     NSPECIES=$3
//     # f=$(ls *.ortho | head -n13 | tail -n1)
//     gene=$(echo ${f%.ortho*})
//     for ortho in $(cut -f2 $f)
//     do
//         # ortho=$(cut -f2 $f | head -n1 | tail -n1)
//         grep "$ortho" ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv > \
//             ${gene}-${ortho}-list_gene_names.tmp
//         for i in $(seq 1 $NSPECIES)
//         do
//             # i=1
//             species=$(cut -f$i species_names.tmp | sed -z "s/\r//g" | sed -z "s/\n//g")
//             col=$(echo $i + 1 | bc)
//             cut -f${col} ${gene}-${ortho}-list_gene_names.tmp | \
//                 sed -z "s/, /\n/g" | \
//                 sed "s/$species|//g" | \
//                 sed "/^$/d" | \
//                 sed "/\r/d" > \
//                 ${species}-${gene}-${ortho}-list_gene_names.tmp
//         if [ $(cat ${species}-${gene}-${ortho}-list_gene_names.tmp | wc -l) -lt 1 ]
//         then
//             rm ${species}-${gene}-${ortho}-list_gene_names.tmp
//         fi
//         done
//         rm  ${gene}-${ortho}-list_gene_names.tmp
//     done
//     ' > extract_gene_names.sh
//     chmod +x extract_gene_names.sh
//     # time \
//     # parallel ./extract_gene_names.sh \
//     #     {} \
//     #     ${DIR_ORTHOGROUPS} \
//     #     ${NSPECIES} \
//     #     ::: $(ls *.ortho)
//     time ./extract_gene_names.sh \
//         EPSPS.ortho \
//         ${DIR_ORTHOGROUPS} \
//         ${NSPECIES}

//     ### Extract gene sequences
//     echo '#!/bin/bash
//     f=$1
//     DIR=$2
//     SRC=$3
//     species=$(echo $f | cut -d"-" -f1)
//     gene=$(echo $f | cut -d"-" -f2)
//     ortho=$(echo $f | cut -d"-" -f3)
//     for query in $(cat $f)
//     do
//         julia ${SRC}/extract_sequence_using_name_query.jl \
//                         ${DIR}/${species}.cds \
//                         ${query} \
//                         ${species}-${gene}-${ortho}-${query}.cds.tmp \
//                         ${species}-${gene}-${ortho}-${query} \
//                         false
//     done
//     cat ${species}-${gene}-${ortho}-*.cds.tmp > ${species}-${gene}-${ortho}.cds
//     rm $f
//     ' > extract_sequences_in_parallel.sh
//     chmod +x extract_sequences_in_parallel.sh
//     # time \
//     # parallel ./extract_sequences_in_parallel.sh \
//     #     {} \
//     #     ${DIR} \
//     #     ${SRC} \
//     #     ::: $(ls *-list_gene_names.tmp)
//     time \
//     parallel ./extract_sequences_in_parallel.sh \
//         {} \
//         ${DIR} \
//         ${SRC} \
//         ::: $(ls *-EPSPS-*-list_gene_names.tmp)

//     '''
// }

workflow {
    // DOWNLOAD_GENE_AASEQ(params.dir, params.genes) | \
        // IDENTIFY_GENES_FROM_ORTHOGROUPS | \
        GENE_EXPANSION_CONTRACTION(params.dir,
                                   params.cafe5_n_gamma_cats,
                                   params.cafe5_pvalue)
}
