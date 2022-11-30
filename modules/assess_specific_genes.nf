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
        val species_of_interest
        val genes
        val cafe5_n_gamma_cats
        val cafe5_pvalue
    output:
        val dir
        val species_of_interest
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
        val species_of_interest
        val cafe5_n_gamma_cats
        val cafe5_pvalue
    output:
        val dir
        val species_of_interest
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
    makeblastdb \
        -in $f \
        -dbtype prot
    ' > prepare_blastdb_per_orthogroup.sh
    chmod +x prepare_blastdb_per_orthogroup.sh
    find ${DIR_ORTHOGROUP_SEQS} -name "*.fa" > orthogroup_faa_list.tmp
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
    out=${query%.faa*}-${ortname%.fa*}.blastout
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
        val species_of_interest
        val cafe5_n_gamma_cats
        val cafe5_pvalue
    output:
        val dir
        val species_of_interest
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

    echo "Extract gene families (Output: {GENE}.orthocounts)"
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
    echo "  (1/3) {GENE}.ortho"
    echo "  (2/3) {GENE}.conex"
    echo "  (3/3) BLASTOUT/"
    '''
}

process EXTRACT_GENE_CDS {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
        val species_of_interest
    output:
        val dir
        val species_of_interest
    shell:
    '''
    #!/usr/bin/env bash
    echo "Extract species names and number of species"
    cd !{dir}/SPECIFIC_GENES/
    head -n1 ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv | cut -f2- > species_names.tmp
    NSPECIES=$(awk -F"\\t" "{print NF}" species_names.tmp)

    echo "Extract gene names"
    echo '#!/bin/bash
    f=$1
    DIR_ORTHOGROUPS=$2
    NSPECIES=$3
    # f=$(ls *.ortho | head -n13 | tail -n1)
    gene=$(echo ${f%.ortho*})
    for ortho in $(cut -f2 $f)
    do
        # ortho=$(cut -f2 $f | head -n1 | tail -n1)
        grep "$ortho" ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv > ${gene}-${ortho}-list_gene_names.tmp
        for i in $(seq 1 $NSPECIES)
        do
            # i=1
            species=$(cut -f$i species_names.tmp | sed -z "s/\\r//g" | sed -z "s/\\n//g")
            col=$(echo $i + 1 | bc)
            cut -f${col} ${gene}-${ortho}-list_gene_names.tmp | \
                sed -z "s/, /\\n/g" | \
                sed "s/$species|//g" | \
                sed "/^$/d" | \
                sed "/\\r/d" > \
                ${species}-${gene}-${ortho}-list_gene_names.tmp
            if [ $(cat ${species}-${gene}-${ortho}-list_gene_names.tmp | wc -l) -lt 1 ]
            then
                rm ${species}-${gene}-${ortho}-list_gene_names.tmp
            fi
        done
        rm  ${gene}-${ortho}-list_gene_names.tmp
    done
    ' > extract_gene_names.sh
    chmod +x extract_gene_names.sh
    parallel -j !{task.cpus} \
        ./extract_gene_names.sh \
            {} \
            ${DIR_ORTHOGROUPS} \
            ${NSPECIES} \
            ::: $(ls *.ortho)

    echo "Extract gene sequences"
    echo '#!/bin/bash
    f=$1
    DIR=$2
    SRC=$3
    species=$(echo $f | cut -d"-" -f1)
    gene=$(echo $f | cut -d"-" -f2)
    ortho=$(echo $f | cut -d"-" -f3)
    for query in $(cat $f)
    do
        julia ${SRC} \
            ${DIR}/${species}.cds \
            ${query} \
            ${species}-${gene}-${ortho}-${query}.cds.tmp \
            ${species}-${gene}-${ortho}-${query} \
            false
    done
    cat ${species}-${gene}-${ortho}-*.cds.tmp > ${species}-${gene}-${ortho}.cds
    rm $f
    ' > extract_sequences_in_parallel.sh
    chmod +x extract_sequences_in_parallel.sh
    for gene in $(ls *-list_gene_names.tmp | cut -d'-' -f2 | sort | uniq)
    do
        parallel -j !{task.cpus} \
            ./extract_sequences_in_parallel.sh \
                {} \
                !{dir}/CDS \
                !{projectDir}/../scripts/extract_sequence_using_name_query.jl \
                ::: $(ls *-${gene}-*-list_gene_names.tmp)
    done

    echo "Cleanup"
    rm extract_gene_names.sh
    rm extract_sequences_in_parallel.sh
    rm *.tmp

    echo "Output:"
    echo "  (1/1) {species}-{gene}-{ortho}.cds"
    '''
}

process ALIGN_GENE_CDS {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
        val species_of_interest
    output:
        val dir
    shell:
    '''
    #!/usr/bin/env bash
    echo "Merge per gene_per orthogropup prior to alignment (Outputs: {gene}-{ortho}.cds)"
    cd !{dir}/SPECIFIC_GENES/
    for gene in $(ls *.cds | cut -d"-" -f2 | sort | uniq)
    do
        for ortho in $(ls *-${gene}-*.cds | cut -d"-" -f3 | cut -d"." -f1 | sort | uniq)
        do
            fname_output=$(ls *-${gene}-${ortho}.cds | cut -d"-" -f2 | sort | uniq | sed -z "s/\\n/-/g")${ortho}.cds.tmp
            cat *-${gene}-${ortho}.cds > ${fname_output}
            rm *-${gene}-${ortho}.cds
            mv ${fname_output} ${fname_output%.tmp*}
        done
        rm *.tmp
    done
    
    echo "Align CDS per orthogroup per gene (Outputs: {gene}-{ortho}.aln)"
    echo '#!/bin/bash
    f=$1
    ext=$2
    macse \
        -prog alignSequences \
        -seq ${f} \
        -out_NT ${f%.${ext}*}.aligned.unsorted.cds.tmp \
        -out_AA ${f%.${ext}*}.aligned.unsorted.prot.tmp
    # Convert stop codons and frameshifts as "---" for compatibility with downstream tools
    macse \
        -prog exportAlignment \
        -align ${f%.${ext}*}.aligned.unsorted.cds.tmp \
        -codonForFinalStop --- \
        -codonForInternalStop NNN \
        -codonForExternalFS --- \
        -codonForInternalFS --- \
        -out_NT ${f%.${ext}*}.aln \
        -out_AA ${f%.${ext}*}.AA.prot
    ' > align_in_parallel.sh
    chmod +x align_in_parallel.sh
    parallel -j !{task.cpus} \
        ./align_in_parallel.sh \
            {} \
            cds \
            ::: $(ls *.cds)
    rm *.tmp *.AA.prot

    echo "Remove alignments and cds without the genes from the species of interest (Outputs: {gene}-{ortho}.aln)"
    for f in $(ls *.aln)
    do
        # f=$(ls *.aln | head -n1 | tail -n1)
        x=$(grep "^>!{species_of_interest}" $f | wc -l)
        if [ $x -eq 0 ]
        then
            rm $f
            rm ${f%.aln*}.cds
        fi
    done

    echo "Create pairwise cds alignments using the first alignment from the speies of interest as the focal alignment per orthogroup per gene (Outputs: {gene}-{ortho}.aln.pw)"
    echo '#!/bin/bash
    f=$1
    species=$2
    for i in $(seq 1 $(grep "^>${species}" $f | wc -l))
    do
        focal_aln=$(grep "^>${species}" $f | head -n${i} | tail -n1 | sed "s/^>//g")
        grep -A1 "^>${species}" $f | head -n$(echo "$i * 2" | bc) | tail -n1 > ${f}.FOCAL_SEQ.tmp
        touch ${f}-${i}.pw.tmp
        ### temp file without the focal alignment
        sed -e "/$focal_aln/,+1d" $f > ${f}.tmp
        for line in $(seq 2 2 $(cat ${f}.tmp | wc -l))
        do
            # line=2
            curr_aln_name=$(head -n$(echo $line -1 | bc) ${f}.tmp | tail -n1 | sed "s/>//g")
            if [ ${focal_aln} != ${curr_aln_name} ]
            then
                name=${focal_aln}--:--${curr_aln_name}
                echo $name >> ${f}-${i}.pw.tmp                           ### alignment pair name
                cat ${f}.FOCAL_SEQ.tmp >> ${f}-${i}.pw.tmp               ### focal alignment
                head -n${line} ${f}.tmp | tail -n1 >> ${f}-${i}.pw.tmp   ### current alignment
                echo "" >> ${f}-${i}.pw.tmp
            fi
        done
        ### Clean-up
        mv ${f}-${i}.pw.tmp ${f}-${i}.pw
        rm ${f}.FOCAL_SEQ.tmp ${f}.tmp
        ### Remove single alignments
        if [ $(cat ${f}-${i}.pw | wc -l) -eq 0 ]
        then
            rm ${f}-${i}.pw
        fi
    done
    ' > prepare_pairwise_alignments_with_species_of_interest_as_focus_in_parallel.sh
    chmod +x prepare_pairwise_alignments_with_species_of_interest_as_focus_in_parallel.sh
    parallel -j !{task.cpus} \
        ./prepare_pairwise_alignments_with_species_of_interest_as_focus_in_parallel.sh \
            {} \
            !{species_of_interest} \
            ::: $(ls *.aln)

    echo "Cleanup"
    rm align_in_parallel.sh
    rm prepare_pairwise_alignments_with_species_of_interest_as_focus_in_parallel.sh

    echo "Output:"
    echo "  (1/3) {gene}-{ortho}.cds"
    echo "  (2/3) {gene}-{ortho}.aln"
    echo "  (3/3) {gene}-{ortho}.aln.pw"
    '''
}

process KAKS {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    echo "KaKs_calculator2 for pairwise orthogroup gene comparisons with sliding 15-bp windows"
    cd !{dir}/SPECIFIC_GENES/
    echo '#!/bin/bash
    f=$1
    SRC_1=$2
    SRC_2=$3
    julia ${SRC_1} \
        ${f} \
        15 \
        15 \
        ${f}.windows.tmp

    KaKs_Calculator \
        -m MS \
        -i ${f}.windows.tmp \
        -o ${f%.tmp*}.kaks.tmp

    Rscript ${SRC_2} \
        ${f%.tmp*}.kaks.tmp\
        0.001
    ' > KaKs_per_window_and_plot_in_parallel.sh
    chmod +x KaKs_per_window_and_plot_in_parallel.sh
    parallel -j !{task.cpus} \
        ./KaKs_per_window_and_plot_in_parallel.sh \
            {} \
            !{projectDir}/../scripts/split_alignment_pairs.jl \
            !{projectDir}/../scripts/plot_KaKs_across_windows.R \
            ::: $(ls *.pw)

    echo "Cleanup"
    rm *.tmp

    echo "Output:"
    echo "  (1/2) {gene}-{ortho}.aln-{i}.kaks.svg"
    echo "  (2/2) {gene}-{ortho}.aln-{i}.kaks-SIGNIFICANT_PEAKS.csv"
    '''
}

workflow {
    DOWNLOAD_GENE_AASEQ(params.dir,
                        params.species_of_interest,
                        params.genes,
                        params.cafe5_n_gamma_cats,
                        params.cafe5_pvalue) | \
        IDENTIFY_GENES_FROM_ORTHOGROUPS | \
        GENE_EXPANSION_CONTRACTION | \
        EXTRACT_GENE_CDS | \
        ALIGN_GENE_CDS | \
        KAKS
}
