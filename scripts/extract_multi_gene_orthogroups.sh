#!/bin/bash
### NOTE: The file: multi_gene_list.geneNames contains the gene names of one species. The gene names are in the second column (space-delimited), and each gene is comma-space-delimited
j=$1            ### Line number corresponding to the current gene
src_julia_4=$2  ### Julia script to extract sequences using sequence names
src_julia_6=$3  ### Julia script to calculate 4DTv
# j=1017
line=$(head -n${j} multi_gene_list.geneNames | tail -n1)
ORTHONAME=$(echo $line | cut -d" " -f1)
# Extract CDS of these multi-copy genes
for name in $(echo $line | cut -d" " -f2- | sed -z "s/, /\n/g")
do
    # name=$(echo $line | cut -d" " -f2- | sed -z "s/, /\n/g" | head -n2 | tail -n1)
    SPECIES=$(echo $name | cut -d"|" -f1)
    GENE_NAME=$(echo $name | cut -d"|" -f2 | sed -z "s/\r//g")
    julia ${src_julia_4} \
        CDS/${SPECIES}.cds \
        ${GENE_NAME} \
        ${ORTHONAME}-${GENE_NAME}.fasta \
        ${GENE_NAME} \
        false
done
# Concatenate the sequences
if [ $(ls ${ORTHONAME}-*.fasta | wc -l) -gt 1 ]
then
    cat ${ORTHONAME}-*.fasta > ${ORTHONAME}.fasta
fi
# Align the CDS across species
macse \
    -prog alignSequences \
    -seq ${ORTHONAME}.fasta \
    -out_NT ${ORTHONAME}.aligned.unsorted.cds.tmp \
    -out_AA ${ORTHONAME}.aligned.unsorted.prot.tmp
# Convert stop codons and frameshifts as "---" for compatibility with downstream tools
macse \
    -prog exportAlignment \
    -align ${ORTHONAME}.aligned.unsorted.cds.tmp \
    -codonForFinalStop --- \
    -codonForInternalStop NNN \
    -codonForExternalFS --- \
    -codonForInternalFS --- \
    -out_NT ${ORTHONAME}.NT.cds \
    -out_AA ${ORTHONAME}.AA.prot
# Calculate 4DTv (i.e. the ratio of the number of 4-fold degenerate codons with transversion and the total number of 4-fold degenerate codons)
julia ${src_julia_6} \
    ${ORTHONAME}.NT.cds \
    ${ORTHONAME}.4DTv.tmp
# Clean-up
rm ${ORTHONAME}*.fasta
rm ${ORTHONAME}.aligned.unsorted*.tmp ${ORTHONAME}.NT.cds ${ORTHONAME}.AA.prot
##### rm ${ORTHONAME}.axt.tmp ${ORTHONAME}.kaks.tmp