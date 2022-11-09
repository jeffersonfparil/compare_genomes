#!/usr/bin/env bash
for line in $(cat urls.txt)
do
    spe=$(echo ${line} | cut -d, -f1)
    url=$(echo ${line} | cut -d, -f2)
    if [ $(echo ${spe} | rev | cut -d. -f1 | rev) = "fna" ]
    then
        DIR=$(echo /data-weedomics-3/GENOMES)
    elif [ $(echo ${spe} | rev | cut -d. -f1 | rev) = "gff" ]
    then
        DIR=$(echo /data-weedomics-3/GFF)
    elif [ $(echo ${spe} | rev | cut -d. -f1 | rev) = "cds" ]
    then
        DIR=$(echo /data-weedomics-3/CDS)
    else
        DIR=$(echo /data-weedomics-3/PROTEOMES)
    fi
    wget ${url} -O - | gunzip -c - > ${DIR}/${spe}
done
