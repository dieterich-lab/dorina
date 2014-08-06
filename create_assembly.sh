#!/bin/bash

GTF_FILE=${1}

ASSEMBLY=${GTF_FILE/\.gtf/}

mkdir ${ASSEMBLY}
gtf2gff3 ${GTF_FILE} > ${ASSEMBLY}/${ASSEMBLY}.gff

pushd ${ASSEMBLY} > /dev/null
grep gene ${ASSEMBLY}.gff > all.gff
grep CDS ${ASSEMBLY}.gff > cds.gff
grep three_prime ${ASSEMBLY}.gff > 3_utr.gff
grep five_prime ${ASSEMBLY}.gff > 5_utr.gff
grep exon ${ASSEMBLY}.gff > exon.gff

bedtools subtract -s -a all.gff -b exon.gff |sed -e "s/\tgene\t/\tintron\t/" > intron.gff

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from ${ASSEMBLY}.chromInfo" > ${ASSEMBLY}.genome

bedtools complement -i all.gff -g ${ASSEMBLY}.genome > intergenic.bed

rm exon.gff
popd > /dev/null

