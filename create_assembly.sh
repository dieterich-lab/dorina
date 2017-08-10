#!/bin/bash

GTF_FILE=${1}
# path to this script from :
#https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ASSEMBLY=${GTF_FILE/\.gtf/}

mkdir ${ASSEMBLY}
#mv ${GTF_FILE} ${ASSEMBLY}/${ASSEMBLY}.gff

perl ${DIR}/gtf2gff3.pl ${GTF_FILE} > ${ASSEMBLY}/${ASSEMBLY}.gff

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

