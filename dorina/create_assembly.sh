#!/usr/bin/env bash
usage="USAGE:

  Bash script to prepare an annotation for doRiNA

  1) Download genePredToGtf from UCSC
    http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/genePredToGtf
    http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf

  2) Download gtf2gff3.pl
    http://genes.mit.edu/burgelab/miso/scripts/gtf2gff3.pl

  3) Have bedtools and mysql installed and accessible

  4) After making sure the exectables are accessible at the current working
     directory, use this script with a single parameter, a valid UCSC assembly,
     such as hg38, hg19 or mm10. It will generate a list of separated files,
     one for each feature.

  $(basename $0) [assembly] [-h]
    -h       : print this help and exit
    assembly : a valid UCSC assembly [default=hg38]"
# default values
ASSEMBLY=${1:-hg38}

mkdir ${ASSEMBLY}
pushd ${ASSEMBLY} > /dev/null
if [ "$(uname)" == "Darwin" ]; then
    wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/genePredToGtf --no-clobber
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf --no-clobber

wget http://genes.mit.edu/burgelab/miso/scripts/gtf2gff3.pl --no-clobber

rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/${ASSEMBLY}/database/refGene.txt.gz ${ASSEMBLY}.txt.gz
gzip -d ${ASSEMBLY}.txt.gz
cut -f 2- ${ASSEMBLY}.txt > ${ASSEMBLY}.input
./genePredToGtf file ${ASSEMBLY}.input ${ASSEMBLY}.gtf
perl gtf2gff3.pl ${ASSEMBLY}.gtf > ${ASSEMBLY}.gff
bedtools sort -i ${ASSEMBLY}.gff > tmp && mv tmp ${ASSEMBLY}.gff

grep gene ${ASSEMBLY}.gff > all.gff
grep CDS ${ASSEMBLY}.gff > cds.gff
grep three_prime ${ASSEMBLY}.gff > 3_utr.gff
grep five_prime ${ASSEMBLY}.gff > 5_utr.gff
grep exon ${ASSEMBLY}.gff > exon.gff

bedtools subtract -s -a all.gff -b exon.gff |sed -e "s/\tgene\t/\tintron\t/" > intron.gff
wget http://hgdownload.soe.ucsc.edu/goldenPath/${ASSEMBLY}/bigZips/${ASSEMBLY}.chrom.sizes --no-clobber
sort -k1,1 -k2,2n  ${ASSEMBLY}.chrom.sizes > tmp &&  mv tmp ${ASSEMBLY}.chrom.sizes
bedtools complement -i all.gff -g ${ASSEMBLY}.chrom.sizes > intergenic.bed

rm ${ASSEMBLY}.txt ${ASSEMBLY}.input ${ASSEMBLY}.gtf ${ASSEMBLY}.all
popd > /dev/null

