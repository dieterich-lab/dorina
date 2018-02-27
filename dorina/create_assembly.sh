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

mysql --host=genome-mysql.soe.ucsc.edu --user=genome \
-Ne "select a.name, a.chrom, a.strand, a.txStart, a.txEnd,\
a.cdsStart, a.cdsEnd, a.exonCount, a.exonStarts, a.exonEnds,\
 0 as score, b.geneSymbol from knownGene a join \
kgXref b on a.name=b.kgID" ${ASSEMBLY} > ${ASSEMBLY}.genePred

genePredToGtf file ${ASSEMBLY}.genePred ${ASSEMBLY}.knownGene.gtf -utr

perl gtf2gff3 ${ASSEMBLY}.knownGene.gtf > ${ASSEMBLY}/${ASSEMBLY}.gff

grep gene ${ASSEMBLY}.gff > all.gff
grep CDS ${ASSEMBLY}.gff > cds.gff
grep three_prime ${ASSEMBLY}.gff > 3_utr.gff
grep five_prime ${ASSEMBLY}.gff > 5_utr.gff
grep exon ${ASSEMBLY}.gff > exon.gff

bedtools subtract -s -a all.gff -b exon.gff |sed -e "s/\tgene\t/\tintron\t/" > intron.gff

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from ${ASSEMBLY}.chromInfo" > ${ASSEMBLY}.genome
# http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

bedtools complement -i all.gff -g ${ASSEMBLY}.genome > intergenic.bed

rm ${ASSEMBLY}.genePred ${ASSEMBLY}.knownGene.gtf
popd > /dev/null

