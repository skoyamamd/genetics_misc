#!/bin/bash

dir=.
infile=testinput.txt.gz
out=testout

dir=$(dirname $(readlink -f $0))
infile=$1
out=$2

{
  echo "##fileformat=VCFv4.2"
  echo "#CHROM POS ID REF ALT QUAL FILTER INFO" | tr " " "\t"
} > $out.input.vcf

zcat $infile | \
  awk -v OFS="\t" '{
    split($1,arr,":");
    print "chr"arr[1],arr[2],$1,arr[3],arr[4],".",".","."
  }' | sort -V >> $out.input.vcf

gzip -f $out.input.vcf

java -Xmx4G -jar $dir/ref/picard.jar LiftoverVcf \
  I=$out.input.vcf.gz \
  O=$out.hg38.vcf.gz \
  REJECT=$out.rejected.vcf.gz \
  CHAIN=$dir/ref/hg19ToHg38.over.chain.gz \
  R=$dir/ref/hg38.fa.gz

zcat $out.hg38.vcf.gz | grep -v "#" | \
  awk '{print $1":"$2":"$4":"$5"\t"$3}' | \
  sed '1s/^/#hg38\thg19\n/' | gzip -c > $out.hg38.txt.gz

rm $out.input.vcf.gz
rm $out.hg38.vcf.{gz,gz.tbi}

