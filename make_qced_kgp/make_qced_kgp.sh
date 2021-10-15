#!/bin/bash

####################################
## Create qced 1KG genotypes      ##
## dependencies plink2,plink,wget ##
####################################

mkdir -p tmp out

dir=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502
wget -O tmp/sample.panel  $dir/integrated_call_samples_v3.20130502.ALL.panel

#######################
## Step 1: QC bfiles ##
#######################

for chr in {1..22}; do

  file=ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

  wget -O tmp/$file $dir/$file

  plink2 \
    --vcf tmp/$file \
    --snps-only just-acgt \
    --max-alleles 2 \
    --make-pgen \
    --maf 0.01 \
    --geno 0.01 \
    --threads 4 \
    --memory 4800 \
    --out tmp/KGP.$chr

  cat tmp/KGP.$chr.pvar | grep -v "#" | \
    cut -f 1,2 | sort | uniq -d | \
    awk -v OFS="\t" '{print $1,$2-1,$2}' > tmp/KGP.$chr.dup.bed

  cat ../tmp/KGP.$chr.pvar | grep -v "#" | \
    awk '($4~/A|T/ && $5~/T|A/) || ($4~/C|G/ && $5~/C|G/)' | \
    awk -v OFS="\t" '{print $1,$2-1,$2}' >> tmp/KGP.$chr.dup.bed

  plink2 \
    --pfile tmp/KGP.$chr \
    --exclude bed0 tmp/KGP.$chr.dup.bed \
    --make-pgen \
    --set-all-var-ids '@:#:$r:$a' \
    --threads 4 \
    --memory 4800 \
    --out tmp/KGP.nondup.$chr

  plink2 \
    --pfile tmp/KGP.nondup.$chr \
    --keep-allele-order \
    --make-bed \
    --threads 4 \
    --memory 4800 \
    --out tmp/KGP.nondup.$chr

  rm tmp/KGP.$chr.{pvar,psam,pgen}
  rm tmp/KGP.nondup.$chr.{pvar,psam,pgen}
  rm tmp/KGP.$chr.dup.bed
  rm tmp/$file

done

##########################
## Step 2: merge bfiles ##
##########################

ls tmp/KGP.nondup.*.bed | sort -V | sed 's/.bed//' | head -n 1 > tmp/merge.head
ls tmp/KGP.nondup.*.bed | sort -V | sed 's/.bed//' | sed '1d' > tmp/merge.list

plink \
  --bfile $(cat tmp/merge.head) \
  --merge-list tmp/merge.list \
  --threads 4 \
  --memory 4800 \
  --out tmp/KGP.merged

###############################################################
## Step 3: population wise Hardy-Weinberg equilibrium filter ##
###############################################################

cat tmp/sample.panel | sed '1d' | cut -f 3 | sort | uniq | while read pop; do

  plink2 \
    --bfile tmp/KGP.merged \
    --hwe 1e-6 \
    --hardy \
    --keep <(cat tmp/sample.panel | sed '1d' | grep -w $pop | cut -f 1) \
    --threads 4 \
    --memory 4800 \
    --out tmp/KGP.merged.$pop > /dev/null

  cat tmp/KGP.merged.$pop.hardy | awk '$10<1e-6' > tmp/$pop.hwd

  Nhwd=$(cat tmp/$pop.hwd | wc -l)

  echo $pop $Nhwd

done

#########################
## Step 4: Final pfile ##
#########################

cat tmp/*.hwd | cut -f 2 | \
  sort | uniq > tmp/merged.hwd

plink2 \
  --bfile tmp/KGP.qced \
  --exclude tmp/merged.hwd \
  --make-pgen \
  --threads 4 \
  --memory 4800 \
  --out out/KGP.qced

rm tmp/KGP.nondup.*
rm tmp/KGP.merged.*
rm tmp/merge.{head,list}
rm tmp/KGP.*.log
rm tmp/*.hwd

