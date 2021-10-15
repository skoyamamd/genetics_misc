#!/bin/bash
set -e
set pipefail

####################################
## Create qced 1KG genotypes      ##
## dependencies plink2,plink,wget ##
####################################

mkdir -p tmp out

dir=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502
wget -O tmp/sample.panel $dir/integrated_call_samples_v3.20130502.ALL.panel

hd=tmp/KGP.tmp

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
    --out $hd.$chr

  cat $hd.$chr.pvar | grep -v "#" | \
    cut -f 1,2 | sort | uniq -d | \
    awk -v OFS="\t" '{print $1,$2-1,$2}' > $hd.$chr.dup.bed

  cat $hd.$chr.pvar | grep -v "#" | \
    awk '($4~/A|T/ && $5~/T|A/) || ($4~/C|G/ && $5~/C|G/)' | \
    awk -v OFS="\t" '{print $1,$2-1,$2}' >> $hd.$chr.dup.bed

  plink2 \
    --pfile $hd.$chr \
    --exclude bed0 $hd.$chr.dup.bed \
    --make-pgen \
    --set-all-var-ids '@:#:$r:$a' \
    --threads 4 \
    --memory 4800 \
    --out $hd.nondup.$chr

  plink2 \
    --pfile $hd.nondup.$chr \
    --keep-allele-order \
    --make-bed \
    --threads 4 \
    --memory 4800 \
    --out $hd.nondup.$chr

  rm tmp/$file

done

##########################
## Step 2: merge bfiles ##
##########################

ls $hd.nondup.*.bed | sort -V | sed 's/.bed//' | head -n 1 > $hd.merge.head
ls $hd.nondup.*.bed | sort -V | sed 's/.bed//' | sed '1d'  > $hd.merge.list

plink \
  --bfile $(cat $hd.merge.head) \
  --merge-list $hd.merge.list \
  --threads 4 \
  --memory 4800 \
  --out $hd.merged

###############################################################
## Step 3: population wise Hardy-Weinberg equilibrium filter ##
###############################################################

cat tmp/sample.panel | sed '1d' | cut -f 3 | sort | uniq | while read pop; do

  plink2 \
    --bfile $hd.merged \
    --hwe 1e-6 \
    --hardy \
    --keep <(cat tmp/sample.panel | sed '1d' | grep -w $pop | cut -f 1) \
    --threads 4 \
    --memory 4800 \
    --out $hd.$pop > /dev/null

  cat $hd.$pop.hardy | awk '$10<1e-6' > $hd.$pop.hwd

  Nhwd=$(cat $hd.$pop.hwd | wc -l)

  echo $pop $Nhwd

done

#########################
## Step 4: Final pfile ##
#########################

cat $hd.*.hwd | cut -f 2 | sort | uniq > $hd.merged.hwd

plink2 \
  --bfile $hd.merged \
  --exclude $hd.merged.hwd \
  --make-pgen \
  --threads 4 \
  --memory 4800 \
  --out out/KGP.qced && rm $hd.*

