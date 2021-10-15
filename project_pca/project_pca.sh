##########################
## Useage               ##
## bash 10_prep.sh \    ##
##  <<kgp.pgen>> \      ##
##  <<target.vcf.gz>> \ ##
##  <<output_prefix>>   ##
##########################

kgp=$1
vcf=$2
out=$3

plink2 \
  --vcf $vcf \
  --make-pgen \
  --geno 0.01 \
  --maf 0.01 \
  --set-all-var-ids '@:#:$r:$a' \
  --rm-dup force-first \
  --snps-only just-acgt \
  --new-id-max-allele-len 10000 \
  --out $out.tmp

plink2 \
  --pfile $kgp \
  --extract <(cat $out.tmp.pvar |cut -f 3) \
  --indep-pairwise 50 10 0.2 \
  --out $out.tmp

plink2 \
  --pfile $kgp \
  --extract $out.tmp.prune.in \
  --make-pgen \
  --out $out.tmp.kgp_pruned

plink2 \
  --pfile $out.tmp.kgp_pruned \
  --freq counts \
  --pca allele-wts \
  --out $out.tmp.kgp_pruned

## Score target population

plink2 \
  --pfile $out.tmp \
  --read-freq $out.tmp.kgp_pruned.acount \
  --score $out.tmp.kgp_pruned.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
  --rm-dup force-first \
  --score-col-nums 6-15 \
  --out $out.tmp.tgt

## Score kgp

plink2 \
  --pfile $kgp \
  --read-freq $out.tmp.kgp_pruned.acount \
  --score $out.tmp.kgp_pruned.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
  --rm-dup force-first \
  --score-col-nums 6-15 \
  --out $out.tmp.kgp

## merge the data

echo ID PC{1..10} | tr " " "\t" | gzip -c > $out.kgp_projected.txt.gz

cat \
  <(cat $out.tmp.kgp.sscore | sed '1d' | cut -f 1,4-) \
  <(cat $out.tmp.tgt.sscore | sed '1d' | cut -f 1,4-) | \
  gzip -c >> $out.kgp_projected.txt.gz

rm $out.tmp.*

Rscript knn.R \
  $out.kgp_projected.txt.gz \
  $out.kgp_projected.pdf

