
mkdir -p ref && cd ref

wget -nc https://github.com/broadinstitute/picard/releases/download/2.26.4/picard.jar
wget -nc https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
wget -nc https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget -nc https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
wget -nc https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz

picardjar=/apps/lib/picard/2.7.1/picard.jar

java -jar picard.jar \
  CreateSequenceDictionary \
  R=hg38.fa.gz \
  O=hg38.fa.gz.dict

java -jar picard.jar \
  CreateSequenceDictionary \
  R=hg19.fa.gz \
  O=hg19.fa.gz.dict

