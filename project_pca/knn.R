library(class)
library(tidyverse)

args    = commandArgs(trailingOnly=T)
infile  = args[1]
outfile = args[2]
outpdf  = args[2]

d       = read.table(infile, fill=T, header=T)
sample  = read.table("tmp/sample.panel", header=T)
tmp     = left_join(d, sample, by=c("ID" = "sample"))
tmp$pop = ifelse(is.na(tmp$super_pop), "tgt", "kgp")

train   = tmp %>% filter(pop=="kgp")
test    = tmp %>% filter(pop=="tgt")

set.seed(1111)
test$knn  = knn(train[,2:11], test[,2:11], train$super_pop)
train$knn = train$super_pop

out = rbind(test, train)

pdf(outpdf, height=5, width=10)

  out %>% ggplot(aes(x=PC1, y=PC2, color=factor(knn))) + geom_point() + facet_grid(pop~knn)
  out %>% ggplot(aes(x=PC2, y=PC3, color=factor(knn))) + geom_point() + facet_grid(pop~knn)
  out %>% ggplot(aes(x=PC3, y=PC4, color=factor(knn))) + geom_point() + facet_grid(pop~knn)

dev.off()

