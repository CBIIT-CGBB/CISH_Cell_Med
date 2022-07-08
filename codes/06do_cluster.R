rm(list=ls());

library(GCluster);
wt    <- 5;

files <- dir("tsne", ".txt", full.names = T);
for (inff in files){
  dat   <- read.table(inff, header=T, row.names=1);
  set.seed(1234);
  out   <- GCluster(dat=dat, wt=wt);
  dat.s <- data.frame(dat, clu=out$membership);
  outf  <- gsub("tsne", "cluster", inff)
  write.table(dat.s, outf, quote=F, sep="\t")
}

