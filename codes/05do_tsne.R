rm(list=ls());

options(stringsAsFactors = F);
library(Rtsne);

dat   <- read.table(gzfile("../data_set2/dat_nor_Seurat.txt.gz"), header=T, row.names=1);
n1     <- strsplit(colnames(dat), "\\.");
n2     <- do.call(rbind, n1);
n.u    <- unique(n2[,1]);

per_    <- 30;
max_    <- 2000;

for (n in n.u){
  dat.i <- grep(n, colnames(dat));
  dat.s <- dat[,dat.i];
  dat.t <- t(dat.s);
  outf  <- paste0("tsne/", n, ".txt");
  tsne       <- Rtsne(dat.t, dims = 2, perplexity=per_, verbose=TRUE, max_iter = max_);
  out        <- tsne$Y;
  row.names(out) <- colnames(dat.s);
  write.table(out, outf, sep="\t", quote=F);
}

