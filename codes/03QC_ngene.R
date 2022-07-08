#!/usr/local/apps/R/3.5/3.5.0_build2/bin/Rscript --slave

rm(list=ls());

options(stringsAsFactors = F);

dat   <- read.table(gzfile("../raw_data/matrix_raw.txt.gz"), header=T, row.names=1);
n1    <- read.table("dup_name.txt", header=T);
n2    <- read.table("mit_name.txt", header=T);
n     <- unique(c(n1[,1], n2[,1]));
n.i   <- which(colnames(dat) %in% n);
dat.s <- dat[,-n.i];

ngene1 <- apply(dat.s, 2, function(x) sum(x>1));
ngene0 <- apply(dat.s, 2, function(x) sum(x>0));
cols   <- rainbow(10, alpha=0.6);

cut1   <- 100;
cut0   <- 500;
g.n1   <- sum(ngene1 < cut1);
g.n0   <- sum(ngene0 < cut0);

g.i    <- which(ngene0 < cut0);
out.s  <- data.frame(name1=names(g.i), name2=names(g.i));
write.table(out.s, "ngene_name.txt", quote=F, sep="\t", row.names=F);

pdf("ngene.pdf", 8, 4);
par(mfrow=c(1,2));
hist(ngene0, breaks=100, col = cols[7], border=cols[7], xlab="gene number per cell", 
     main="Number of genes (reads>0)");
abline(v=cut0, col=cols[1], lwd=2);
t.s0 <- paste0(g.n0, " of ", dim(dat.s)[1], "\n", "cutoff: ", cut0); 
text(3000, 1500, t.s0);
hist(ngene1, breaks=100, col = cols[7], border=cols[7], xlab="gene number per cell", 
     main="Number of genes (reads>1)");
abline(v=cut1, col=cols[1], lwd=2);
t.s1 <- paste0(g.n1, " of ", dim(dat.s)[1], "\n", "cutoff: ", cut1) 
text(1200, 1500, t.s1);
dev.off();

