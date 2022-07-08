#!/usr/local/apps/R/3.5/3.5.0_build2/bin/Rscript --slave

rm(list=ls());

options(stringsAsFactors = F);

dat   <- read.table(gzfile("../raw_data/matrix_raw.txt.gz"), header=T, row.names=1);
dat   <- as.matrix(dat);
anno  <- read.table("features.tsv");
row.names(anno) <- anno[,1];

an.i  <- match(row.names(dat), row.names(anno));
anno  <- anno[an.i,];
sum(row.names(dat)==row.names(anno))==nrow(dat);

mit.i <- grep("^MT-", anno[,2]);

anno[mit.i,];
dat.s <- dat[mit.i,];

mit.r <- nrow(dat.s)/nrow(dat);

col.v <- apply(dat,   2, sum);
mit.v <- apply(dat.s, 2, sum);

dat.f  <- mit.v/col.v;
cols  <- rainbow(10, alpha=0.5);
#hist(dat.f, breaks=100, col = cols[7], border=cols[7], xlab="Reads Proportion", 
#     main="Proportions of \nMitochondria Reads", xlim=c(0, 0.05));

cutoff <- 0.1;
dat.i  <- which(dat.f > cutoff);
length(dat.i);
mit.n <- colnames(dat)[dat.i];

dat.j <- which(colnames(dat) %in% mit.n);
length(dat.j);
mit.o <- colnames(dat)[dat.j];
write.table(data.frame(name1=mit.o, name2=mit.o), "mit_name.txt", quote=F, sep="\t", row.names=F);

pdf("mit.pdf", 4, 4);
hist(dat.f, breaks=100, col = cols[7], border=cols[7], xlab="Reads Proportion", 
     main="Proportions of \nMitochondria Reads");
abline(v=cutoff, col=cols[1], lwd=2);
text.s <- paste0("There are \n", length(dat.j), " cells that\nwere signed into\nthe death cells.")
text(0.6, 6000, text.s)
dev.off();
