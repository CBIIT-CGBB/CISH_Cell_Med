#!/usr/local/apps/R/3.5/3.5.0_build2/bin/Rscript --slave

rm(list=ls());

options(stringasFoctors=F);
dat    <- read.table(gzfile("../raw_data/matrix_raw.txt.gz"), header=T, row.names=1);
mit    <- read.table("mit_name.txt", header=T);
mit.i  <- which(colnames(dat) %in% mit[,1]);
dat.s  <- dat[,-mit.i];
##
sum.l  <- apply(dat.s, 2, function(x) sum(x, na.rm=T));
cols <- rainbow(10, alpha=0.7);
#hist(sum.l, breaks=100, col = cols[7], border=cols[7], xlab="Total Reads Per Cell", 
#     main="");

cutoff <- 20000;
col.i  <- which(sum.l > cutoff);
out.s  <- data.frame(name1=colnames(dat.s)[col.i], name2=colnames(dat.s)[col.i]);
write.table(out.s, "dup_name.txt", quote=F, sep="\t", row.names=F);
num    <- length(col.i);

pdf("check_double.pdf", 4, 4);
hist(sum.l, breaks=100, col = cols[7], border=cols[7], xlab="Total Reads Per Cell", 
     main="");
abline(v=cutoff, col=cols[1], lwd=2);
text.s <- paste0("There are \n", num, " cells that\nwere signed into\nthe duplicated cells.")
text(55000, 8000, text.s)
dev.off();
