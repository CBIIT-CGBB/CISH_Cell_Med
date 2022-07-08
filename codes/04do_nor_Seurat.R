#!/usr/local/apps/R/3.5/3.5.0_build2/bin/Rscript --slave

rm(list=ls());

library(Seurat);
options(stringsAsFactors = F);
dat   <- read.table(gzfile("../raw_data/matrix_raw.txt.gz"), header=T, sep="\t");

dat1   <- read.table("mit_name.txt", header=T);
dat2   <- read.table("dup_name.txt", header=T);
dat3   <- read.table("ngene_name.txt", header=T);
dat.n  <- unique(c(dat1[,1], dat2[,1], dat3[,1]));
dat.i  <- which(colnames(dat) %in% dat.n);
dat.s  <- dat[,-dat.i];

dat.v  <- apply(dat.s, 1, var);
dat.i  <- which(dat.v == 0);
if(length(dat.i) > 0){
  dat.s  <- dat.s[-dat.i,];
}

dat.o <- CreateSeuratObject(raw.data = dat.s, min.cells = 0, min.genes = 0, 
                            project = "NCI");
dat.n <- NormalizeData(object = dat.o);
dat.m <- dat.n@data;
dat.m <- as.matrix(dat.m);

write.table(dat.m, gzfile("dat_nor_Seurat.txt.gz"), sep="\t", quote=F);
save(dat.m, "dat_nor_Seurat.RData");

cols <- rainbow(10, alpha=0.8);
pdf("matrix_nor_Seurat.pdf", 12, 6);
set.seed(1234);
dat.j <- sample(1:ncol(dat.m), 50);
boxplot(dat.m[,dat.j], col=cols[7]);
dev.off();


