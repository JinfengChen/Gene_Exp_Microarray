#http://www.bioconductor.org/help/course-materials/2005/BioC2005/labs/lab01/estrogen/
#http://bioinf.wehi.edu.au/marray/ibc2004/lab3/lab3.html#EstrogenToptable
#read data and preprocess
library(limma)
library(affy)
targets <- readTargets("estrogen.txt", sep="")
abatch <- ReadAffy(filenames=targets$filename)
eset <- rma(abatch)

#design matrix
f <- paste(targets$estrogen,targets$time.h,sep="")
f <- factor(f)
f
design <- model.matrix(~0+f)
colnames(design) <- levels(f)

#fit
fit <- lmFit(eset, design)

#contrast matrix
cont.matrix <- makeContrasts(E10="present10-absent10",E48="present48-absent48",Time="absent48-absent10",levels=design)

#extract linear model fit for the contrasts
fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

#Assessing differential expression
numGenes <- nrow(eset@assayData$exprs)
completelist <- topTable(fit2, coef=2, adjust="fdr", number=numGenes)
write.table(completelist, file="diffgene.list",sep="\t",quote=FALSE,col.names=NA)



