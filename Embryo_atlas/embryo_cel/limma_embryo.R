#read data and preprocess
library(limma)
library(affy)
targets <- readTargets("embryo.txt", sep="\t")
abatch <- ReadAffy(filenames=targets$filename)
eset <- rma(abatch)

#design matrix
f <- targets$embryo
f <- factor(f)
design <- model.matrix(~0+f)
colnames(design) <- levels(f)

#fit
fit <- lmFit(eset, design)

#contrast matrix
cont.matrix <- makeContrasts(E21="S2_3_4DAP-S1_0_2DAP", E31="S3_5_10DAP-S1_0_2DAP", E41="S4_11_20DAP-S1_0_2DAP", E51="S5_21_29DAP-S1_0_2DAP", levels=design)

#extract linear model fit for the contrasts
fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

#Assessing differential expression
numGenes <- nrow(eset@assayData$exprs)
completelist <- topTable(fit2, coef=2, adjust="fdr", number=numGenes)
write.table(completelist, file="S3_5_10DAP-S1_0_2DAP.limma.list",sep="\t",quote=FALSE,col.names=NA)
completelist <- topTable(fit2, coef=1, adjust="fdr", number=numGenes)
write.table(completelist, file="S2_3_4DAP-S1_0_2DAP.limma.list",sep="\t",quote=FALSE,col.names=NA)
completelist <- topTable(fit2, coef=3, adjust="fdr", number=numGenes)
write.table(completelist, file="S4_11_20DAP-S1_0_2DAP.limma.list",sep="\t",quote=FALSE,col.names=NA)
completelist <- topTable(fit2, coef=4, adjust="fdr", number=numGenes)
write.table(completelist, file="S5_21_29DAP-S1_0_2DAP.limma.list",sep="\t",quote=FALSE,col.names=NA)


