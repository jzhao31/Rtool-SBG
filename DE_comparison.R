library(limma)

args <- commandArgs(trailingOnly=TRUE)

lfc_cutoff <- as.numeric(args[5])

p_cutoff <- as.numeric(args[6])

adjp_cutoff <- as.numeric(args[7])

Cuff_all <- read.table(args[1],header=T)

#TNBC vs HER2
Cuff_TNBC_HER2 <- Cuff_all[which(Cuff_all$sample_1=="TNBC_Tumor" & Cuff_all$sample_2=="HER2_Tumor"),]

DESeq_TNBC_HER2 <- read.table(args[2],sep=",",header=T)

Merge_TNBC_HER2 <- merge(Cuff_TNBC_HER2,DESeq_TNBC_HER2,by.x="gene",by.y="X")

Merge_TNBC_HER2_NoEmpty <- Merge_TNBC_HER2[complete.cases(Merge_TNBC_HER2[,c('log2.fold_change.','p_value','q_value','log2FoldChange','pvalue','padj')]),]

Sig_TNBC_HER2 <- matrix(rep(0),nrow(Merge_TNBC_HER2_NoEmpty),2)
rownames(Sig_TNBC_HER2) <- Merge_TNBC_HER2_NoEmpty$gene

for(i in 1:nrow(Merge_TNBC_HER2_NoEmpty)){
	if((abs(Merge_TNBC_HER2_NoEmpty$log2.fold_change.[i])>lfc_cutoff & Merge_TNBC_HER2_NoEmpty$p_value[i]<p_cutoff & Merge_TNBC_HER2_NoEmpty$q_value[i]<adjp_cutoff) & (abs(Merge_TNBC_HER2_NoEmpty$log2FoldChange[i])>lfc_cutoff & Merge_TNBC_HER2_NoEmpty$pvalue[i]<p_cutoff & Merge_TNBC_HER2_NoEmpty$padj[i]<adjp_cutoff)){
		Sig_TNBC_HER2[i,1] <- 1
		Sig_TNBC_HER2[i,2] <- 1
	}
	else if((abs(Merge_TNBC_HER2_NoEmpty$log2.fold_change.[i])>lfc_cutoff & Merge_TNBC_HER2_NoEmpty$p_value[i]<p_cutoff & Merge_TNBC_HER2_NoEmpty$q_value[i]<adjp_cutoff) & !(abs(Merge_TNBC_HER2_NoEmpty$log2FoldChange[i])>lfc_cutoff & Merge_TNBC_HER2_NoEmpty$pvalue[i]<p_cutoff & Merge_TNBC_HER2_NoEmpty$padj[i]<adjp_cutoff)){
		Sig_TNBC_HER2[i,1] <- 1
		Sig_TNBC_HER2[i,2] <- 0
	}
	else if(!(abs(Merge_TNBC_HER2_NoEmpty$log2.fold_change.[i])>lfc_cutoff & Merge_TNBC_HER2_NoEmpty$p_value[i]<p_cutoff & Merge_TNBC_HER2_NoEmpty$q_value[i]<adjp_cutoff) & (abs(Merge_TNBC_HER2_NoEmpty$log2FoldChange[i])>lfc_cutoff & Merge_TNBC_HER2_NoEmpty$pvalue[i]<p_cutoff & Merge_TNBC_HER2_NoEmpty$padj[i]<adjp_cutoff)){
		Sig_TNBC_HER2[i,1] <- 0
		Sig_TNBC_HER2[i,2] <- 1
	}
		else {
		Sig_TNBC_HER2[i,1] <- 0
		Sig_TNBC_HER2[i,2] <- 0
	}
}

Venn_TNBC_HER2 <- vennCounts(Sig_TNBC_HER2)



#HER2 vs NonTNBC
Cuff_HER2_NonTNBC <- Cuff_all[which(Cuff_all$sample_1=="HER2_Tumor" & Cuff_all$sample_2=="NonTNBC_Tumor"),]

DESeq_HER2_NonTNBC <- read.table(args[3],sep=",",header=T)

Merge_HER2_NonTNBC <- merge(Cuff_HER2_NonTNBC,DESeq_HER2_NonTNBC,by.x="gene",by.y="X")

Merge_HER2_NonTNBC_NoEmpty <- Merge_HER2_NonTNBC[complete.cases(Merge_HER2_NonTNBC[,c('log2.fold_change.','p_value','q_value','log2FoldChange','pvalue','padj')]),]

Sig_HER2_NonTNBC <- matrix(rep(0),nrow(Merge_HER2_NonTNBC_NoEmpty),2)
rownames(Sig_HER2_NonTNBC) <- Merge_HER2_NonTNBC_NoEmpty$gene

for(i in 1:nrow(Merge_HER2_NonTNBC_NoEmpty)){
	if((abs(Merge_HER2_NonTNBC_NoEmpty$log2.fold_change.[i])>lfc_cutoff & Merge_HER2_NonTNBC_NoEmpty$p_value[i]<p_cutoff & Merge_HER2_NonTNBC_NoEmpty$q_value[i]<adjp_cutoff) & (abs(Merge_HER2_NonTNBC_NoEmpty$log2FoldChange[i])>lfc_cutoff & Merge_HER2_NonTNBC_NoEmpty$pvalue[i]<p_cutoff & Merge_HER2_NonTNBC_NoEmpty$padj[i]<adjp_cutoff)){
		Sig_HER2_NonTNBC[i,1] <- 1
		Sig_HER2_NonTNBC[i,2] <- 1
	}
	else if((abs(Merge_HER2_NonTNBC_NoEmpty$log2.fold_change.[i])>lfc_cutoff & Merge_HER2_NonTNBC_NoEmpty$p_value[i]<p_cutoff & Merge_HER2_NonTNBC_NoEmpty$q_value[i]<adjp_cutoff) & !(abs(Merge_HER2_NonTNBC_NoEmpty$log2FoldChange[i])>lfc_cutoff & Merge_HER2_NonTNBC_NoEmpty$pvalue[i]<p_cutoff & Merge_HER2_NonTNBC_NoEmpty$padj[i]<adjp_cutoff)){
		Sig_HER2_NonTNBC[i,1] <- 1
		Sig_HER2_NonTNBC[i,2] <- 0
	}
	else if(!(abs(Merge_HER2_NonTNBC_NoEmpty$log2.fold_change.[i])>lfc_cutoff & Merge_HER2_NonTNBC_NoEmpty$p_value[i]<p_cutoff & Merge_HER2_NonTNBC_NoEmpty$q_value[i]<adjp_cutoff) & (abs(Merge_HER2_NonTNBC_NoEmpty$log2FoldChange[i])>lfc_cutoff & Merge_HER2_NonTNBC_NoEmpty$pvalue[i]<p_cutoff & Merge_HER2_NonTNBC_NoEmpty$padj[i]<adjp_cutoff)){
		Sig_HER2_NonTNBC[i,1] <- 0
		Sig_HER2_NonTNBC[i,2] <- 1
	}
		else {
		Sig_HER2_NonTNBC[i,1] <- 0
		Sig_HER2_NonTNBC[i,2] <- 0
	}
}

Venn_HER2_NonTNBC <- vennCounts(Sig_HER2_NonTNBC)




#TNBC vs nonTNBC
Cuff_TNBC_NonTNBC <- Cuff_all[which(Cuff_all$sample_1=="TNBC_Tumor" & Cuff_all$sample_2=="NonTNBC_Tumor"),]

DESeq_TNBC_NonTNBC <- read.table(args[4],sep=",",header=T)

Merge_TNBC_NonTNBC <- merge(Cuff_TNBC_NonTNBC,DESeq_TNBC_NonTNBC,by.x="gene",by.y="X")

Merge_TNBC_NonTNBC_NoEmpty <- Merge_TNBC_NonTNBC[complete.cases(Merge_TNBC_NonTNBC[,c('log2.fold_change.','p_value','q_value','log2FoldChange','pvalue','padj')]),]

Sig_TNBC_NonTNBC <- matrix(rep(0),nrow(Merge_TNBC_NonTNBC_NoEmpty),2)
rownames(Sig_TNBC_NonTNBC) <- Merge_TNBC_NonTNBC_NoEmpty$gene

for(i in 1:nrow(Merge_TNBC_NonTNBC_NoEmpty)){
	if((abs(Merge_TNBC_NonTNBC_NoEmpty$log2.fold_change.[i])>lfc_cutoff & Merge_TNBC_NonTNBC_NoEmpty$p_value[i]<p_cutoff & Merge_TNBC_NonTNBC_NoEmpty$q_value[i]<adjp_cutoff) & (abs(Merge_TNBC_NonTNBC_NoEmpty$log2FoldChange[i])>lfc_cutoff & Merge_TNBC_NonTNBC_NoEmpty$pvalue[i]<p_cutoff & Merge_TNBC_NonTNBC_NoEmpty$padj[i]<adjp_cutoff)){
		Sig_TNBC_NonTNBC[i,1] <- 1
		Sig_TNBC_NonTNBC[i,2] <- 1
	}
	else if((abs(Merge_TNBC_NonTNBC_NoEmpty$log2.fold_change.[i])>lfc_cutoff & Merge_TNBC_NonTNBC_NoEmpty$p_value[i]<p_cutoff & Merge_TNBC_NonTNBC_NoEmpty$q_value[i]<adjp_cutoff) & !(abs(Merge_TNBC_NonTNBC_NoEmpty$log2FoldChange[i])>lfc_cutoff & Merge_TNBC_NonTNBC_NoEmpty$pvalue[i]<p_cutoff & Merge_TNBC_NonTNBC_NoEmpty$padj[i]<adjp_cutoff)){
		Sig_TNBC_NonTNBC[i,1] <- 1
		Sig_TNBC_NonTNBC[i,2] <- 0
	}
	else if(!(abs(Merge_TNBC_NonTNBC_NoEmpty$log2.fold_change.[i])>lfc_cutoff & Merge_TNBC_NonTNBC_NoEmpty$p_value[i]<p_cutoff & Merge_TNBC_NonTNBC_NoEmpty$q_value[i]<adjp_cutoff) & (abs(Merge_TNBC_NonTNBC_NoEmpty$log2FoldChange[i])>lfc_cutoff & Merge_TNBC_NonTNBC_NoEmpty$pvalue[i]<p_cutoff & Merge_TNBC_NonTNBC_NoEmpty$padj[i]<adjp_cutoff)){
		Sig_TNBC_NonTNBC[i,1] <- 0
		Sig_TNBC_NonTNBC[i,2] <- 1
	}
		else {
		Sig_TNBC_NonTNBC[i,1] <- 0
		Sig_TNBC_NonTNBC[i,2] <- 0
	}
}

Venn_TNBC_NonTNBC <- vennCounts(Sig_TNBC_NonTNBC)






# results

pdf()

vennDiagram(Venn_TNBC_HER2, names=c("STAR-TopHat-Cufflinks","STAR-HTSeq-DESeq2"),cex=0.8)

mtext(paste("TNBC vs HER2\n","\n","Significance level --log2 fold change:",args[5],", p-value:",args[6],", adjusted p-value:",args[7]),side=3,line=-3)

mtext(paste("The number of DEG by STAR-TopHat-Cufflinks is",Venn_TNBC_HER2[3,3]+Venn_TNBC_HER2[4,3],"\n","The number of DEG by STAR-HTSeq-DESeq2 is",Venn_TNBC_HER2[2,3]+Venn_TNBC_HER2[4,3],"\n","There are",Venn_TNBC_HER2[4,3],"overlapped differential expressed genes"),side=1)

vennDiagram(Venn_HER2_NonTNBC, names=c("STAR-TopHat-Cufflinks","STAR-HTSeq-DESeq2"),cex=0.8)

mtext(paste("HER2 vs NonTNBC\n","\n","Significance level --log2 fold change:",args[5],", p-value:",args[6],", adjusted p-value:",args[7]),side=3,line=-3)

mtext(paste("The number of DEG by STAR-TopHat-Cufflinks is",Venn_HER2_NonTNBC[3,3]+Venn_HER2_NonTNBC[4,3],"\n","The number of DEG by STAR-HTSeq-DESeq2 is",Venn_HER2_NonTNBC[2,3]+Venn_HER2_NonTNBC[4,3],"\n","There are",Venn_HER2_NonTNBC[4,3],"overlapped differential expressed genes"),side=1)

vennDiagram(Venn_TNBC_NonTNBC, names=c("STAR-TopHat-Cufflinks","STAR-HTSeq-DESeq2"),cex=0.8)

mtext(paste("TNBC vs NonTNBC\n","\n","Significance level --log2 fold change:",args[5],", p-value:",args[6],", adjusted p-value:",args[7]),side=3,line=-3)

mtext(paste("The number of DEG by STAR-TopHat-Cufflinks is",Venn_TNBC_NonTNBC[3,3]+Venn_TNBC_NonTNBC[4,3],"\n","The number of DEG by STAR-HTSeq-DESeq2 is",Venn_TNBC_NonTNBC[2,3]+Venn_TNBC_NonTNBC[4,3],"\n","There are",Venn_TNBC_NonTNBC[4,3],"overlapped differential expressed genes"),side=1)

dev.off()



Merge_cutoff_TNBC_HER2 <- Merge_TNBC_HER2_NoEmpty[which(Sig_TNBC_HER2[,1]==1 & Sig_TNBC_HER2[,2]==1),]

write.table(file="TNBC-HER2-overlapped_DEG.txt",Merge_cutoff_TNBC_HER2)

Merge_cutoff_HER2_NonTNBC <- Merge_HER2_NonTNBC_NoEmpty[which(Sig_HER2_NonTNBC[,1]==1 & Sig_HER2_NonTNBC[,2]==1),]

write.table(file="HER2-NonTNBC-overlapped_DEG.txt",Merge_cutoff_HER2_NonTNBC)

Merge_cutoff_TNBC_NonTNBC <- Merge_TNBC_NonTNBC_NoEmpty[which(Sig_TNBC_NonTNBC[,1]==1 & Sig_TNBC_NonTNBC[,2]==1),]

write.table(file="TNBC-NonTNBC-overlapped_DEG.txt",Merge_cutoff_TNBC_NonTNBC)




