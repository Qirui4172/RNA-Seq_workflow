#!/usr/bin/env Rscript

#================================================================================================================================
Usage<-function(){
  cat("Usage: Rscript DiffExprComparison.r [reads_count] [sample_info] [adj_pvalue] [fold_change] [gene_list] [rm_sample]\n\n",

  "Parameters:\n",
  "[reads_count]      reads count matrix, \"Comb_PR_MI_readscount.txt\"\n",
  "[sample_info]      sample information file, \"sample.info\"\n",
  "[adj_pvalue]       adjusted Pvalue threshold for finding DEGs, i.e. 0.05, 0.1\n",
  "[fold_change]      fold change threshold for finding DEGs, i.e. 2, 1.5\n",
  "[gene_list]        gene list of interest, \"genelist.txt\"\n",
  "[rm_sample]        sample to be removed, i.e. \"none\"\n\n",

  "Example:\n",
  "Rscript DiffExprComparison.r Comb_PR_MI_readscount.txt sample.info 0.05 2 genelist.txt \"none\" \n\n",

  "The following files will be generated:\n",
  "Comb_PR_MI_deseq2norm.tsv: DESeq2-normalized RNA-seq reads count matrix, used for statistic comparison\n",
  "DEGs_CombinationvsDMSO_up.tsv, DEGs_CombinationvsDMSO_down.tsv, DEGs_MIvsDMSO_up.tsv, DEGs_MIvsDMSO_down.tsv: up-/down-regulated DEG lists of Combination/MI vs DMSO\n",
  "GeneNames_in_heatmap.tsv: gene names in the heatmap\n",
  "diffExpress_plots.pdf: plots file\n\n",

  "Function: Use DESeq2 to find differentially expressed genes (DEGs) and generate plots.\n",
  "Qirui Zhang (qirui.zhang@uni-greifswald.de)\n",
  "11-01-2022\n\n"
  )
}

args<-commandArgs(TRUE)
if(length(args)!=6){Usage();quit();}


cat("\n#===========================================================================================\n")
# Load libraries and arguments
time<-format(Sys.time(), format='%H:%M:%S %Y-%m-%d')
cat(time, "Loading libraries and arguments...\n\n")

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
pdf("diffExpress_plots.pdf", useDingbats=FALSE)
options(scipen=20)

readscount.file<-args[1]  # Comb_PR_MI_readscount.txt
sampleinfo.file<-args[2]  # sample.info
adjust.pvalue<-as.numeric(args[3])  # 0.05
fold.change<-as.numeric(args[4])  # 2
genelist.file<-args[5]  # "genelist.txt"
rm.sample<-as.character(args[6])  # "none"

# Read files
cat("Reading files...\n\n")
reads.count<-read.table(readscount.file, header=T, stringsAsFactors=F)
rownames(reads.count)<-reads.count$Gene
reads.count<-subset(reads.count, select=-c(Gene))
if(rm.sample!="none"){reads.count<-reads.count[,!(colnames(reads.count) %in% rm.sample)]}

sample.info<-read.table(sampleinfo.file, header=T, stringsAsFactors=F)
colnames(sample.info)<-c("Sample", "Group")
sample.info$Group<-as.factor(sample.info$Group)
if(rm.sample!="none"){sample.info<-sample.info[which(sample.info$Sample!=rm.sample),]}

genelist<-scan(genelist.file, what=(f1=""))


cat("\n#===========================================================================================\n")
# Normalize reads count

cat("\nNormalizing reads count...\n\n")
dds.full<-DESeqDataSetFromMatrix(countData=reads.count, colData=sample.info, design= ~ Group)
dds<-dds.full[rowSums(counts(dds.full))>=10,]
  # here should be corresponding to downstream step "dds.sub=dds.sub[rowSums(counts(dds.sub))>=10,]"; if removed too many genes here (i.e. rowSums(counts(dds.sub))>=50) and removed fewer genes in the later step (rowSums(counts(dds.sub))>=10), then some genes that are identified DEGs may be not existing in original dds and vst dataset, it gives warning for plotting in the downstream steps.

dds<-DESeq(dds)
normalized.counts<-as.data.frame(counts(dds, normalized=TRUE))
vst<-varianceStabilizingTransformation(dds, blind=FALSE)
write.table(as.data.frame(normalized.counts), "Comb_PR_MI_deseq2norm.tsv", row.names=T, col.names=T, quote=F, sep="\t")


cat("\n#===========================================================================================\n")
# Differential comparison
DiffComparison<-function(dds.full, treatment){
  cat("---------------------------------------------------------------\n")
  cat("Comparing", treatment, "vs DMSO...\n\n")
  dds.sub=dds.full[,dds.full$Group==treatment | dds.full$Group=="DMSO"]
  dds.sub=dds.sub[rowSums(counts(dds.sub))>=10,]
    # here should be corresponding to upstream step "dds<-dds.full[rowSums(counts(dds.full))>=10,]".
  dds.sub$Group=droplevels(dds.sub$Group)
  dds.sub$Group=relevel(dds.sub$Group, ref="DMSO")
  dds.sub=DESeq(dds.sub)
  res.sub=results(dds.sub, contrast=c("Group", treatment, "DMSO"), alpha=adjust.pvalue)
  summary(res.sub)
  title.MAplot=paste(treatment, " vs DMSO", sep="")
  plotMA(res.sub, main={title.MAplot}, ylim=c(-10,10))

  deg.sub.up=res.sub[!is.na(res.sub$padj) & res.sub$padj<adjust.pvalue & res.sub$log2FoldChange>=log2(fold.change),]
  deg.sub.down=res.sub[!is.na(res.sub$padj) & res.sub$padj<adjust.pvalue & res.sub$log2FoldChange <= -log2(fold.change),]
  deg.sub.up=as.data.frame(deg.sub.up)
  deg.sub.down=as.data.frame(deg.sub.down)
  res.sub=as.data.frame(res.sub)

  up.num=nrow(deg.sub.up)
  down.num=nrow(deg.sub.down)
  cat("For genes with adjusted Pvalue<", adjust.pvalue, "; foldchange>=", fold.change, "\n", sep="")
  cat("Up-regulated:", up.num, "\n")
  cat("Down-regulated:", down.num, "\n\n")

  title.deg.up=paste("DEGs_", treatment, "vsDMSO_up.tsv", sep="")
  title.deg.down=paste("DEGs_", treatment, "vsDMSO_down.tsv", sep="")
  title.res.sub=paste("AllGenes_", treatment, "vsDMSO.tsv", sep="")
  write.table(deg.sub.up, {title.deg.up}, col.names=T, row.names=T, quote=F, sep="\t")
  write.table(deg.sub.down, {title.deg.down}, col.names=T, row.names=T, quote=F, sep="\t")
  write.table(res.sub, {title.res.sub}, col.names=T, row.names=T, quote=F, sep="\t")

  out.list=list(res.sub, deg.sub.up, deg.sub.down)
  return(out.list)
}

diffresult.comb<-DiffComparison(dds.full, "Combination")  # Combination vs DMSO
res.comb<-diffresult.comb[[1]]
deg.comb.up<-diffresult.comb[[2]]
deg.comb.down<-diffresult.comb[[3]]

diffresult.mi<-DiffComparison(dds.full, "MI")  # MI vs DMSO
res.mi<-diffresult.mi[[1]]
deg.mi.up<-diffresult.mi[[2]]
deg.mi.down<-diffresult.mi[[3]]


cat("\n#===========================================================================================\n")
cat("Sample correlation heatmap...\n")
sampleDist<-dist(t(assay(vst)))
sampleDist.mx<-as.matrix(sampleDist)
rownames(sampleDist.mx)<-NULL
colnames(sampleDist.mx)<-sample.info$Sample
colors<-rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(100))
pheatmap(as.data.frame(sampleDist.mx), clustering_distance_rows=sampleDist, clustering_distance_cols=sampleDist, border_color=NA, col=colors, show_rownames=T)

cat("PCA plot...\n")
data<-plotPCA(vst, intgroup="Group", returnData=TRUE)
percentVar<-round(100*attr(data, "percentVar"))
pca.plot1<-ggplot(data, aes(PC1, PC2, color=Group))+geom_point(size=3)+theme_bw()+theme(panel.grid=element_blank(), panel.border=element_rect(size=0.5), text=element_text(size=18, family="sans", color="black"), axis.text=element_text(size=15, color="black"), legend.text=element_text(size=15))+xlab(paste0("PC1: ",percentVar[1],"% variance"))+ylab(paste0("PC2: ",percentVar[2],"% variance"))+coord_fixed(ratio=2)
pca.plot1

# use in-house script
pca.results<-prcomp(t(as.data.frame(assay(vst))), scale=TRUE)
percentVar<-round(((pca.results$sdev^2)/sum(pca.results$sdev^2))*100, 1)
data<-as.data.frame(pca.results$x)
data$Group<-sample.info$Group
pca.plot2<-ggplot(data, aes(PC1, PC2, colour=Group))+geom_point(size=4, shape=19)+geom_text_repel(aes(label=rownames(data)), size=3, color="black")+scale_color_manual(values=c("Combination"="#F8766D", "MI"="#00BA38", "DMSO"="#619CFF"))

pca.plot2+xlab(paste("PC1: ", percentVar[1], "% variance", sep=""))+ylab(paste("PC2: ", percentVar[2], "% variance", sep=""))+theme_bw()+theme(panel.grid=element_blank(), panel.border=element_rect(colour="black", size=0.5), text=element_text(size=18, color="black"), axis.title=element_text(size=18), axis.text=element_text(size=15, color="black"), legend.title=element_text(size=15), legend.text=element_text(size=15))+coord_fixed(ratio=1.5)


cat("\n#===========================================================================================\n")
cat("Volcano plot...\n")
VolcanoPlot<-function(res, volcano.title){
  volcano=as.data.frame(res)
#  volcano$significant=as.factor(ifelse(!is.na(volcano$padj) & volcano$padj < adjust.pvalue & abs(volcano$log2FoldChange)>=log2(fold.change), ifelse(volcano$log2FoldChange>=log2(fold.change), "Up", "Down"), "No"))
  volcano$padj=ifelse(is.na(volcano$padj),1,volcano$padj)

  # seperate groups
  volcano$group=rep(0)

##### NOT USE ###########################################################
if(FALSE){
  # group1: down-regulated within normal fc and padj range (blue dots)
  volcano[which(-log10(volcano$padj) > -log10(adjust.pvalue) & -log10(volcano$padj) <= y.max & volcano$log2FoldChange >= -fc_max.abs & volcano$log2FoldChange <= -log2(fold.change)), 8]=rep(1)

  # group2: up-regulated within normal fc and padj range (red dots)
  volcano[which(-log10(volcano$padj) > -log10(adjust.pvalue) & -log10(volcano$padj) <= y.max & volcano$log2FoldChange >= log2(fold.change) & volcano$log2FoldChange <= fc_max.abs), 8]=rep(2)

  # group3: non-significant within normal fc or padj range (grey dots)
  volcano[which((-log10(volcano$padj) <= -log10(adjust.pvalue) & abs(volcano$log2FoldChange) <= fc_max.abs) | (-log10(volcano$padj) > -log10(adjust.pvalue) & -log10(volcano$padj) < y.max & abs(volcano$log2FoldChange) <= log2(fold.change))), 8]=rep(3)

  # group4: non-significant outside of normal fc or padj range (grey trangle)
  volcano[which(-log10(volcano$padj) <= -log10(adjust.pvalue) & volcano$log2FoldChange < -fc_max.abs), c("log2FoldChange", "group")]=list(log2FoldChange = -fc_max.abs, group=4)
  volcano[which(-log10(volcano$padj) >= -log10(adjust.pvalue) & volcano$log2FoldChange > fc_max.abs), c("log2FoldChange", "group")]=list(log2FoldChange = fc_max.abs, group=4)
  volcano[which(-log10(volcano$padj) > y.max & abs(volcano$log2FoldChange) < log2(fold.change)), c("padj", "group")]=list(padj = 10^-y.max, group=4)

  # group5: down-regulated outside of normal fc or padj range (blue trangle)
  volcano[which(-log10(volcano$padj) > -log10(adjust.pvalue) & -log10(volcano$padj) <= y.max & volcano$log2FoldChange < -fc_max.abs), c("log2FoldChange", "group")]=list(log2FoldChange = -fc_max.abs, group=5)
  volcano[which(-log10(volcano$padj) > y.max & volcano$log2FoldChange < -fc_max.abs), c("log2FoldChange", "padj", "group")]=list(log2FoldChange = -fc_max.abs, padj = 10^-y.max, group=5)
  volcano[which(-log10(volcano$padj) > y.max & volcano$log2FoldChange <= -log2(fold.change) & volcano$log2FoldChange >= -fc_max.abs), c("padj", "group")]=list(padj = 10^-y.max, group=5)

  # group6: up-regulated outside of normal fc or padj range (red trangle)
  volcano[which(-log10(volcano$padj) > -log10(adjust.pvalue) & -log10(volcano$padj) <= y.max & volcano$log2FoldChange > fc_max.abs), c("log2FoldChange", "group")]=list(log2FoldChange = fc_max.abs, group=6)
  volcano[which(-log10(volcano$padj) > y.max & volcano$log2FoldChange > fc_max.abs), c("log2FoldChange", "padj", "group")]=list(log2FoldChange = fc_max.abs, padj = 10^-y.max, group=6)
  volcano[which(-log10(volcano$padj) > y.max & volcano$log2FoldChange >= log2(fold.change) & volcano$log2FoldChange <= fc_max.abs), c("padj", "group")]=list(padj = 10^-y.max, group=6)

  volcano.label<-volcano[which(volcano$group!=3 & volcano$group!=4),]
}
####################################################################################
  # seperate groups
  volcano[which(-log10(volcano$padj) > -log10(adjust.pvalue) & volcano$log2FoldChange <= -log2(fold.change)), 8]=rep(1)  # down-regulated
  volcano[which(-log10(volcano$padj) > -log10(adjust.pvalue) & volcano$log2FoldChange >= log2(fold.change)), 8]=rep(2)  # up-regulated
  volcano[which(-log10(volcano$padj) <= -log10(adjust.pvalue) | (-log10(volcano$padj) > -log10(adjust.pvalue) & abs(volcano$log2FoldChange) < log2(fold.change))), 8]=rep(3)  # non-significant
  volcano.label<-volcano[which(volcano$group!=3),]

  p<-ggplot(volcano, aes(log2FoldChange, -log10(padj)))+geom_point(data=volcano[which(volcano$group==1),], color="royalblue3", alpha=0.7)+geom_point(data=volcano[which(volcano$group==2),], color=brewer.pal(11, "RdYlBu")[2], alpha=0.7)+geom_point(data=volcano[which(volcano$group==3),], color="grey55", alpha=0.7)+geom_point(data=volcano[which(volcano$group==4),], shape=2, color="gray55", alpha=0.7)+geom_point(data=volcano[which(volcano$group==5),], shape=2, color="royalblue3", alpha=0.7)+geom_point(data=volcano[which(volcano$group==6),], shape=2, color=brewer.pal(11,"RdYlBu")[2], alpha=0.7)+geom_text_repel(data=volcano.label, aes(log2FoldChange, -log10(padj), label=rownames(volcano.label)), size=2, color="black", max.overlaps=getOption("ggrepel.max.overlaps", default=10))

  p+labs(title={volcano.title}, x="log2FoldChange", y="-log10(padj)")+geom_hline(yintercept=-log10(adjust.pvalue), linetype=2, color="black", size=0.3)+geom_vline(xintercept=c(-log2(fold.change), log2(fold.change)), linetype=2, color="black", size=0.3)+theme_bw()+theme(text=element_text(size=20, color="black"), panel.grid=element_blank(), panel.border=element_blank(), axis.line=element_line(colour="black", size=0.5), axis.ticks=element_line(colour="black"), axis.text=element_text(colour="black"), legend.position = c(.95, .95))
}

VolcanoPlot(res.comb, "Combination vs DMSO")
VolcanoPlot(res.mi, "MI vs DMSO")


cat("\n#===========================================================================================\n")
cat("Heatmap of DEGs...\n")
vst.df<-as.data.frame(assay(vst))
anno.label<-sample.info$Group
names(anno.label)<-sample.info$Sample
anno.label<-as.data.frame(anno.label)
names(anno.label)<-"Group"
anno.color<-list(Group=c("Combination"="#F8766D", "MI"="#00BFC4","DMSO"="#619CFF"))
colors<-colorRampPalette(c("blue4","white","red3"))(100)
#colors<-colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(3:9)]), brewer.pal(9, "OrRd")[c(3:9)]))(100)

### NOT use anymore ##################################################
if(FALSE){
subset_vst<-function(deg_sub, vst_full, prefix){
  vst_sub<-vst_full[rownames(deg_sub),]
  if(nrow(deg_sub)!=0){rownames(vst_sub)<-paste(prefix, "_", rownames(vst_sub), sep="")}
  return(vst_sub)
}
vst.comb.up<-subset_vst(deg.comb.up, vst.df, "comb_up")
vst.comb.down<-subset_vst(deg.comb.down, vst.df, "comb_down")
vst.mi.up<-subset_vst(deg.mi.up, vst.df, "mi_up")
vst.mi.down<-subset_vst(deg.mi.down, vst.df, "mi_down")
vst.4heatmap<-rbind(vst.comb.up, vst.comb.down, vst.mi.up, vst.mi.down)
row_gap1<-nrow(vst.comb.up)
row_gap2<-row_gap1+nrow(vst.comb.down)
row_gap3<-row_gap2+nrow(vst.mi.up)
row_gap4<-row_gap3+nrow(vst.mi.down)
pheatmap(vst.4heatmap, main="Heatmap of DEGs", scale="row", color=colors, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, fontsize=12, annotation_col=anno.label, annotation_colors=anno.color, annotation_names_col=F, border_color=NA, gaps_row=c(row_gap1, row_gap2, row_gap3, row_gap4))
}
#####################################################################

deg.names<-unique(c(rownames(deg.comb.up), rownames(deg.comb.down), rownames(deg.mi.up), rownames(deg.mi.down)))
vst.4heatmap2<-vst.df[deg.names,]
heatmap.p<-pheatmap(vst.4heatmap2, main="Heatmap of DEGs", scale="row", color=colors, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, fontsize=12, fontsize_row=0.5, annotation_col=anno.label, annotation_colors=anno.color, annotation_names_col=F, border_color=NA)
heatmap.p
gene_names.sorted<-rownames(vst.4heatmap2[heatmap.p$tree_row[["order"]],])

# dendrogram plot
plot(heatmap.p$tree_row, hang=-1, cex=0.2)
abline(h=6.2, col="red", lty=2, lwd=2)
gene_names.clust<-cutree(heatmap.p$tree_row, h=6.2)
gene_names.df<-data.frame(gene=gene_names.sorted, cluster=NA, stringsAsFactors=F)
gene_names.df[,"cluster"]<-paste("cluster", gene_names.clust[gene_names.df$gene], sep="")
write.table(gene_names.df, "GeneNames_in_heatmap.tsv", row.names=F, col.names=F, quote=F, sep="\t")


cat("\n#===========================================================================================\n")
cat("Dot plot of only interested genes...\n") 
intrsted_gene.comb<-res.comb[genelist,]
intrsted_gene.mi<-res.mi[genelist,]
intrsted_gene.df<-data.frame(log2FoldChange_comb=intrsted_gene.comb$log2FoldChange, log2FoldChange_mi=intrsted_gene.mi$log2FoldChange, padj_comb=intrsted_gene.comb$padj, padj_mi=intrsted_gene.mi$padj, row.names=rownames(intrsted_gene.comb))

intrsted_gene.df$Significance<-ifelse(abs(intrsted_gene.df$log2FoldChange_comb) >= log2(fold.change) & intrsted_gene.df$padj_comb < adjust.pvalue, "Combination", "No")
intrsted_gene.df$Significance<-ifelse(abs(intrsted_gene.df$log2FoldChange_mi) >= log2(fold.change) & intrsted_gene.df$padj_mi < adjust.pvalue, ifelse(intrsted_gene.df$Significance=="Combination", "Both", "MI-503"), intrsted_gene.df$Significance)

# plot
intrsted_gene.p<-ggplot(intrsted_gene.df, aes(log2FoldChange_comb, log2FoldChange_mi, color=Significance))+geom_point()+geom_text_repel(aes(label=rownames(intrsted_gene.df)), color="black", size=3, max.overlaps=getOption("ggrepel.max.overlaps", default=100))+scale_color_manual(values=c("Both"="#F8766D", "Combination"="#00BFC4", "MI-503"="#619CFF", "No"="#B79F00"))

intrsted_gene.p<-intrsted_gene.p+theme_bw()+theme(panel.grid=element_blank(), panel.border=element_rect(color="black", size=0.5), text=element_text(color="black", size=12), axis.title=element_text(size=12), axis.text=element_text(color="black", size=12), legend.title=element_text(size=12), legend.text=element_text(size=12))

intrsted_gene.p<-intrsted_gene.p+xlab("Combination (log2foldchange)")+ylab("MI-503 (log2foldchange)")+scale_x_continuous(breaks=c(-3, -2, -1, 0, 1, 2), labels=c("-3", "-2", "-1", "0", "1", "2"), limits=c(-3.5, 2.5))+scale_y_continuous(breaks=c(-3, -2, -1, 0, 1, 2), labels=c("-3", "-2", "-1", "0", "1", "2"), limits=c(-3.5, 2.5))+geom_hline(yintercept=-log2(fold.change), linetype=2, color="black", size=0.3)+geom_hline(yintercept=log2(fold.change), linetype=2, color="black", size=0.3)+geom_vline(xintercept=-log2(fold.change), linetype=2, color="black", size=0.3)+geom_vline(xintercept=log2(fold.change), linetype=2, color="black", size=0.3)+coord_fixed(ratio=1)

intrsted_gene.p


#------------------------------------------------------------------
cat("Dot plot of all genes...\n")
# extract "log2FoldChange" & "padj" columns
res.comb.tmp<-res.comb[,c("log2FoldChange", "padj")]
colnames(res.comb.tmp)<-c("log2fc_comb", "padj_comb")
res.mi.tmp<-res.mi[,c("log2FoldChange", "padj")]
colnames(res.mi.tmp)<-c("log2fc_mi", "padj_mi")

# use all genes that appeared in either res.comb or res.mi
allgene.names<-union(rownames(res.comb), rownames(res.mi))
dot.data<-data.frame(log2fc_comb=res.comb.tmp[allgene.names, "log2fc_comb"], padj_comb=res.comb.tmp[allgene.names, "padj_comb"], log2fc_mi=res.mi.tmp[allgene.names, "log2fc_mi"], padj_mi=res.mi.tmp[allgene.names, "padj_mi"], row.names=allgene.names)

# replace "log2fc=NA" with "log2fc=0"; replace "padj=NA" with "padj=1"
dot.data$log2fc_comb<-ifelse(is.na(dot.data$log2fc_comb), 0, dot.data$log2fc_comb)
dot.data$padj_comb<-ifelse(is.na(dot.data$padj_comb), 1, dot.data$padj_comb)
dot.data$log2fc_mi<-ifelse(is.na(dot.data$log2fc_mi), 0, dot.data$log2fc_mi)
dot.data$padj_mi<-ifelse(is.na(dot.data$padj_mi), 1, dot.data$padj_mi)

### NOT use anymore ##################################################
if(FALSE){
# judge whether genes significant or not
dot.data$Significance<-ifelse(dot.data$padj_comb<adjust.pvalue & abs(dot.data$log2fc_comb)>=log2(fold.change), "Combination", "No")
dot.data$Significance<-ifelse(dot.data$padj_mi<adjust.pvalue & abs(dot.data$log2fc_mi)>=log2(fold.change), ifelse(dot.data$Significance=="Combination", "Both", "MI-503"), dot.data$Significance)

# plot
intrsted_gene.p2<-ggplot(dot.data, aes(log2fc_comb, log2fc_mi))+geom_point(data=dot.data[!rownames(dot.data) %in% genelist,], aes(color=Significance), size=1, alpha=0.8)+scale_color_manual(values=c("Both"="#FF3333", "Combination"="#00CC66", "MI-503"="#0066FF", "No"="#666666"))+geom_point(data=dot.data[rownames(dot.data) %in% genelist,], color="black", size=1)+geom_text_repel(data=dot.data[genelist,], aes(label=genelist), color="black", size=3, max.overlaps=getOption("ggrepel.max.overlaps", default=100))

intrsted_gene.p2<-intrsted_gene.p2+theme_bw()+theme(panel.grid=element_blank(), panel.border=element_rect(color="black", size=0.5), text=element_text(color="black", size=12), axis.title=element_text(size=12), axis.text=element_text(color="black", size=12), legend.title=element_text(size=12), legend.text=element_text(size=12))

intrsted_gene.p2<-intrsted_gene.p2+xlab("Combination (log2foldchange)")+ylab("MI-503 (log2foldchange)")+scale_x_continuous(breaks=c(-6, -4, -2, 0, 2, 4, 6, 8), labels=c("-6", "-4", "-2", "0", "2", "4", "6", "8"), limits=c(-8, 8))+scale_y_continuous(breaks=c(-6, -4, -2, 0, 2, 4, 6, 8), labels=c("-6", "-4", "-2", "0", "2", "4", "6", "8"), limits=c(-8, 8))+geom_hline(yintercept=-log2(fold.change), linetype=2, color="black", size=0.3)+geom_hline(yintercept=log2(fold.change), linetype=2, color="black", size=0.3)+geom_vline(xintercept=-log2(fold.change), linetype=2, color="black", size=0.3)+geom_vline(xintercept=log2(fold.change), linetype=2, color="black", size=0.3)+geom_abline(intercept=0, slope=1, linetype=2, color="black", size=0.3)+coord_fixed(ratio=1)
}
####################################################################

# seperate genes to different areas (clusters) according to their log2fc
dot.data$Cluster<-rep(NA)
dot.data[which(dot.data$log2fc_mi<0 & dot.data$log2fc_comb<0 & dot.data$log2fc_mi>dot.data$log2fc_comb),"Cluster"]<-"leftdown"
dot.data[which(dot.data$log2fc_mi>0 & dot.data$log2fc_comb>=0 & dot.data$log2fc_mi>dot.data$log2fc_comb),"Cluster"]<-"leftup"
dot.data[which(dot.data$log2fc_mi<0 & dot.data$log2fc_comb<0 & dot.data$log2fc_mi<dot.data$log2fc_comb),"Cluster"]<-"rightdown"
dot.data[which(dot.data$log2fc_mi>=0 & dot.data$log2fc_comb>0 & dot.data$log2fc_mi<dot.data$log2fc_comb),"Cluster"]<-"rightup"
dot.data[which((dot.data$log2fc_mi>0 & dot.data$log2fc_comb<0) | (dot.data$log2fc_mi<0 & dot.data$log2fc_comb>0)),"Cluster"]<-"diagnal"
dot.data[is.na(dot.data$Cluster),"Cluster"]<-"diagnal"

# plot
Dotplot_Allgenes<-function(data, label_or_no){
  allgenes.p=ggplot(data, aes(log2fc_comb, log2fc_mi))+geom_point(data=dot.data[!rownames(dot.data) %in% genelist,], aes(color=Cluster), size=1, alpha=0.6)+scale_color_manual("Clusters", values=c("leftdown"="blue", "rightdown"="blue4", "leftup"="red", "rightup"="red4", "diagnal"="gray55"), breaks=c("leftdown", "rightdown", "leftup", "rightup"), labels=c("Comb down more", "MI down more", "MI up more", "Comb up more"))+geom_point(data=dot.data[rownames(dot.data) %in% genelist,], color="black", size=1)

  # add label or not
  if(label_or_no=="yes"){allgenes.p=allgenes.p+geom_text_repel(data=dot.data[genelist,], aes(label=genelist), color="black", size=3, max.overlaps=getOption("ggrepel.max.overlaps", default=100))}

  allgenes.p=allgenes.p+theme_bw()+theme(panel.grid=element_blank(), panel.border=element_rect(color="black", size=0.5), text=element_text(color="black", size=12), axis.title=element_text(size=12), axis.text=element_text(color="black", size=12), legend.title=element_text(size=10), legend.text=element_text(size=10))

  allgenes.p=allgenes.p+xlab("Combination (log2foldchange)")+ylab("MI-503 (log2foldchange)")+scale_x_continuous(breaks=c(-6, -4, -2, 0, 2, 4, 6, 8), labels=c("-6", "-4", "-2", "0", "2", "4", "6", "8"), limits=c(-7, 8))+scale_y_continuous(breaks=c(-6, -4, -2, 0, 2, 4, 6, 8), labels=c("-6", "-4", "-2", "0", "2", "4", "6", "8"), limits=c(-7, 8))+geom_hline(yintercept=0, linetype=2, color="black", size=0.3)+geom_vline(xintercept=0, linetype=2, color="black", size=0.3)+geom_abline(intercept=0, slope=1, linetype=2, color="black", size=0.3)+coord_fixed(ratio=1)

  allgenes.p
}
Dotplot_Allgenes(dot.data, "yes")
Dotplot_Allgenes(dot.data, "no")


#=================================================================================================================
cat("Done with analysis!\n\n")

