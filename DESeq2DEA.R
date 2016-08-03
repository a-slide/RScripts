#! /usr/bin/Rscript

###################################################################################################
# VERIFY ARG AND LOAD LIBRAIRIES
###################################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length (args) != 1)
    stop("USAGE: RNA_Seq.R <csv file containing : SampleName FileName Condition>")

cat ("### LOAD LIBRARIES ###\n\n")

options(width = 160) # wider terminal size
library("DESeq2")
library("gplots")
library("ggtern")
library("RColorBrewer")
library("FactoMineR")

###################################################################################################
# FUNCTIONS
###################################################################################################

run_analysis <- function(pdata, group) {

    # Start analysis of the entire group
    cat (paste ("\nANALYSING GROUP: ", group, "\n"))

    # Drop the unused conditions due to the subset function and extract the condictions or DESeq2
    pdata$Condition <- factor(pdata$Condition)
    conds = pdata$Condition
    nconds = nlevels (pdata$Condition)
    nind = nrow (pdata)

    cat (paste ("\tNumber of conditions", nconds, "\n"))

    if (nconds < 2) {
        cat ("\t\tNot enough conditions \n")
        return()
        }

    # Create a directory to store files related to this group
    dir.create(file.path("./", group), showWarnings = FALSE)
    prefix = paste0("./", group, "/", group)

    # Make unique colors to correspond to each condition
    cond_colors = rainbow(nconds, alpha = 0.6)[conds]
    names(cond_colors) = conds

    # Make A DESEQ2 object direct from the HTSeq counts tables
    cds = DESeqDataSetFromHTSeqCount(sampleTable=pdata, directory = ".", design= ~ Condition)

    # Perform the DE analysis
    dds = DESeq(cds, fitType="mean")
    res = results(dds)
    # Remove NA and sort by adjusted p value
    res = na.omit(res)
    res = res[order(res$padj),]

    # Normalize count and visualize the effect by read count before / after normalization
    raw_counts = counts(dds,normalized=FALSE)
    norm_counts = counts(dds,normalized=TRUE)

    # Dispersion Plot of Fit
    cat("  Create a dispersion plot of fit\n")
    svg (filename = paste(prefix, "Dispersion_Plot.svg", collapse="_"), width = 10, height = 10, pointsize = 12)
    plotDispEsts(dds, main = paste(group, "- Dispersion_Plot", collapse=" "))
    dev.off()

#    # Write the significantly differentially expressed genes and the group level
#    cat("  Extract significant genes and write a table\n")
#    significant = res[abs(res$log2FoldChange) >= 1 & res$padj <= 0.05 & !is.na(res$padj),]
#    sub_significant = res[abs(res$log2FoldChange) < 1 & abs(res$log2FoldChange) >= 0.584963 & res$padj > 0.05 & !is.na(res$padj) & res$padj <= 0.1 & !is.na(res$padj) , ]

#    table_name = paste(prefix, "overall_significant_DE.csv", collapse="_")
#    write("SIGNIFICANT (Fold Change >= 2 AND ajusted pValue <= 0.05)\n", file = table_name)
#    if(nrow(significant)>0) {
#        write.table(as.data.frame(significant), file = table_name, sep="\t", col.names = T, append=T)
#    }
#    write("\nSUB SIGNIFICANT (Fold Change >= 1.5 AND ajusted pValue <= 0.1)", file = table_name, append=T)
#    if(nrow(sub_significant)>0) {
#        write.table(as.data.frame(sub_significant), file = table_name, sep="\t", col.names = T, append=T)
#    }
#    write("\nALL", file = table_name, append=T)
#    write.table(as.data.frame(res), file = table_name, sep="\t", col.names = T, append=T)

#    ## MA Plot => log2 fold changes attributable to a given variable over the mean of normalized counts
#    cat("  Create a MA plot\n")
#    svg (filename = paste(prefix, "MA_Plot.svg", collapse="_"), width = 10, height = 10, pointsize = 12)
#    plotMA(res, main = paste(group, "- MA_Plot", collapse=" "), ylim=c(-5,5))
#    if(nrow(significant)>0) {
#        text (significant$baseMean, significant$log2FoldChange, labels=rownames(significant), col="tomato4")
#    }
#    if(nrow(sub_significant)>0) {
#        text (sub_significant$baseMean, sub_significant$log2FoldChange, labels=rownames(sub_significant), col="tomato")
#    }
#    abline(h=c(-1,1),lty=2,col="grey")
#    dev.off()

    # Regularized log transformation
    # The function rlog, stands for regularized log, transforming the original count data to the log2 scale by fitting
    # a model with a term for each sample and a prior distribution on the coefficients which is estimated from the
    # data
    rld = rlog(dds, fitType="mean")
    rlogMat = assay(rld)
    # Remove null values
    rlogMat = rlogMat[rowSums(rlogMat) > 0, ]

    # Heatmap showing the expression data of the 30 most highly expressed genes from regularized log transformatio data
    cat("  Create a heatmap of the 30 most highly expressed genes\n")
    select <- order(rowMeans(rlogMat), decreasing=TRUE)[1:30]
    svg (filename = paste(prefix, "Heatmap_count.svg", collapse="_"), width = 10, height = 10, pointsize = 10)
    heatmap.2(
        rlogMat[select,],
        cexCol = 2, cexRow = 1,
        main = paste(group, "- 30 most highly expressed genes", collapse=" "),
        Rowv = FALSE,
        Colv = FALSE,
        col = colorRampPalette(c("royalblue4","gray95","red4")),
        scale="none",
        dendrogram="none",
        trace="none",
        density.info="none",
        margin=c(15, 15))
    dev.off()

    # Barplot of unnormalized and normalized counts
    cat("  Create a Barplot showing the effect of the normalization\n")
    svg (filename = paste(prefix, "Read count.svg", collapse="_"), width = 15, height = 8, pointsize = 12)
    par(mfrow=c(1,3), mar=c(10, 6, 4, 0))
    barplot(apply(raw_counts,2,sum), col=cond_colors, las=2, main="Read count not normalized")
    barplot(apply(norm_counts,2,sum), col=cond_colors, las=2, main="Read count Normalized")
    barplot(apply(rlogMat,2,sum), col=cond_colors, las=2, main="Read count log Normalized")
    dev.off()

    # If enough samples in the group only
    if (nind >= 3) {

        # Heatmap of the sample-to-sample correlation
        distsRL = dist (t(rlogMat))
        cat("  Create a heatmap of the sample-to-sample correlation\n")
        svg (filename = paste(prefix, "Heatmap_sample_cor.svg", collapse="_"), width = 10, height = 10, pointsize = 10)
        heatmap.2(
            as.matrix(distsRL),
            cexCol = 2, cexRow = 2,
            main = paste(group, "- Sample correlation", collapse=" "),
            Rowv=as.dendrogram(hclust(distsRL)),
            symm=TRUE,
            trace="none",
            density.info="none",
            col = colorRampPalette(c("gray90","gray10")),
            margin=c(15, 15),
            ColSideColors=cond_colors,
            RowSideColors=cond_colors)
        dev.off()

        # Perform PCA with FactoMineR package
        cat("  Perform a principal component analysis\n")
        pca = PCA(t(rlogMat), scale.unit=T, graph = F) # transpose matrix for correct pca orientation

        svg (filename = paste(prefix, "PCA_Samples.svg", collapse="_"), width = 10, height = 10, pointsize = 12)
        plot(pca, axes = c(1, 2), choix = "ind", title = paste(group, "- PCA Samples", collapse=" "), col.ind = cond_colors)
        dev.off()

        svg (filename = paste(prefix, "PCA_genes.svg", collapse="_"), width = 10, height = 10, pointsize = 12)
        plot(pca, axes = c(1, 2), choix = "var", title = paste(group, "- PCA Genes", collapse=" "), cex = 0.5, lim.cos2.var = 0.75)
        dev.off()

        Dim = dimdesc(pca, axes = c(1, 2))
        output = rbind ( c("Dim1.Correlation", "Dim1.p-Value"), Dim$Dim.1$quanti, c("Dim2.Correlation", "Dim2.p-Value"), Dim$Dim.2$quanti)
        write.table(output, paste(prefix, "_PCA_var_dim.csv", collapse="_"), sep="\t", col.names = F)
    }

    # Compare all conditions pair per pair
    combi = t(combn(levels(conds),2))

    for (i in 1:nrow(combi)) {
        cond1 = combi[i,1]
        cond2 = combi[i,2]

        cat(paste("  Comparing conditions", cond1, "and", cond2, "\n", collapse=" "))

        # Results only for the 2 considered conditions
        sub_res = results(dds, contrast=c("Condition",cond1, cond2))

        # Remove NA and sort by adjusted p value
        sub_res = na.omit(sub_res)
        sub_res = sub_res[order(sub_res$padj),]

        # Write the significantly differentially expressed genes and the pair level
        cat("    Extract significant genes and write a table\n")
        significant = sub_res[abs(sub_res$log2FoldChange) >= 1 & sub_res$padj <= 0.05 & !is.na(sub_res$padj),]
        #cat ("SIGNIFICANT: ", rownames(significant), "\n")
        sub_significant = sub_res [abs(sub_res$log2FoldChange) < 1 & abs(sub_res$log2FoldChange) >= 0.584963 & sub_res$padj > 0.05 & sub_res$padj <= 0.1, ]
        #cat ("SUB SIGNIFICANT: ", rownames(sub_significant), "\n")

        table_name = paste(prefix, cond1, "vs", cond2, "significant_DE.csv", collapse="_")
        write("SIGNIFICANT (Fold Change >= 2 AND ajusted pValue <= 0.05)\n", file = table_name)
        if(nrow(significant)>0) {
            write.table(as.data.frame(significant), file = table_name, sep="\t", col.names = T, append=T)
        }
        write("\nSUB SIGNIFICANT (Fold Change >= 1.5 AND ajusted pValue <= 0.1)", file = table_name, append=T)
        if(nrow(sub_significant)>0) {
            write.table(as.data.frame(sub_significant), file = table_name, sep="\t", col.names = T, append=T)
        }
        write("\nALL", file = table_name, append=T)
        write.table(as.data.frame(sub_res), file = table_name, sep="\t", col.names = T, append=T)

        # MA-PLOT=> log2 fold changes attributable to a given variable over the mean of normalized counts

        cat("    Create a MA plot\n")
        svg (filename = paste(prefix, "MA_plot", cond1,"vs",cond2,".svg" , collapse="_"), width = 10, height = 10, pointsize = 12)
        plotMA(sub_res, main = paste(group, "- MA-plot", cond1, "vs", cond2 , collapse=" "), ylim=c(-5,5))
        if(nrow(significant)>0) {
            text (significant$baseMean, significant$log2FoldChange, labels=rownames(significant), col="tomato4")
        }
        if(nrow(sub_significant)>0) {
            text (sub_significant$baseMean, sub_significant$log2FoldChange, labels=rownames(sub_significant), col="tomato")
        }
        abline(h=c(-1,1),lty=2,col="grey")
        dev.off()

        ## VOLCANO PLOT
        #out_of_range = significant[sub_significant$padj
        out_of_range = significant[significant$padj < 0.00000001,]
        out_of_range$padj = replace(out_of_range$padj, out_of_range$padj < 0.00000001, 0.00000001)
        print (out_of_range)

        # Look for RNAs changing significantly between conditions
        cat("    Create a volcano plot\n")
        svg (filename = paste(prefix, "Volcano_plot", cond1,"vs",cond2,".svg" , collapse="_"), width = 10, height = 10, pointsize = 12)
        plot(
            sub_res$log2FoldChange,
            -log10(sub_res$padj),
            main=paste(group, "- Volcano plot", cond1, "vs", cond2 , collapse=" "),
            xlim=c(-5,5),
            ylim=c(0,8),
            pch=20,
            cex=0.5,
            xlab="Log2 Fold change",
            ylab="Log10 Ajusted p Value")

        abline(h=-log10(0.05),lty=2,col="grey")
        abline(v=c(-1,1),lty=2,col="grey")
        legend("topleft",cond2,cex=0.5)
        legend("topright",cond1,cex=0.5)

        if(nrow(significant)>0) {
            text (significant$log2FoldChange, -log10(significant$padj), labels=rownames(significant), col="tomato4")
        }
        if(nrow(sub_significant)>0) {
            text (sub_significant$log2FoldChange, -log10(sub_significant$padj), labels=rownames(sub_significant), col="tomato")
        }
        if(nrow(out_of_range)>0) {
            text( out_of_range$log2FoldChange, -log10(out_of_range$padj), labels=rownames(out_of_range), col="tomato4")
            points (out_of_range$log2FoldChange, -log10(out_of_range$padj), col = "black", cex=0.5, pch=24)
        }

        dev.off()

    }
    return()
}

###################################################################################################
# MAIN SCRIPT
###################################################################################################

cat ("### INITIALIZE ###\n\n")

# Parse args
sample_list = args[1]
gff = args[2]

# Read in a sample data table
pdata = read.table(sample_list, header=TRUE, strip.white=TRUE)

# Analyse per Time point groups
for (group in levels(pdata$TimePoint)) {
    sub_pdata = as.data.frame(pdata[pdata$TimePoint == group, ])
    run_analysis (sub_pdata, group)
    }

# Analyse per series groups
for (group in levels(pdata$Serie)) {
    sub_pdata = as.data.frame(pdata[pdata$Serie == group, ])
    run_analysis (sub_pdata, group)
    }

#warnings()
