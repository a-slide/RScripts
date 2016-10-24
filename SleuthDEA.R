#!/usr/bin/env Rscript
options( warn = -1 )

###################################################################################################
# FUNCTIONS
###################################################################################################

sleuth_analysis = function(s2c, t2g, full_model, outdir, max_bootstrap, norm_fun, filter_fun, aggregation_column, sig_level, wald, wald_plot,
    lrt, var_graph, rds, count_table, fit_table) {
    
    cat("Extract the data from the kalisto results and normalise\n")
    so = sleuth_prep(s2c, full_model, filter_fun, t2g, max_bootstrap, norm_fun, norm_fun, aggregation_column)
    
    cat("Estimate parameters for the sleuth response error measurement (full) model\n")
    so = sleuth_fit(so)
    
    cat("Estimate parameters for the sleuth reduced model (shrinkage)\n")
    so = sleuth_fit(so, ~1, 'reduced')
    
    ## Wald tests    
    if (wald == TRUE){
        cat("Performing Wald test for all cofactors\n")
        plots = list()

        for (beta_factor in colnames(so$design_matrix)) {
            if (beta_factor != "(Intercept)") {
                
                cat(paste("\tAnalysing cofactor:", beta_factor, "\n"))
                so = sleuth_wt(so, which_beta=beta_factor, which_model="full")
                
                cat("\t\tGenerate tables\n")
                # Extract and save the results in a table
                res = sleuth_results(so, beta_factor, test_type="wt", which_model="full", show_all=FALSE)
                fname = paste0(outdir, beta_factor, "_wald_test.tsv")
                write.table(res, file = fname, row.names=TRUE, na="",col.names=TRUE, sep="\t", quote=FALSE)
                
                # Count number of significant points
                cat(paste("\t\tSignificant genes: ", nrow(filter(res, qval <= sig_level)), "\n"))
                cat(paste("\t\tSignificant genes enriched: ", nrow(filter(res, qval<= sig_level & b<0)), "\n"))
                cat(paste("\t\tSignificant genes depleted: ", nrow(filter(res, qval<= sig_level & b>0)), "\n"))
                
                if (wald_plot ==TRUE) {
                    cat("\t\tGenerate plots\n")
                    # Extract the ERCC results if available
                    ERCC_res = filter(res, grepl("ERCC",target_id) & !is.na(b) & !is.na(qval))
                    ERCC_id = select(ERCC_res, target_id)
                    
                    # Volcano plot
                    p = plot_volcano(so, beta_factor, sig_level=sig_level)
                    p = p+ggtitle(paste("Volcano Plot, Wald test, Full Model, FDR" ,sig_level ,"Condition", beta_factor))
                    p = p+geom_point(data=ERCC_res, aes(b, -log10(qval)), colour="lightgreen")
                    plots =  append(plots, list(p))
                    
                    # MA plot
                    p = plot_ma(so, beta_factor, sig_level=sig_level, highlight=ERCC_id, highlight_color = "lightgreen")
                    p = p+ggtitle(paste("MA Plot, Wald test, Full Model, FDR" ,sig_level ,"Condition", beta_factor))
                    plots =  append(plots, list(p))
                    
                    # QQ plot
                    p = plot_qq(so, beta_factor, sig_level=sig_level)
                    p = p+ggtitle(paste("MA Plot, Wald test, Full Model, FDR" ,sig_level ,"Condition", beta_factor))
                    plots =  append(plots, list(p))
                }
            }
        }
        if (wald_plot ==TRUE) {
            cat("\tExport plots\n")
            # Plot all the wald test graphs in the same file
            gs = marrangeGrob(grobs=plots, nrow=3, ncol=1)
            ggsave(paste0(outdir, "Wald_DE_plot.pdf"), gs, width=15, height=15)
        }
    }
    
    ## Likelyhood ratio tests
    if (lrt == TRUE){
        cat("Perform likelyhood ratio testing\n")
        so = sleuth_lrt(so, 'reduced', 'full')
        
        cat("\t\tExport tables\n")
        # Save all the results in a table
        res = sleuth_results(so, 'reduced:full', test_type = 'lrt')
        fname = paste0(outdir, "likelyhood_ratio_test.tsv")
        write.table(res, file = fname, row.names=TRUE, na="",col.names=TRUE, sep="\t", quote=FALSE)
        
        # Count number of significant points
        cat(paste("\t\tSignificant genes: ", nrow(filter(res, qval<= sig_level)), "\n"))
    }
    
    ## Variance graphs
    if (var_graph == TRUE){
        cat("Plot variance graphs\n")
        # plot PCA graph
        p1 = plot_pca(so, text_labels=T, use_filtered = T, units = 'tpm')
        p1 = p1+ggtitle(paste("PCA plot, PC1/PC2, Unit tpm"))
        # plot principal component contribution to the variance
        p2 = plot_pc_variance(so, scale=T, use_filtered = T, units = 'tpm', )
        p2 = p2+ggtitle(paste("% contribution to PC variances, Unit tpm, Scaled, Filtered values"))
        # plot the mean variance graph with fitting line
        p3 = plot_mean_var(so)
        p3 = p3+ggtitle(paste("Mean Variance plot"))
        
        # Plot all the graphs in the same file
        lay = rbind(c(1,1,1,NA,NA), c(1,1,1,2,2), c(1,1,1,2,2), c(3,3,3,3,3), c(3,3,3,3,3))
        gs = grid.arrange(p1, p2, p3, layout_matrix = lay)
        ggsave(paste0(outdir, "variance_plots.pdf"), gs, width=15, height=15)
    }
    
    # Table of normalized counts
    if (count_table == TRUE){
        cat("Export count table\n")
        fname = paste0(outdir, "obs_norm.tsv")
        write.table(so$obs_norm, file = fname, row.names=TRUE, na="",col.names=TRUE, sep="\t", quote=FALSE)
    }
    
    # Table of fit parameters counts
    if (fit_table == TRUE) {
        cat("Export fit table\n")
        fname = paste0(outdir, "obs_norm_filt.tsv")
        write.table(so$fits[["full"]]$summary, file = fname, row.names=TRUE, na="",col.names=TRUE, sep="\t", quote=FALSE)
    }
    
    #Save the Sleuth object for further analysis if requested
    if (rds == TRUE) {
        cat("Export Sleuth rds results\n")
        fname = paste0(outdir, "results.rds")
        saveRDS(so, file = fname)
    }
    return(so)
}

###################################################################################################
# MAIN SCRIPT
###################################################################################################

# Lib import
options(width = 160) # wider terminal size
suppressWarnings(suppressMessages(library("dplyr")))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("methods"))
suppressWarnings(suppressMessages(library("sleuth")))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))

cat ("\n### INITIALIZE ###\n")

# Define options for the command line arguments
option_list = list(
    make_option(c("--s2c"), type="character", metavar="character", default=NULL, 
        help="Sample list file containing the sample name (colname = sample), grouping covariate column(s) and the path (colname = path) to the kallisto dir containing results [mandatory argument]"),
    make_option(c("--t2g"), type="character", metavar="character", default=NULL,
        help="File containing the transcript (colname = target_id) to gene (ens_gene) information [default \"%default\"]"),
    make_option(c("-o", "--outdir"), type="character", metavar="character", default="./out/",
        help="Directory where to output the result files [default \"%default\"]"),
    make_option(c("-c", "--select_col"), type="character", metavar="character", default=NULL,
        help="if needed will select row with *select_val* from this column [default \"%default\"]"),
    make_option(c("-v", "--select_val"), type="character", metavar="character", default=NULL,
        help="if needed will select row with this value from *select_col* column [default \"%default\"]"),
    make_option(c("-m", "--full_model"), type="character", metavar="character", default="1",
        help="Define the model to fit the data based on the covariates. Example: cell+localization [default \"%default\"]"),
    make_option(c("-b", "--max_bootstrap"), type="integer", metavar="number", default=NULL,
        help="Max number of bootstrap to analyse [default All, Minimun 2]"),
    make_option(c("-n", "--no_norm"), action="store_true", default=FALSE,
        help="If given, will deactivate the default geometric mean normalisation [default \"%default\"]"),
    make_option(c("-e", "--min_est"), type="integer", metavar="number", default=5,
        help="Minimal number of estimated counts per transcript for the filtering step [default \"%default\"]"),
    make_option(c("-p", "--min_prop"), type="double", metavar="number", default=0.47,
        help="Minimal proportion of sample with *min_est* counts for the filtering step [default \"%default\"]"),
    make_option(c("-a", "--aggreg_col"), action="store_true", default=FALSE,
        help=" [If given, will aggregate genes based on the *ens_gene* column of the file indicated by the t2g arg [default \"%default\"]"),
    make_option(c("-t", "--threads"), type="integer", metavar="number", default=1,
        help="Number of core to be used (for gene level aggregation only) [default \"%default\"]"),
    make_option(c("-s", "--sig_level"), type="double", metavar="number", default=0.1,
        help="FDR threshold fo significance in the plots [default \"%default\"]"),
    make_option(c("-w", "--wald"), action="store_true", default=FALSE,
        help="If given, will perform differential expression tests with wald method [default \"%default\"]"),
    make_option(c("--wald_plot"), action="store_true", default=FALSE,
        help="If given, will plot graphs associated with wald test (volcano, MA and QQ plots) [default \"%default\"]"),
    make_option(c("-l", "--lrt"), action="store_true", default=FALSE,
        help="If given, will perform differential expression tests with likelyhood method (no graph) [default \"%default\"]"),
    make_option(c("-g", "--var_graph"), action="store_true", default=FALSE,
        help="If given, will plot variance graphs (PCA, Mean variance plot) [default \"%default\"]"),
    make_option(c("-r", "--rds"), action="store_true", default=FALSE,
        help="If given, output an rds file for further analysis with sleuth_live (takes time) [default \"%default\"]"),
    make_option(c("--count_table"), action="store_true", default=FALSE,
        help="If given, output tsv tables containing the normalized + filtered counts [default \"%default\"]"),
    make_option(c("--fit_table"), action="store_true", default=FALSE,
        help="If given, will output a table summarizing the fit performed by sleuth [default \"%default\"]"))

# Define usage string
usage = "SleuthDEA.T --s2c Sample_Sheet_file.tsv [--t2g transcript_to_gene.tsv, ...]"

# Parse Arguments
opt_parser = OptionParser(usage=usage, option_list=option_list)
opt = parse_args(opt_parser)

cat("Parse arguments and import files\n")

# Verify that the sample to covariate file is given and valid
if (!is.null(opt$s2c)){
    if (file.access(opt$s2c)== 0) {
        # Import the file and select row if needed
        s2c = read.table(opt$s2c, header=T, stringsAsFactors=F)
        if (!is.null(opt$select_col)) { s2c = s2c[s2c[,opt$select_col] == opt$select_val,] }
        cat("Select sample:\n")
        cat (paste("\t", s2c$sample, "\n"))
    }else {
        print_help(opt_parser)
        stop("Invalid s2c file\n", call.=FALSE)
    }
} else {
    print_help(opt_parser)
    stop("Please provide a sample to covariate tabulated file with kallisto directory result path(s)\n", call.=FALSE)
}

# Verify that the sample to covariate file is given and valid
if (!is.null(opt$t2g)){
    if (file.access(opt$t2g) == 0) {
        # Import the gene info in table and rename the fields for Sleuth compatibility
        t2g = read.table(opt$t2g, header=T, stringsAsFactors=F)
        t2g = select(t2g, target_id, ens_gene, transcript_name, gene_name)
    }
} else {
    cat("No gene to transcript information provided or invalid file\n")
    t2g = NULL
}

# Create the output directory
dir.create(opt$outdir, showWarnings=F, recursive=T)
# Convert the model to a formula
full_model = as.formula(paste0("~", opt$full_model))
# Define the norm function or deactivate it
if (opt$no_norm == TRUE) {
    norm_fun = function(mat) { setNames(rep(1.0,ncol(mat)),colnames(mat)) }
} else {
    norm_fun=norm_factors
}
# Define a filter function
filter_fun = function(row, ...) {mean(row>=opt$min_est)>=opt$min_prop}
# Redefine aggregation to Null if no transcript to gene file was provided
if (is.null(opt$t2g) | opt$aggreg_col == FALSE) {
    aggregation_column = NULL
}else {
    aggregation_column = "ens_gene"
}
# Multicore option
options(mc.cores = opt$threads)

cat ("\n### SLEUTH ANALYSIS ###\n")
# Launch the Sleuth Pipeline
so = sleuth_analysis(s2c, t2g, full_model, opt$outdir, opt$max_bootstrap, norm_fun, filter_fun, aggregation_column, opt$sig_level,
     opt$wald, opt$wald_plot, opt$lrt, opt$var_graph, opt$rds, opt$count_table, opt$fit_table)

cat ("\n### DONE ###\n")
