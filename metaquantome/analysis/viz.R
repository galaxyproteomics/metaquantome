# options
options(stringsAsFactors = FALSE, message=FALSE)

####### ==================== #######
#              LIBRARIES           #
####### ==================== #######
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(gplots))
suppressMessages(library(jsonlite))

####### ==================== #######
#         UTILITY FUNCTIONS        #
####### ==================== #######
read_result <- function(file){
    df <- read.delim(file, sep="\t", stringsAsFactors=FALSE)
    return(df)
}

impute <- function(mat) {
    mat[mat == 0] <- NA
    mat[is.na(mat)] <- min(mat, na.rm = TRUE) * 1e-3
    mat
}

pad <- function(data, multiple){
    mind <- min(data)
    maxd <- max(data)
    if (mind >= 0){
        pmind <- mind - mind * (multiple - 1)
    } else {
        pmind <- mind + mind * (multiple - 1)
    }
    if (maxd >= 0){
        pmaxd <- maxd * multiple
    } else {
        pmaxd <- maxd - maxd * (multiple - 1)
    }
    c(pmind, pmaxd)
}

####### ==================== #######
#              BARPLOT             #
####### ==================== #######
mq_barplot <- function(df, img, mode, meancol,
                       nterms, width, height, target_rank){
    if (!(meancol %in% names(df))){
        stop('Mean column name not found in dataframe. Check spelling and try again.',
             call. = FALSE)
    }
    # filter to desired rank, if taxonomy
    if (mode == "tax"){
        # todo: add check for target rank in
        # list of available target ranks
        df <- df[df[, "rank"] == target_rank, ]
    }
    png(img, height=height, width=width, units="in", res=300)
    df[, meancol] <- 2^df[, meancol]
    reord <- df[order(df[, meancol], decreasing=TRUE), ]
    sub_reord <- reord[1:nterms, ]
    if (mode == "tax"){
        barnames = sub_reord[, "taxon_name"]
    } else {
        barnames = sub_reord[, "name"]
    }
    barplot(names.arg = barnames,
            height = sub_reord[, meancol],
            las = 1, cex.names = 0.5,
            xlab = "Taxon",
            ylab = "Total Peptide Intensity")
    grid()
    # send the message to the ether
    ether <- dev.off()
}

barplot_cli <- function(args){
    # args are as follows
    # 1. plot type (guaranteed to be 'bar')
    # 2. pltfile - output image file
    # 3. input tabular file
    # 4. mode - tax or fn
    # 5. name of mean column
    # 6. number of taxa
    # 7. image width (default 5)
    # 8. image height (default 5)
    # 9. target rank (taxonomy only)
    img <- args[2]
    infile <- args[3]
    mode <- args[4]
    meancol <- args[5]
    nterms <- args[6]
    width <- as.numeric(args[7])
    height <- as.numeric(args[8])
    target_rank <- args[9]
    df <- read.delim(infile, sep="\t", stringsAsFactors=FALSE)
    plt <- mq_barplot(df, img=img, mode=mode,
                      meancol=meancol,
                      nterms=nterms,
                      height=height, width=width,
                      target_rank=target_rank)
}

####### ==================== #######
#              HEATMAP             #
####### ==================== #######
cor.dist <- function(x){
    as.dist(1-cor(t(x)))
}
hclust.ward <- function(x) {
    hclust(x,method="ward.D")
}

grp_color_values <- c("dodgerblue", "darkorange", "yellow2", "red2", "darkviolet", "black")


# from brewer.pal(name = 'PuBu', n = 9)
heatmap_colors <- c("#FFF7FB", "#ECE7F2", "#D0D1E6", "#A6BDDB", "#74A9CF",
                    "#3690C0", "#0570B0", "#045A8D", "#023858")

mq_heatmap <- function(img, df, all_intcols, colSideColors, filter_to_sig, alpha, width, height){
    # df is the output from either expand, stat, or filter
    # samp_columns is a vector of all columns with the term intensities
    # colSide colors is a vector of colors for the groups. Must be in the same order as samp_columns
    # filter to sig
    if (filter_to_sig){
        print(df)
        df <- df[df$corrected_p < alpha, ]
        print(df)
    }

    # impute
    mat <- impute(data.matrix(df[, all_intcols]))

    # scale rows
    datmat.scale <- t(apply(mat, 1, scale))
    rownames(datmat.scale) <- df$id
    colnames(datmat.scale) <- colnames(df[, all_intcols])
    par(mar = rep(5, 4))
    feature.dend <- as.dendrogram(hclust.ward(cor.dist(datmat.scale)))
    sample.dend <- as.dendrogram(hclust.ward(cor.dist(t(datmat.scale))))
    png(img, width=width, height=height, res=500, units="in")
    heatmap.2(datmat.scale,
              Colv = sample.dend,
              Rowv = feature.dend,
              distfun = cor.dist,
              hclustfun = hclust.ward,
              trace="none",
              col = heatmap_colors,
              margins = c(10, 10),
              cexRow = 0.3,
              ColSideColors = colSideColors,
              density.info = "none")
    ether <- dev.off()
}

heatmap_cli <- function(args){
    # 1. plot type (guaranteed to be 'heatmap')
    # 2. pltfile - output image file
    # 3. input tabular file
    # 4. all intcols from SampleGroups, as comma-separated list
    # 5. json dump from SampleGroups() top-level dictionary
    # 6. filter to sig ("True" if TRUE)
    # 7. alpha - significance level
    # 8. image width (default 5)
    # 9. image height (default 5)

    img <- args[2]
    infile <- args[3]
    df <- read.delim(infile, sep="\t", stringsAsFactors=FALSE)

    # split all_intcols from SampleGroups(), for samp_columns vector
    all_intcols <- unlist(strsplit(args[4], ","))
    nsamps <- length(all_intcols)

    # read in json dump from SampleGroups() basic dictionary
    # will be a list
    jsonDump <- fromJSON(args[5])
    grps <- names(jsonDump)
    ngrps <- length(grps)

    # create grps to color mapping
    grp_col_mapping <- grp_color_values[1:ngrps]
    names(grp_col_mapping) <- grps

    # we need sample -> group mapping
    colSideColors <- rep(0, nsamps)
    for (i in 1:nsamps){
        this_intcol <- all_intcols[i]
        for (j in 1:ngrps){
            samps_in_grp <- jsonDump[[grps[j]]]
            if (this_intcol %in% samps_in_grp){
                colSideColors[i] <- grp_col_mapping[grps[j]]
            }
        }
    }
    filter_to_sig <- (args[6] == "True")
    alpha <- as.numeric(args[7])
    width <- as.numeric(args[8])
    height <- as.numeric(args[9])
    mq_heatmap(img, df, all_intcols, colSideColors, filter_to_sig, alpha, width, height)
}

####### ==================== #######
#              PCA             #
####### ==================== #######



####### ==================== #######
#              VOLCANO             #
####### ==================== #######
mq_volcano <- function(df, img, fc_name, width, height, textannot, gosplit){
    # df is the dataframe after stat
    # fc_name is the name of the column with the fold change data
    # textcol is the name of the column with the text describing the term
    df$fc <- df[, fc_name]
    df$neglog10p <- -log10(df[, "corrected_p"])
    df$de <- abs(df$fc) > 1 & df$corrected_p < 0.05
    xmax <- max(df$fc) * 1.2
    xmin <- min(df$fc) * 1.2
    ymax <- max(df$neglog10p) * 1.2
    ymin <- 0
    volcano_colors <- scale_color_manual(values = c("grey50", "seagreen3"), guide=FALSE)
    if (gosplit){
        vplt <- ggplot(df, aes(x = fc, y = neglog10p)) +
            geom_point(aes(color = de)) +
            facet_grid(.~namespace)
    } else {
        vplt <- ggplot(df, aes(x = fc, y = neglog10p)) +
            geom_point(aes(color = de))
    }
    if (!(textannot=="None")){
        df$id <- df[, textannot]
        vplt <- vplt +
            geom_text(data = subset(df, de),
                aes(label = id), check_overlap = TRUE, nudge_y = 0.15,
                    size = 3)
    }
    vplt <- vplt +
        geom_vline(xintercept = -1, lty = 2, alpha = 0.7) +
        geom_vline(xintercept = 1, lty = 2, alpha = 0.7) +
        geom_hline(yintercept = -log10(0.05), lty = 3) +
        theme_bw(base_size=16) +
        volcano_colors +
        scale_x_continuous(limits = c(xmin, xmax)) +
        scale_y_continuous(limits = c(ymin, ymax)) +
        labs(x = "Log2 Fold Change", y = "-Log10 P Value")
    vplt
    ggsave(file=img, width=width, height=height, units="in", dpi=500)
}

volcano_cli <- function(args){
    # args are as follows
    # 1. plot type (guaranteed to be 'volcano')
    # 2. pltfile - output image file
    # 3. input tabular file
    # 4. name of text annotation column
    # 5. name of fold change column
    # 6. whether to split GO by ontology/namespace
    # 7. image width (default 5)
    # 8. image height (default 5)
    img <- args[2]
    infile <- args[3]
    textannot <- args[4]
    fc_name <- args[5]
    gosplit <- (args[6] == "True")
    width <- as.numeric(args[7])
    height <- as.numeric(args[8])
    df <- read.delim(infile, sep="\t", stringsAsFactors=FALSE)
    plt <- mq_volcano(df, img=img, textannot=textannot, fc_name=fc_name, gosplit=gosplit, width=width, height=height)
}


####### ==================== #######
#              MAIN                #
####### ==================== #######
main <- function(){
    args <- commandArgs(trailingOnly=TRUE)
    plttype <- args[1]
    if (plttype == "bar"){
        barplot_cli(args)
    }
    if (plttype == "volcano"){
        volcano_cli(args)
    }
    if (plttype == "heatmap"){
        heatmap_cli(args)
    }
}

main()
