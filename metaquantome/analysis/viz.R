# options
options(stringsAsFactors = FALSE, message=FALSE)

# libraries
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))


# constants
grp_color_values <- c("dodgerblue", "darkorange2")

# ns_ws_colors <- rep(ns_ws_color_values, 11)

# from brewer.pal(name = 'PuBu', n = 9)
heatmap_colors <- c("#FFF7FB", "#ECE7F2", "#D0D1E6", "#A6BDDB", "#74A9CF",
                    "#3690C0", "#0570B0", "#045A8D", "#023858")

# read in result table
read_result <- function(file){
    df <- read.delim(file, sep="\t", stringsAsFactors=FALSE)
    return(df)
}

df <- data.frame(taxon_name = c('a', 'b', 'c'),
                      NS_mean = c(1, 10, 5), stringsAsFactors = FALSE)

# barplot function
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

# mq_barplot('test.png', grpmean = "NS_mean", df=df, ntaxa=5)

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

# heatmap function

# clust sep function

# pca function

# volcano
volcano_colors <- scale_color_manual(values = c("grey50", "seagreen3"), guide = FALSE)


# main arg parsing
main <- function(){
    args <- commandArgs(trailingOnly=TRUE)
    plttype <- args[1]
    if (plttype == "bar"){
        barplot_cli(args)
    }
}

main()
