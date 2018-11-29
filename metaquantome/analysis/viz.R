# options
options(stringsAsFactors = FALSE, message=FALSE, warnings=FALSE)

####### ==================== #######
#              LIBRARIES           #
####### ==================== #######
suppressWarnings(suppressMessages(library(ggplot2)))
suppressMessages(library(gplots))
suppressMessages(library(jsonlite))
suppressMessages(library(stringr))

####### ==================== #######
#             CONSTANTS            #
####### ==================== #######
grp_color_values <- c("dodgerblue", "darkorange",
    "yellow2", "red2", "darkviolet", "black")

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

get_all_intcols <- function(str_intcols){
   # split all_intcols from SampleGroups(), for samp_columns vector
    all_intcols <- unlist(strsplit(str_intcols, ","))
    return(all_intcols)
}

get_colors_from_groups <- function(json_dump, all_intcols){
    # read in json dump from SampleGroups() basic dictionary
    # will be a list
    grp_list <- fromJSON(json_dump)
    grps <- sort(names(grp_list))
    ngrps <- length(grps)

    # create grps to color mapping
    grp_col_mapping <- grp_color_values[1:ngrps]
    names(grp_col_mapping) <- grps

    # we need sample -> group mapping
    nsamps <- length(all_intcols)
    colSideColors <- rep(0, nsamps)
    for (i in 1:nsamps){
        this_intcol <- all_intcols[i]
        for (j in 1:ngrps){
            samps_in_grp <- grp_list[[grps[j]]]
            if (this_intcol %in% samps_in_grp){
                colSideColors[i] <- grp_col_mapping[grps[j]]
            }
        }
    }
    return(colSideColors)
}


####### ==================== #######
#              BARPLOT             #
####### ==================== #######
mq_barplot <- function(df, img, mode, meancol,
                       nterms, width, height, target_rank, int_barcol){
    if (!(meancol %in% names(df))){
        stop('Mean column name not found in dataframe. Check spelling and try again.',
             call. = FALSE)
    }
    # filter to desired rank, if taxonomy
    if (mode == "t"){
        # todo: add check for target rank in
        # list of available target ranks
        df <- df[df[, "rank"] == target_rank, ]
    }

    df[, meancol] <- 2^df[, meancol]
    reord <- df[order(df[, meancol], decreasing=TRUE), ]
    sub_reord <- reord[1:nterms, ]
    if (mode == "t"){
        barnames = sub_reord[, "taxon_name"]
    } else {
        barnames = sub_reord[, "name"]
    }
    # color mapping
    barcol <- grp_color_values[int_barcol]
    # plot
    png(file=img, height=height, width=width, units="in", res=500)
    par(mar=c(5.1, 6.1, 4.1, 2.1))
    barplot(names.arg = barnames, col=barcol,
            height = sub_reord[, meancol],
            las = 1, cex.names = 0.7,
            xlab = "Taxon",
            ylab = "")
    mtext(text="Total Peptide Intensity", side=2, line=5)
    grid()
    # send the message to the ether
    ether <- dev.off()
}

barplot_cli <- function(args){
    # args are as follows
    # 1. plot type (guaranteed to be 'bar')
    # 2. pltfile - output image file
    # 3. input tabular file
    # 4. mode - t or f
    # 5. name of mean column
    # 6. number of taxa
    # 7. image width (default 5)
    # 8. image height (default 5)
    # 9. target rank (taxonomy only)
    # 10. bar color - integer from 1 to 6
    img <- args[2]
    infile <- args[3]
    df <- read_result(infile)
    mode <- args[4]
    meancol <- args[5]
    nterms <- args[6]
    width <- as.numeric(args[7])
    height <- as.numeric(args[8])
    target_rank <- args[9]
    barcol <- as.numeric(args[10])
    plt <- mq_barplot(df, img=img, mode=mode,
                      meancol=meancol,
                      nterms=nterms,
                      height=height, width=width,
                      target_rank=target_rank, int_barcol=barcol)
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


# from brewer.pal(name = 'PuBu', n = 9)
heatmap_colors <- c("#FFF7FB", "#ECE7F2", "#D0D1E6", "#A6BDDB", "#74A9CF",
                    "#3690C0", "#0570B0", "#045A8D", "#023858")

mq_heatmap <- function(img, df, all_intcols, colSideColors, filter_to_sig, alpha, width, height, strip){
    # df is the output from either expand, stat, or filter
    # samp_columns is a vector of all columns with the term intensities
    # colSide colors is a vector of colors for the groups. Must be in the same order as samp_columns
    # filter to sig
    if (filter_to_sig){
        df <- df[df$corrected_p < alpha, ]
    }

    # impute
    mat <- impute(data.matrix(df[, all_intcols]))

    # scale rows
    datmat.scale <- t(apply(mat, 1, scale))
    rownames(datmat.scale) <- df$id
    
    if (strip != "None"){
    	colnames(datmat.scale) <- str_replace(colnames(df[, all_intcols]), strip, "")
    } else {
    	colnames(datmat.scale) <- colnames(df[, all_intcols])
    }

    feature.dend <- as.dendrogram(hclust.ward(cor.dist(datmat.scale)))
    sample.dend <- as.dendrogram(hclust.ward(cor.dist(t(datmat.scale))))

    png(filename=img, width=width, height=height, res=500, units="in")
    par(mar = rep(5, 4))
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
    df <- read_result(infile)

    # split all_intcols from SampleGroups(), for samp_columns vector
    all_intcols <- get_all_intcols(args[4])

    # get color mapping
    colSideColors <- get_colors_from_groups(args[5], all_intcols)
    filter_to_sig <- (args[6] == "True")
    alpha <- as.numeric(args[7])
    width <- as.numeric(args[8])
    height <- as.numeric(args[9])
    strip <- args[10]
    mq_heatmap(img, df, all_intcols, colSideColors, filter_to_sig, alpha, width, height, strip)
}

####### ==================== #######
#                PCA               #
####### ==================== #######

make_df_list_from_pca <- function(pca_res, ind_list){
    # pca_res is the result of prcomp
    # ind_list is a list of indices for each group
    pointmat <- data.frame(pca_res$rotation[, 1:2])
    grps <- lapply(ind_list, function(i) pointmat[i, ])
    return(grps)
}

sep_n <- function(clust){
    # clust is a list of dataframes, where the columns in each dataframe are the
    # intensities for each experimental condition
    nclust <- length(clust)
    means <- lapply(clust, colMeans)
    possible_dists <- combn(1:nclust, 2)
    n_possible_dists <- ncol(possible_dists)
    dists <- rep(0, n_possible_dists)
    for (i in 1:n_possible_dists){
        comb <- possible_dists[, i]
        dists[i] <- sum((means[[comb[1]]] - means[[comb[2]]])^2)
    }
    avg_dist <- mean(dists)
    within_variance <- sapply(1:nclust, function(i) mean((clust[[i]] - means[[i]])^2))
    avg_dist / sum(within_variance)
}

mq_prcomp <- function(img, df, all_intcols, json_dump, colors, calculate_sep, width, height, strip){
    # function(img, df, all_intcols, colSideColors, filter_to_sig, alpha, width, height)
    # df is the result from filter
    # samp_columns is a vector of all intensity column names
    # cols is a vector of group colors
    mat <- impute(data.matrix(df[, all_intcols]))
    pr <- prcomp(mat, scale=TRUE, center=TRUE)


    # calculate separation
    # make list of indices
    if (calculate_sep){
        grp_list <- fromJSON(json_dump)
        grps <- names(grp_list)
        ngrps <- length(grps)
        list_indices <- vector(length=ngrps, mode="list")
        mat_colnames <- colnames(mat)
        for (i in 1:ngrps){
            names_in_grp <- grp_list[[i]]
            list_indices[[i]] <- sapply(names_in_grp, function(g) which(mat_colnames == g))
        }
        df_list <- make_df_list_from_pca(pr, list_indices)
        sep <- sep_n(df_list)
        put_sep_in_title <- paste0("Cluster Separation: ", format(sep, digits=5))
    } else {
        put_sep_in_title <- NA
    }

    # plot
    xpadding <- 1.1
    ypadding <- 1.3
    
    if (strip != "None"){
    	point_names <- str_replace_all(colnames(mat), strip, "")
    } else {
    	point_names <- colnames(mat)
    }

    png(filename=img, width=width, height=height, res=500, units="in")
    plot(pr$rotation,
         xlab = paste("PCA1 (",
                  format(summary(pr)$importance[2, 1]*100, digits = 3), "%)",
                  sep = ""),
         ylab = paste("PCA2 (",
                  format(summary(pr)$importance[2, 2]*100, digits = 3), "%)",
                  sep = ""),
         main=put_sep_in_title,
         col = colors, pch = 20,
         xlim = pad(pr$rotation[, 1], xpadding),
         ylim = pad(pr$rotation[, 2], ypadding),
         cex = 1.5)
    text(pr$rotation,
        labels = point_names,
        col = colors,
        pos = 3,
        cex = 1)
    grid()
    ether <- dev.off()
}

prcomp_cli <- function(args){
    # 1. plot type (guaranteed to be 'pca')
    # 2. pltfile - output image file
    # 3. input tabular file
    # 4. all intcols from SampleGroups, as comma-separated list
    # 5. json dump from SampleGroups() top-level dictionary
    # 6. calculate sep
    # 7. image width (default 5)
    # 8. image height (default 5)
    img <- args[2]
    infile <- args[3]
    df <- read_result(infile)
    # split all_intcols from SampleGroups(), for samp_columns vector
    all_intcols <- get_all_intcols(args[4])
    # get color mapping
    json_dump <- args[5]
    colors <- get_colors_from_groups(json_dump, all_intcols)
    calculate_sep <- (args[6] == "True")
    width <- as.numeric(args[7])
    height <- as.numeric(args[8])
    strip <- args[9]
    mq_prcomp(img=img, df=df, all_intcols=all_intcols, json_dump=json_dump,
        colors=colors, calculate_sep=calculate_sep, width=width, height=height,
        strip=strip)
}


####### ==================== #######
#              VOLCANO             #
####### ==================== #######
mq_volcano <- function(df, img, fc_name, flip_fc, width, height, textannot, gosplit){
    # df is the dataframe after stat
    # fc_name is the name of the column with the fold change data
    # textcol is the name of the column with the text describing the term
    df$fc <- df[, fc_name]
    if (flip_fc){
        df$fc <- (-1)*df$fc
    }
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
            facet_grid(namespace~.)
    } else {
        vplt <- ggplot(df, aes(x = fc, y = neglog10p)) +
            geom_point(aes(color = de))
    }
    if (!(textannot=="None")){
        df$id <- df[, textannot]
        vplt <- vplt +
            geom_text(data = subset(df, de),
                aes(label = id), check_overlap = TRUE, nudge_y = 0.05,
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
    # 6. whether to flip fc
    # 7. whether to split GO by ontology/namespace
    # 8. image width (default 5)
    # 9. image height (default 5)
    img <- args[2]
    infile <- args[3]
    df <- read_result(infile)
    textannot <- args[4]
    fc_name <- args[5]
    flip_fc <- (args[6] == "True")
    gosplit <- (args[7] == "True")
    width <- as.numeric(args[8])
    height <- as.numeric(args[9])
    plt <- mq_volcano(df, img=img, textannot=textannot, fc_name=fc_name, flip_fc=flip_fc, gosplit=gosplit, width=width, height=height)
}

####### ==================== #######
#              FT DIST             #
####### ==================== #######
short_onto_to_long <- c("bp" = "biological_process",
                        "mf"="molecular_function",
                        "cc"="cellular_component")

mq_ft_dist <- function(df, img, whichway, name, id, meancol,
                       nterms, width, height, target_rank, target_onto,
                       int_barcol){
    if (!(meancol %in% names(df))){
        stop('Mean column name not found in dataframe. Check spelling and try again.',
             call. = FALSE)
    }
    # filter to name or id, depending on which is not missing
	if (name != "None"){
		if (whichway == "t_dist"){
			df <- df[df[, "name"] == name, ]
		}
		if (whichway == "f_dist"){
			df <- df[df[, "taxon_name"] == name, ]
		}
	}
	if (id != "None"){
		if (whichway == "t_dist"){
			df <- df[df[, "go_id"] == id, ]
		}
		if (whichway == "f_dist"){
			df <- df[df[, "tax_id"] == id, ]
		}
	}
	
    # filter to desired rank, if taxonomic distribution
    if (whichway == "t_dist"){
        # list of available target ranks
        df <- df[df[, "rank"] == target_rank, ]
    } else if (whichway == "f_dist"){
        long_onto <- short_onto_to_long[target_onto]
        df <- df[df[, "namespace"] == long_onto, ]
        # also, remove any bp, cc, or mf rows
        df <- df[!(df[, "name"] %in% short_onto_to_long), ]
    } else {
    	stop("Wrong whichway - should be t_dist or f_dist.")
    }
	
	# number of rows left
	if (nrow(df) == 0){
		stop("No rows remaining in dataframe.")
	}
    # calculate proportions
    df[, meancol] <- 2^df[, meancol]

    df$props <- df[, meancol] / sum(df[, meancol], na.rm=TRUE)
    reord <- df[order(df[, "props"], decreasing=TRUE), ]
	
    if (nterms == "all"){
        sub_reord <- reord
    } else {
        sub_reord <- reord[1:nterms, ]
    }

    if (whichway == "t_dist"){
        barnames = sub_reord[, "taxon_name"]
    } else {
        barnames = sub_reord[, "name"]
    }

    # color mapping
    barcol <- grp_color_values[int_barcol]

    # other names and stuff
    if (whichway == "t_dist"){
        xlab = "Taxon"
    }
    if (whichway == "f_dist"){
        xlab = "Functional Term"
    }

    # plot
    png(file=img, height=height, width=width, units="in", res=500)
    if (whichway == "t_dist"){
    	par(mar=c(10, 6.1, 4.1, 2.1)) # bottom, left, top and right margins
    	yline <- 8
    }
    if (whichway == "f_dist"){
    	par(mar=c(15, 6.1, 4.1, 2.1)) # bottom, left, top and right margins
    	yline <- 13
    }
   
    barplot(names.arg = barnames, col=barcol,
            height = sub_reord[, "props"],
            las = 1, cex.names = 0.7,
            xlab = "",
            ylab = "",
    		las = 2)
    mtext(text="Proportion of Peptide Intensity", side=2, line=5)
    mtext(text=xlab, side=1, line=yline)
    grid()
    # send the message to the ether
    ether <- dev.off()
}

ft_dist_cli <- function(args){
    # args are as follows
    # 1. plot type (guaranteed to be 'bar')
    # 2. pltfile - output image file
    # 3. input tabular file
    # 4. whichway - t_dist or f_dist
	# 5. name
	# 6. id
    # 7. name of mean column
    # 8. number of terms
    # 9. image width (default 5)
    # 10. image height (default 5)
    # 11. target rank (t_dist only)
	# 12. target onto (f_dist only)
    # 13. bar color - integer from 1 to 6
    # function(df, img, whichway, meancol,
    #                    nterms, width, height, target_rank, int_barcol)
    img <- args[2]
    infile <- args[3]
    df <- read_result(infile)
    whichway <- args[4]
    name <- args[5]
    id <- args[6]
    meancol <- args[7]
    nterms <- args[8]
    width <- as.numeric(args[9])
    height <- as.numeric(args[10])
    target_rank <- args[11]
    target_onto <- args[12]
    barcol <- as.numeric(args[13])
    plt <- mq_ft_dist(df, img=img, whichway=whichway,
    				  name=name, id=id,
                      meancol=meancol,
                      nterms=nterms,
                      height=height, width=width,
                      target_rank=target_rank, target_onto=target_onto,
    				  int_barcol=barcol)
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
    if (plttype == "pca"){
        prcomp_cli(args)
    }
    if (plttype == "ft_dist"){
    	ft_dist_cli(args)
    }
}

main()
