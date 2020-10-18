# options
options(stringsAsFactors = FALSE, message=FALSE, warnings=FALSE)

####### ==================== #######
#              LIBRARIES           #
####### ==================== #######
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
suppressMessages(library(jsonlite))
suppressMessages(library(stringr))

####### ==================== #######
#             CONSTANTS            #
####### ==================== #######
# mapping of integers to colors
grp_color_values <- c("dodgerblue", "darkorange",
    "yellow2", "red2", "darkviolet", "black")

X_AXIS_ROT <- 60

####### ==================== #######
#         UTILITY FUNCTIONS        #
####### ==================== #######
# read tab-separated file
read_result <- function(file){
    df <- read.delim(file, sep="\t", stringsAsFactors=FALSE)
    return(df)
}

# write tab-separated file of plot data
write_plot_table <- function(df, path){
    write.table(df, file=path, sep="\t", quote=FALSE, row.names=FALSE)
}

# for PCA and heatmaps, replace
# missing values with 1/1000 of minimum intensity in
# whole dataframe
impute <- function(mat) {
    mat[mat == 0] <- NA
    mat[is.na(mat)] <- min(mat, na.rm = TRUE) * 1e-3
    mat
}

# for pca plot: pad the plot area so that the labels
# display properly
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

check_ranks <- function(df, target_rank) {
    ranks_in_df <- unique(df[, "rank"])
    if (!(target_rank %in% ranks_in_df)) {
        stop('There are no annotations with the desired target_rank. Check data and try again.')
    }
}

short_onto_to_long <- c("bp" = "biological_process",
                        "mf"="molecular_function",
                        "cc"="cellular_component")

filter_to_desired_onto <- function(df, target_onto) {
    if (is.null(target_onto)) {
        stop('must provide target ontology for function or function-taxonomy mode',
             call. = FALSE)
    }
    long_onto <- short_onto_to_long[target_onto]
    df <- df[df[, "namespace"] == long_onto, ]
    # also, remove any bp, cc, or mf rows
    df <- df[!(df[, "name"] %in% short_onto_to_long), ]
    if (nrow(df) == 0) {
        stop(paste0("No terms in the dataframe come from the desired ontology (", target_onto, ")"))
    }
    return(df)
}

####### ==================== #######
#              BARPLOT             #
####### ==================== #######
mq_barplot <- function(df, img, mode, meancol,
                       nterms, width, height, target_rank, target_onto,
                       int_barcol, tabfile){
    if (!(meancol %in% names(df))){
        stop('Mean column name not found in dataframe. Check spelling and try again.',
             call. = FALSE)
    }
    # filter to desired rank, if taxonomy
    if (mode == "t"){
        check_ranks(df, target_rank)
        # filter
        df <- df[df[, "rank"] == target_rank, ]
    }
    if (mode == "f") {
        df <- filter_to_desired_onto(df, target_onto)
    }
    # exponentiate
    df[, meancol] <- 2^df[, meancol]

    # reorder, for taking the top N terms
    reord <- df[order(df[, meancol], decreasing=TRUE), ]
    # take top N terms or number of rows, whichever is less
    sub_reord <- reord[1:min(nterms, nrow(reord)), ]
    if (mode == "t"){
      barnamecol = "taxon_name"
      xlab = target_rank
    } else {
      barnamecol = "name"
      xlab = "Term"
    }
    # color mapping
    barcol <- grp_color_values[int_barcol]

    ggplot(sub_reord) +
      geom_bar(aes_(x = reorder(sub_reord[, barnamecol], -sub_reord[, meancol]),
                    y = as.name(meancol)), stat = "identity", fill = barcol, col = "black", position = "dodge") +
      theme_bw() +
      labs(x = xlab, y = "Total Peptide Abundance") +
      theme(axis.text.x = element_text(angle = X_AXIS_ROT, hjust = 1))
    ggsave(img, height = height, width = width, units = "in")

    # write tabfile
    if (!is.null(tabfile)) {
      write_plot_table(sub_reord, tabfile)
    }
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
    # 11. tabfile - path to tabular file to write plot data
    img <- args[2]
    infile <- args[3]
    df <- read_result(infile)
    mode <- args[4]
    meancol <- args[5]
    nterms <- as.numeric(args[6])
    width <- as.numeric(args[7])
    height <- as.numeric(args[8])
    target_rank <- args[9]
    if (target_rank == "None") target_rank <- NULL
    target_onto <- args[10]
    if (target_onto == "None") target_onto <- NULL
    barcol <- as.numeric(args[11])
    tabfile <- args[12]
    if (tabfile == "None") tabfile <- NULL
    plt <- mq_barplot(df, img=img, mode=mode,
                      meancol=meancol,
                      nterms=nterms,
                      height=height, width=width,
                      target_rank=target_rank, target_onto = target_onto,
                      int_barcol=barcol,
                      tabfile=tabfile)
}

####### ==================== #######
#              HEATMAP             #
####### ==================== #######
# correlation distance (not absolute)
cor.dist <- function(x){
    as.dist(1-cor(t(x)))
}

# ward method for hierarchical clustering
hclust.ward <- function(x) {
    hclust(x,method="ward.D")
}


library(scico)
heatmap_colors <- scico(30, palette = 'vik')

mq_heatmap <- function(img, df, all_intcols, colSideColors, filter_to_sig, alpha, width, height, strip, feature_cluster_size, sample_cluster_size, fc_corr_p, infilename){
    # df is the output from either expand, stat, or filter
    # samp_columns is a vector of all columns with the term intensities
    # colSide colors is a vector of colors for the groups. Must be in the same order as samp_columns
    # filter to sig
    if (filter_to_sig){
        if (fc_corr_p == "None"){
            stop("corrected p-value column not defined. did you run metaquantome stat?",
                 call. = FALSE)
        }
        df$corrected_p <- df[, fc_corr_p]
        pvals <- df$corrected_p
        if (is.null(pvals)) {
            stop("the dataset does not have a column named 'corrected_p'. did you run metaquantome stat?",
                 call. = FALSE)
        }
        df <- df[df$corrected_p < alpha, ]
        if (nrow(df) == 0) {
            stop("after filtering to provided value of alpha, no rows remain. please increase alpha",
                 call. = FALSE)
        }
    }

    # impute
    mat <- impute(data.matrix(df[, all_intcols]))

    # scale rows
    datmat.scale <- t(apply(mat, 1, scale))
    rownames(datmat.scale) <- df$id

    # replace strip with "" in column names
    if (strip != "None"){
    	colnames(datmat.scale) <- str_replace(colnames(df[, all_intcols]), strip, "")
    } else {
    	colnames(datmat.scale) <- colnames(df[, all_intcols])
    }

    # build dendrograms
    feature.dend <- as.dendrogram(hclust.ward(cor.dist(datmat.scale)))
    sample.dend <- as.dendrogram(hclust.ward(cor.dist(t(datmat.scale))))
    
    # Output cluster file for features and samples
    feature_cluster = cutree(hclust.ward(cor.dist(datmat.scale)),k=feature_cluster_size);
    write.table(feature_cluster, file = paste("feature_cluster_",infilename,'.txt', sep=""), sep = "\t", col.names=FALSE, quote=FALSE)
    sample_cluster = cutree(hclust.ward(cor.dist(t(datmat.scale))),k=sample_cluster_size);
    write.table(sample_cluster, file = paste("sample_cluster_",infilename,'.txt', sep=""), sep = "\t", col.names=FALSE, quote=FALSE)
    
    # write plot to img path
    png(filename=img, width=width, height=height, res=300, units="in")
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
              density.info = "none",
              labRow = "")
    # send message to the ether
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
    
    infilename = basename(infile)
    
    # split all_intcols from SampleGroups(), for samp_columns vector
    all_intcols <- get_all_intcols(args[4])

    # get color mapping
    colSideColors <- get_colors_from_groups(args[5], all_intcols)
    filter_to_sig <- (args[6] == "True")
    alpha <- as.numeric(args[7])
    width <- as.numeric(args[8])
    height <- as.numeric(args[9])
    strip <- args[10]
    feature_cluster_size = args[11]
    sample_cluster_size = args[12]
    fc_corr_p <- args[13]
    
    mq_heatmap(img, df, all_intcols, colSideColors, filter_to_sig, alpha, width, height, strip, feature_cluster_size, sample_cluster_size, fc_corr_p, infilename)
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
    within_variance <- sapply(1:nclust, function(i) mean(data.matrix((clust[[i]] - means[[i]])^2)))
    avg_dist / sum(within_variance)
}

mq_prcomp <- function(img, df, all_intcols, json_dump, colors, calculate_sep, width, height, strip, infilename){
    # function(img, df, all_intcols, colSideColors, filter_to_sig, alpha, width, height)
    # df is the result from filter
    # samp_columns is a vector of all intensity column names
    # cols is a vector of group colors
    mat <- impute(data.matrix(df[, all_intcols]))
    pr <- prcomp(mat, scale=TRUE, center=TRUE)
    
    
    # write PC1 and PC2 rotation data to file
    pc_data = pr$rotation[,1:2]
    pc_data = data.frame("Samples"=rownames(pc_data), "PC1"=pc_data[,1], "PC2"=pc_data[,2])
    write.table(pc_data, file = paste("PC_data_",infilename,'.txt', sep=""), sep = "\t", quote=FALSE, row.names=FALSE)
    
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

    png(filename=img, width=width, height=height, res=300, units="in")
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
    
    infilename = basename(infile)
    
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
        strip=strip, infilename)
}


####### ==================== #######
#              VOLCANO             #
####### ==================== #######
mq_volcano <- function(df, img, fc_name, fc_corr_p, flip_fc, width, height, textannot, gosplit, tabfile){
    # df is the dataframe after stat
    # fc_name is the name of the column with the fold change data
    # textcol is the name of the column with the text describing the term
    df$fc <- df[, fc_name]
    if (flip_fc){
        df$fc <- (-1)*df$fc
    }
    df$corrected_p <- df[, fc_corr_p]
    #df$neglog10p <- -log10(df[, "corrected_p"])
    df$neglog10p <- -log10(df$corrected_p)
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
        labs(x = "Log2 Fold Change", y = "-Log10 FDR-Corrected P Value")
    vplt
    ggsave(file=img, width=width, height=height, units="in", dpi=300)

    # write tabular data
    if (!is.null(tabfile)) {
        write_plot_table(df, tabfile)
    }
}

volcano_cli <- function(args){
    # args are as follows
    # 1. plot type (guaranteed to be 'volcano')
    # 2. pltfile - output image file
    # 3. input tabular file
    # 4. name of text annotation column
    # 5. name of fold change column
    # 6. name of corrected p-value column
    # 7. whether to flip fc
    # 8. whether to split GO by ontology/namespace
    # 9. image width (default 5)
    # 10. image height (default 5)
    
    img <- args[2]
    infile <- args[3]
    df <- read_result(infile)
    textannot <- args[4]
    fc_name <- args[5]
    fc_corr_p <- args[6]
    flip_fc <- (args[7] == "True")
    gosplit <- (args[8] == "True")
    width <- as.numeric(args[9])
    height <- as.numeric(args[10])
    tabfile <- args[11]
    if (tabfile == "None") tabfile <- NULL
    plt <- mq_volcano(df, img=img, textannot=textannot, fc_name=fc_name,
                      fc_corr_p=fc_corr_p, flip_fc=flip_fc, gosplit=gosplit, width=width,
                      height=height, tabfile=tabfile)
}

####### ==================== #######
#              FT DIST             #
####### ==================== #######
mq_ft_dist <- function(df, img, whichway, name, id, meancol,
                       nterms, width, height, target_rank, target_onto,
                       int_barcol, tabfile){
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
	    if (nrow(df) == 0) {
	        stop(paste0("The term ", name, " was not found in the dataframe."))
	    }
	}
	if (id != "None"){
		if (whichway == "t_dist"){
			df <- df[df[, "go_id"] == id, ]
		}
		if (whichway == "f_dist"){
			df <- df[df[, "tax_id"] == id, ]
		}
	    if (nrow(df) == 0) {
	        stop(paste0("The id ", id, " was not found in the dataframe."))
	    }
	}

    # filter to desired rank, if taxonomic distribution
    if (whichway == "t_dist"){
        check_ranks(df, target_rank)
        df <- df[df[, "rank"] == target_rank, ]
    } else if (whichway == "f_dist"){
        df <- filter_to_desired_onto(df, target_onto)
    } else {
    	stop("Wrong whichway - should be t_dist or f_dist.")
    }

    # calculate proportions
    df[, meancol] <- 2^df[, meancol]

    df$props <- df[, meancol] / sum(df[, meancol], na.rm=TRUE)
    reord <- df[order(df[, "props"], decreasing=TRUE), ]

    if (nterms == "all"){
        sub_reord <- reord
    } else {
        sub_reord <- reord[1:min(nterms, nrow(reord)), ]
    }

    if (whichway == "t_dist"){
        barnamecol = "taxon_name"
    } else {
        barnamecol = "name"
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

    ggplot(sub_reord) +
        geom_bar(aes_(x = reorder(sub_reord[, barnamecol], -sub_reord[, "props"]),
                      y = as.name("props")), stat = "identity", fill = barcol, col = "black") +
        theme_bw() +
        labs(x = xlab, y = "Proportion of Peptide Abundance") +
        theme(axis.text.x = element_text(angle = X_AXIS_ROT, hjust = 1))
    ggsave(img, height = height, width = width, units = "in")

    # write tabular data
    if (!is.null(tabfile)) {
        write_plot_table(sub_reord, tabfile)
    }
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
    tabfile <- args[14]
    if (tabfile == "None") tabfile <- NULL
    plt <- mq_ft_dist(df, img=img, whichway=whichway,
    				  name=name, id=id,
                      meancol=meancol,
                      nterms=nterms,
                      height=height, width=width,
                      target_rank=target_rank, target_onto=target_onto,
    				  int_barcol=barcol, tabfile=tabfile)
}

####### ==================== #######
#           STACKED BAR            #
####### ==================== #######

mq_stacked <- function(img, df, all_intcols, json_dump, nterms, target_rank, width, height, tabfile){
  # df is the dataframe after stat
  # nterms is the number of taxa to show
  
  grp_list <- fromJSON(json_dump)
  
  grp_df <- grp_list %>% 
    as.data.frame() %>% 
    pivot_longer(cols = 1:ncol(.), names_to = "samplegroup", values_to = "sample")
  
  # parse out sample groups, exponentiate, calculate relative abundance
  dat <- df %>% 
    pivot_longer(all_intcols, names_to = "sample", values_to = "abundance") %>% 
    full_join(grp_df) %>% 
    mutate(replicate = str_replace(string = sample, pattern = samplegroup, replacement = "")) %>% 
    replace_na(list(abundance = 0)) %>%
    mutate(abundance = 2^abundance) %>%
    select(sample, samplegroup, replicate, id, name, rank, abundance) %>% 
    filter(rank == target_rank) %>% 
    group_by(sample) %>% 
    mutate(abundance = 100*abundance/sum(abundance))
  
  # reorder taxa levels for plotting
  taxa_levels <- names(sort(tapply(dat$abundance, dat$name, sum)))
  
  # collapse less abundang terms into "Other" if terms exceed desired terms
  if(length(unique(dat$name))>nterms){
    topn <- tail(taxa_levels, nterms)
    dat<-dat %>%
      mutate(name = factor(ifelse(name %in% topn, name, "Other"))) %>% 
      group_by(sample, samplegroup, replicate, name) %>% 
      summarise(abundance = sum(abundance)) %>% 
      group_by(sample)
    taxa_levels <- names(sort(tapply(dat$abundance, dat$name, sum)))
    # reorder taxa levels for plotting, to have "Other" on top
    other_index <- as.numeric(which(taxa_levels == "Other"))
    taxa_levels <- c("Other", taxa_levels[1:(other_index-1)], taxa_levels[(other_index+1):length(taxa_levels)])
  }
  
  # make stacked bar plot
  fig <- dat %>%     
    ggplot(aes(x=replicate, y=abundance, fill=factor(name, levels = taxa_levels)))+
    geom_bar(position="stack", stat = "identity") +
    facet_grid(cols = vars(samplegroup)) +
    labs(x= "Sample", y="Relative Abundance")+
    scale_fill_brewer(name="Taxa", palette = "Set1")
  
  # write tabular file
  write.table(x = dat, file = tabfile, quote = FALSE, row.names = FALSE)

  # save plot
  ggsave(file=img, width=width, height=height, units="in", dpi=300)
  
}

stacked_cli <- function(args){
  
  img <- args[2]
  infile <- args[3]
  df <- read_result(infile)
  
  # split all_intcols from SampleGroups(), for samp_columns vector
  all_intcols <- unlist(strsplit(get_all_intcols(args[4]), split=","))
  
  # other args
  json_dump <- args[5]
  nterms <- as.numeric(args[6])
  target_rank <- toString(args[7])
  width <- as.numeric(args[8])
  height <- as.numeric(args[9])
  tabfile <- args[10]
  
  mq_stacked(img, df, all_intcols, json_dump, nterms, target_rank, width, height, tabfile)
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
