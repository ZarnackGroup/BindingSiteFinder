
#' Diagnostic plots for the differential binding background
#'
#' To perform differential binding analysis between two conditions the
#' \code{\link{calculateBsBackground}} function groups crosslinks per gene
#' into those from binding sites and those from background regions.
#' The \code{\link{filterBsBackground}} function perfroms certain
#' filtering operations on that background to ensure that it's suitable
#' for differential testing. This function visually displays the effect
#' of these filtering operations.
#'
#'
#' @param object a \code{\link{BSFDataSet}} object with background counts filtered
#' by \code{\link{filterBsBackground}}
#' @param filter character; which filter to display in the plot (one of:
#' 'minCounts', 'balanceBackground', 'balanceCondition')
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{calculateBsBackground}}
#' @seealso \code{\link{filterBsBackground}}
#'
#' @import ggplot2
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#'
#' # make testset
#' bds = makeBindingSites(bds, bsSize = 7)
#' bds = assignToGenes(bds, anno.genes = gns)
#' bds = imputeBsDifferencesForTestdata(bds)
#' bds = calculateBsBackground(bds, anno.genes = gns, use.offset = FALSE)
#'
#' # use all filters and remove binding sites that fail (default settings)
#' bds = filterBsBackground(bds)
#'
#' # display minCount filter
#' plotBsBackgroundFilter(bds, filter = "minCounts")
#'
#' # display balance background filter
#' plotBsBackgroundFilter(bds, filter = "balanceBackground")
#'
#' # display balance condition filter
#' plotBsBackgroundFilter(bds, filter = "balanceCondition")
#'
#' @export
plotBsBackgroundFilter <- function(object, filter = c("minCounts", "balanceBackground" ,"balanceCondition")) {

    # bind locally used variables
    s <- r <- name <- value <- ratio.ref <- ratio.comp <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    if (is.null(object@params$filterBsBackground)) {
        msg0 = paste0("Global filter was not applied yet. Run filterBsBackground() to compute background first. \n")
        stop(msg0)
    }

    # handle display plot options
    filter = match.arg(filter, choices = c("minCounts", "balanceBackground", "balanceCondition"))

    # check for plot specific data
    if (is.null(object@plotData$filterBsBackground$data.minCounts) & filter == "minCounts") {
        msg0 = paste0("Plot selection: '", filter, "' is not possible.")
        msg1 = paste0("Make sure to run filterBsBackground with option '", filter, "' first.")
        stop(c(msg0,msg1))
    }
    if (is.null(object@plotData$filterBsBackground$data.balanceBackground) & filter == "balanceBackground") {
        msg0 = paste0("Plot selection: '", filter, "' is not possible.")
        msg1 = paste0("Make sure to run filterBsBackground with option '", filter, "' first.")
        stop(c(msg0,msg1))
    }
    if (is.null(object@plotData$filterBsBackground$data.balanceCondition) & filter == "balanceCondition") {
        msg0 = paste0("Plot selection: '", filter, "' is not possible.")
        msg1 = paste0("Make sure to run filterBsBackground with option '", filter, "' first.")
        stop(c(msg0,msg1))
    }


    # MAIN
    # --------------------------------------------------------------------------

    # make plot for minCounts filter
    if (filter == "minCounts") {
        # get options used
        optstr = object@params$filterBsBackground
        optstrNice = paste0("Cutoff=", optstr$minCounts.cutoff)

        # get data
        df = object@plotData$filterBsBackground$data.minCounts %>%
            mutate(r = rank(s))

        p = ggplot(df, aes(x = r, y = s)) +
            ggrastr::rasterise(ggpointdensity::geom_pointdensity(size = 2), dpi = 300) +
            theme_bw() +
            geom_hline(yintercept = optstr$minCounts.cutoff, linetype = "dashed") +
            scale_y_log10() +
            viridis::scale_color_viridis(option = "rocket") +
            theme(legend.position = "top") +
            guides(color = guide_colorbar(title.position = 'top', title.hjust = 0.5,
                                          barwidth = unit(20, 'lines'), barheight = unit(.5, 'lines'))) +
            annotate("text", x = max(df$r)/2, y = optstr$minCounts.cutoff,
                     label = paste0("Cutoff=", optstr$minCounts.cutoff),
                     vjust = -1.1) +
            labs(
                title = "filterBsBackground() - minCounts",
                subtitle = optstrNice,
                x = "Rank",
                y = "#Crosslinks per gene",
                color = "#Genes"
            )
    }

    # make plot for balanceBackground filter
    if (filter == "balanceBackground") {
        # get options used
        optstr = object@params$filterBsBackground
        optstrNice = paste0("Cutoff.bs=", optstr$balanceBackground.cutoff.bs,
                            " cutoff.bg=", optstr$balanceBackground.cutoff.bg)

        # get data
        df = object@plotData$filterBsBackground$data.balanceBackground %>%
            mutate(name = ifelse(name == "ratio.bg", "Background", "Binding site"))

        df = df[!is.na(df$value),]

        p = ggplot(df, aes(x = name, y = value)) +
            ggdist::stat_halfeye(aes(fill = name), adjust = 0.5, width = .6, .width = 0, justification = -.2, point_colour = NA) +
            geom_boxplot(aes(color = name), width = .12, outlier.size = 0.5) +
            geom_point(data = subset(df, name == "ratio.bg" & value < optstr$balanceBackground.cutoff.bg), aes(fill = name), shape = 21, alpha = 0.3) +
            geom_point(data = subset(df, name == "ratio.bs" & value > optstr$balanceBackground.cutoff.bs), aes(fill = name), shape = 21, alpha = 0.3) +
            annotate("rect", xmin = 1.75, xmax = 2.25, ymin = optstr$balanceBackground.cutoff.bs, ymax = 1, alpha = .3, fill = "#435B66") +
            annotate("rect", xmin = 0.75, xmax = 1.25, ymin = 0, ymax = optstr$balanceBackground.cutoff.bg, alpha = .3, fill = "#A76F6F") +
            annotate("text", x = 2.25, y = 0.5, vjust = -1, label = paste0("Remove genes with bs ratio above: ", optstr$balanceBackground.cutoff.bs)) +
            annotate("text", x = 1.25, y = 0.5, vjust = -1, label = paste0("Remove genes with bg ratio below: ", optstr$balanceBackground.cutoff.bg)) +
            coord_flip() +
            theme_bw() +
            scale_fill_manual(values = c("#A76F6F", "#435B66")) +
            scale_color_manual(values = c("#A76F6F", "#435B66")) +
            theme(legend.position = "none") +
            labs(
                title = "filterBsBackground() - balanceBackground",
                subtitle = optstrNice,
                x = NULL,
                y = "Background to binding site count ratio"
            )
    }

    # make plot for balanceCondition filter
    if (filter == "balanceCondition") {
        # get options used
        optstr = object@params$filterBsBackground
        optstrNice = paste0("Cutoff=", optstr$balanceCondition.cutoff)

        # get reference condition
        this.meta = getMeta(object)
        this.condition.reference = levels(this.meta$condition)[1]
        this.condition.comp = levels(this.meta$condition)[2]

        # get data
        df = object@plotData$filterBsBackground$data.balanceCondition %>% as.data.frame()
        df$ratio.ref = df[,3] / df[,2]
        df$ratio.comp = df[,4] / df[,2]
        df = df %>% select(ratio.ref, ratio.comp) %>% pivot_longer(everything())

        p = ggplot(df, aes(x = name, y = value)) +
            ggdist::stat_halfeye(aes(fill = name), adjust = 0.5, width = .6, .width = 0, justification = -.2, point_colour = NA) +
            geom_boxplot(aes(color = name), width = .12, outlier.size = 0.5) +
            geom_point(data = subset(df, name == "ratio.ref" & value < optstr$balanceCondition.cutoff), aes(fill = name), shape = 21, alpha = 0.3) +
            geom_point(data = subset(df, name == "ratio.comp" & value < optstr$balanceCondition.cutoff), aes(fill = name), shape = 21, alpha = 0.3) +
            annotate("rect", xmin = 0.05, xmax = 2.95, ymin = 0, ymax = optstr$balanceCondition.cutoff, alpha = .3, fill = "darkgrey") +
            # annotate("rect", xmin = 1.75, xmax = 2.25, ymin = 0, ymax = 0.02, alpha = .3, fill = "#435B66") +
            # annotate("rect", xmin = 0.75, xmax = 1.25, ymin = 0.98, ymax = 1, alpha = .3, fill = "#A76F6F") +
            coord_flip() +
            theme_bw() +
            scale_fill_manual(values = c("#A76F6F", "#435B66")) +
            scale_color_manual(values = c("#A76F6F", "#435B66")) +
            theme(legend.position = "none") +
            labs(
                title = "filterBsBackground() - balanceCondition",
                subtitle = optstrNice,
                x = NULL,
                y = "Condition count ratio"
            ) +
            scale_x_discrete(labels = c(paste0(this.condition.comp, "\n (comparison)"),
                                        paste0(this.condition.reference, "\n (reference)")))


    }

    return(p)
}


#' MA style plot
#'
#' Wrapper that plots differential binding results as MA plot. For each binding
#' site the estimated baseMean (log2) is shown on X and the fold-change (log2)
#' is shown on Y.
#'
#' @param object a \code{\link{BSFDataSet}} object with results calculated by
#' \code{\link{calculateBsFoldChange}}
#' @param what character; whether to show results for binding sites or the
#' background (one of: 'bs', 'bg')
#' @param sig.threshold numeric; what P value significance level to use
#' (default = 0.05)
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{calculateBsFoldChange}}
#'
#' @import ggplot2
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#'
#' # make testset
#' bds = makeBindingSites(bds, bsSize = 7)
#' bds = assignToGenes(bds, anno.genes = gns)
#' bds = imputeBsDifferencesForTestdata(bds)
#' bds = calculateBsBackground(bds, anno.genes = gns, use.offset = FALSE)
#'
#' # use all filters and remove binding sites that fail (default settings)
#' bds = filterBsBackground(bds)
#'
#' # calculate fold-changes
#' bds = calculateBsFoldChange(bds)
#'
#' # make MA plot
#' plotBsMA(bds)
#'
#' @export
plotBsMA <- function(object,
                     what = c("bs", "bg"),
                     sig.threshold = 0.05) {

    # bind locally used variables
    bs.padj <- bs.log2FoldChange <- sig <- bs.baseMean <- bg.padj <- NULL
    bg.log2FoldChange <- bg.baseMean <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    if (is.null(object@params$calculateBsFoldChange)) {
        msg0 = paste0("Fold-changes were not calculated yet. Run calculateBsFoldChange() first. \n")
        stop(msg0)
    }

    # handle display plot options
    what = match.arg(what, choices = c("bs", "bg"))

    this.ranges = getRanges(object)
    expected.cols = c("bs.padj", "bs.log2FoldChange", "bs.baseMean",
                      "bg.padj", "bg.log2FoldChange", "bg.baseMean")
    if (!all(expected.cols %in% colnames(as.data.frame(mcols(this.ranges))))) {
        msg0 = paste0("MA plot not possible, results columns: ", paste(expected.cols, collapse = ','), " not present.\n")
        msg1 = paste0("Make sure to run calculateBsFoldChange() first.\n")
        stop(c(msg0, msg1))
    }

    # MAIN
    # --------------------------------------------------------------------------

    optstr = object@params$calculateBsFoldChange
    optstrNice = paste0("alpha=", optstr$alpha)

    bright_up_down_not = c("#999999", "#68b1a5", "#874C62")
    dark_up_down_not = c("#4d4d4d", "#2b544d", "#623747")

    if (what == "bs") {
        df = as.data.frame(mcols(this.ranges)) %>%
            mutate(sig = ifelse(bs.padj < sig.threshold & bs.log2FoldChange > 0, "Up",
                                ifelse(bs.padj < sig.threshold & bs.log2FoldChange < 0, "Down", "Not"))) %>%
            mutate(sig = factor(sig, levels = c("Not", "Up", "Down"))) %>%
            arrange(sig)

        p = ggplot(df, aes(x = log2(bs.baseMean), y = bs.log2FoldChange, color = sig, fill = sig)) +
            ggrastr::rasterise(geom_point(shape = 21, stroke = 0.5, size = 1.5), dpi = 300) +
            geom_hline(yintercept = 0, color = "black", alpha = .5) +
            theme_bw() +
            scale_fill_manual(values = bright_up_down_not) +
            scale_color_manual(values = dark_up_down_not) +
            guides(color = guide_legend(override.aes = list(size = 4))) +
            theme(legend.key.size = unit(1, 'cm'), legend.position = "top") +
            labs(
                title = "plotBsMA() - binding sites",
                subtitle = optstrNice,
                x = "Base mean",
                y = "Fold-change (log2)",
                color = "Regulation",
                fill = "Regulation")
    }
    if (what == "bg") {
        df = as.data.frame(mcols(this.ranges)) %>%
            mutate(sig = ifelse(bg.padj < sig.threshold & bg.log2FoldChange > 0, "Up",
                                ifelse(bg.padj < sig.threshold & bg.log2FoldChange < 0, "Down", "Not"))) %>%
            mutate(sig = factor(sig, levels = c("Not", "Up", "Down"))) %>%
            arrange(sig)

        p = ggplot(df, aes(x = log2(bg.baseMean), y = bg.log2FoldChange, color = sig, fill = sig)) +
            ggrastr::rasterise(geom_point(shape = 21, stroke = 0.5, size = 1.5), dpi = 300) +
            geom_hline(yintercept = 0, color = "black", alpha = .5) +
            theme_bw() +
            scale_fill_manual(values = bright_up_down_not) +
            scale_color_manual(values = dark_up_down_not) +
            guides(color = guide_legend(override.aes = list(size = 4))) +
            theme(legend.key.size = unit(1, 'cm'), legend.position = "top") +
            labs(
                title = "plotBsMA() - background",
                subtitle = optstrNice,
                x = "Base mean",
                y = "Fold-change (log2)",
                color = "Regulation",
                fill = "Regulation")
    }


    return(p)
}

#' Volcano style plot
#'
#' Wrapper that plots differential binding results as volcano plot. For each binding
#' site the estimated fold-change (log2) is shown on X and the adjusted
#' P value (-log10) is shown on Y.
#'
#' @param object a \code{\link{BSFDataSet}} object with results calculated by
#' \code{\link{calculateBsFoldChange}}
#' @param what character; whether to show results for binding sites or the
#' background (one of: 'bs', 'bg')
#' @param sig.threshold numeric; what P value significance level to use
#' (default = 0.05)
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{calculateBsFoldChange}}
#'
#' @import ggplot2
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#'
#' # make testset
#' bds = makeBindingSites(bds, bsSize = 7)
#' bds = assignToGenes(bds, anno.genes = gns)
#' bds = imputeBsDifferencesForTestdata(bds)
#' bds = calculateBsBackground(bds, anno.genes = gns, use.offset = FALSE)
#'
#' # use all filters and remove binding sites that fail (default settings)
#' bds = filterBsBackground(bds)
#'
#' # calculate fold-changes
#' bds = calculateBsFoldChange(bds)
#'
#' # make volcano plot
#' plotBsVolcano(bds)
#'
#' @export
plotBsVolcano <- function(object,
                          what = c("bs", "bg"),
                          sig.threshold = 0.05) {

    # bind locally used variables
    bs.padj <- bs.log2FoldChange <- sig <- bg.padj <- bg.log2FoldChange <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    if (is.null(object@params$calculateBsFoldChange)) {
        msg0 = paste0("Fold-changes were not calculated yet. Run calculateBsFoldChange() first. \n")
        stop(msg0)
    }

    # handle display plot options
    what = match.arg(what, choices = c("bs", "bg"))

    this.ranges = getRanges(object)
    expected.cols = c("bs.padj", "bs.log2FoldChange", "bs.baseMean",
                      "bg.padj", "bg.log2FoldChange", "bg.baseMean")
    if (!all(expected.cols %in% colnames(as.data.frame(mcols(this.ranges))))) {
        msg0 = paste0("MA plot not possible, results columns: ", paste(expected.cols, collapse = ','), " not present.\n")
        msg1 = paste0("Make sure to run calculateBsFoldChange() first.\n")
        stop(c(msg0, msg1))
    }

    # MAIN
    # --------------------------------------------------------------------------

    optstr = object@params$calculateBsFoldChange
    optstrNice = paste0("alpha=", optstr$alpha)

    bright_up_down_not = c("#999999", "#68b1a5", "#874C62")
    dark_up_down_not = c("#4d4d4d", "#2b544d", "#623747")

    if (what == "bs") {
        df = as.data.frame(mcols(this.ranges)) %>%
            mutate(sig = ifelse(bs.padj < sig.threshold & bs.log2FoldChange > 0, "Up",
                                ifelse(bs.padj < sig.threshold & bs.log2FoldChange < 0, "Down", "Not"))) %>%
            mutate(sig = factor(sig, levels = c("Not", "Up", "Down"))) %>%
            arrange(sig)

        p = ggplot(df, aes(x = bs.log2FoldChange, y = -log10(bs.padj), color = sig, fill = sig)) +
            ggrastr::rasterise(geom_point(shape = 21, stroke = 0.5, size = 1.5), dpi = 300) +
            geom_vline(xintercept = 0, color = "black", alpha = .5) +
            theme_bw() +
            scale_fill_manual(values = bright_up_down_not) +
            scale_color_manual(values = dark_up_down_not) +
            guides(color = guide_legend(override.aes = list(size = 4))) +
            theme(legend.key.size = unit(1, 'cm'), legend.position = "top") +
            labs(
                title = "plotBsVolcano() - binding sites",
                subtitle = optstrNice,
                x = "Fold-change (log2)",
                y = "Adjusted P value (-log10)",
                color = "Regulation",
                fill = "Regulation")
    }
    if (what == "bg") {
        df = as.data.frame(mcols(this.ranges)) %>%
            mutate(sig = ifelse(bg.padj < sig.threshold & bg.log2FoldChange > 0, "Up",
                                ifelse(bg.padj < sig.threshold & bg.log2FoldChange < 0, "Down", "Not"))) %>%
            mutate(sig = factor(sig, levels = c("Not", "Up", "Down"))) %>%
            arrange(sig)

        p = ggplot(df, aes(x = bg.log2FoldChange, y = -log10(bg.padj), color = sig, fill = sig)) +
            ggrastr::rasterise(geom_point(shape = 21, stroke = 0.5, size = 1.5), dpi = 300) +
            geom_vline(xintercept = 0, color = "black", alpha = .5) +
            theme_bw() +
            scale_fill_manual(values = bright_up_down_not) +
            scale_color_manual(values = dark_up_down_not) +
            guides(color = guide_legend(override.aes = list(size = 4))) +
            theme(legend.key.size = unit(1, 'cm'), legend.position = "top") +
            labs(
                title = "plotBsVolcano() - background",
                subtitle = optstrNice,
                x = "Fold-change (log2)",
                y = "Adjusted P value (-log10)",
                color = "Regulation",
                fill = "Regulation")
    }

    return(p)
}


