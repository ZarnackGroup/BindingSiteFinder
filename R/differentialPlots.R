
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


#' Gene Regulation Plot
#'
#' Display the fold-change of all binding sites from a given gene on a relative
#' per-nucleotide scale. Binding sites are displayed as dots and with increasing
#' log2 fold-change, they deviate stronger from the center line.
#'
#' For this function to work, binding sites must be assigned to hosting genes
#' using \code{\link{assignToGenes}}. It is also recommended to assing
#' binding sites to transcript regions with \code{\link{assignToTranscriptRegions}}.
#'
#' It is also necessary to calculate the log2 fold-change of binding sites between
#' two conditions using the differential binding workflow \code{\link{calculateBsFoldChange}}.
#'
#' If in addition the transcript regions of the binding sites are given, then
#' shapes are changed accordingly. An edge case can arise from the merging of two
#' \code{\link{BSFDataSet}} objects. If binding sites are overlapping and
#' slightly offset close to the end of a particular transcript region annotation,
#' they might be assigned to different regions in both objects. This results in
#' some ambiguity after the merge, where for instance a binding site can be
#' assigned to CDS and 3'UTR. To handle how such edge cases are displayed, the
#' \code{transcript.regions.outlier.handle} exists. As default, simply the
#' region of the object that was merged first is shown. If one is interested in
#' showing all regions, then the options \code{both} displays both annotations
#' at the same time and labels them accordingly.
#'
#' @param object object; a \code{\link{BSFDataSet}} object
#' @param plot.geneID character;  the gene id of the gene to display. The id must
#' match with the gene ids given in the annotation object.
#' @param anno.annoDB object; an object of class \code{OrganismDbi} that contains
#' the gene annotation (!!! Experimental !!!).
#' @param anno.genes object; an object of class \code{\link{GenomicRanges}}
#' that represents the gene ranges directly.
#' @param match.geneID character; meta column name of the gene ID
#' @param match.geneName character; meta column name of the gene name
#' @param plot.gene.n.tiles numeric; number of tiles the gene should be split in
#' @param alpha numeric; the alpha value to show significantly regulated
#' binding sites. This should match the alpha value used in \code{\link{calculateBsFoldChange}}.
#' @param lfc.cutoff numeric; log2 fold-change cutoff to show significantly
#' regulated binding sites. This should match the lfc.cutoff value used in
#' \code{\link{calculateBsFoldChange}}.
#' @param transcript.regions.outlier.handle character; the option how to handle
#' multiple transcript region annotations being present for the same binding
#' site.
#' @param quiet logical; whether to print messages
#'
#' @return an object of class \code{ggplot2}
#'
#' @seealso \code{\link{BSFind}},
#' \code{\link{calculateBsFoldChange}}
#' \code{\link{assignToGenes}}
#' \code{\link{assignToTranscriptRegions}}
#'
#' @import ggplot2
#'
#' @examples
#'
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[2])
#'
#' # make testset
#' bds = makeBindingSites(bds, bsSize = 7)
#' bds = assignToGenes(bds, anno.genes = gns)
#' bds = assignToTranscriptRegions(object = bds, anno.transcriptRegionList = regions)
#' bds = imputeBsDifferencesForTestdata(bds)
#' bds = calculateBsBackground(bds, anno.genes = gns, use.offset = FALSE)
#'
#' # use all filters and remove binding sites that fail (default settings)
#' bds = filterBsBackground(bds)
#'
#' # calculate fold-changes
#' bds = calculateBsFoldChange(bds)
#'
#' # make example plot
#' exampleGeneId = "ENSG00000253352.10"
#' geneRegulationPlot(bds, plot.geneID = exampleGeneId, anno.genes = gns)
#'
#' @export
geneRegulationPlot <- function(object,
                               plot.geneID = NULL,
                               anno.annoDB = NULL,
                               anno.genes = NULL,
                               match.geneID = "gene_id",
                               match.geneName = "gene_name",
                               plot.gene.n.tiles = 100,
                               alpha = 0.05,
                               lfc.cutoff = 2,
                               transcript.regions.outlier.handle = c("first", "second", "both", "remove"),
                               quiet = FALSE) {

    # bind locally used variables
    bs.padj <- bs.log2FoldChange <- sig <- sigDir <- plot.position <- transcriptRegion <-  NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    # type checks
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is.logical(quiet))

    # check presence of fold-changes -> must be present
    if (is.null(object@params$calculateBsFoldChange)) {
        msg0 = paste0("Fold-changes were not calculated yet. Run calculateBsFoldChange() first. \n")
        stop(msg0)
    }
    # check presence of gene assignment -> must be present
    if (is.null(getRanges(object)$geneID)) {
        msg0 = paste0("Gene assignment was not calculated. Run assignToGenes() first. \n")
        stop(msg0)
        # TODO  ideally this should test for the param present in object@params$assignToGenes,
        #       but this gets overwritten by the mergeing of two bds objects. Change this in the future.
    }
    # check presence of transcript regions -> can be missing
    if (is.null(getRanges(object)$transcriptRegion)) {
        msg0 = paste0("Transcript region assignment was not calculated. Run assignToTranscriptRegions() first to link binding sites to transcript regions. \n")
        warning(msg0)
        # TODO  ideally this should test for the param present in object@params$assignToTranscriptRegions,
        #       but this gets overwritten by the mergeing of two bds objects. Change this in the future.
        # TODO make a version of the plot that does not encode the transcript regions as shape information
    }

    # check presence of correct annotation resource -> one of both must be present
    # Check if none is specified
    if (is.null(anno.annoDB) & is.null(anno.genes)) {
        msg = paste0("None of the required annotation sources anno.annoDB or anno.genes was specified. ")
        stop(msg)
    }
    # Check if both are specified
    if (!is.null(anno.annoDB) & !is.null(anno.genes)) {
        msg = paste0("Both of the required annotation sources anno.annoDB or anno.genes are specified. Please provide only one of the two. ")
        stop(msg)
    }
    # Checks if anno.annoDB should be used
    if (!is.null(anno.annoDB) & is.null(anno.genes)) {
        stopifnot(is(anno.annoDB, "OrganismDb"))
        if (!is.null(anno.genes)) {
            msg = paste0("Parameter anno.annoDB and anno.genes are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            datasource = "anno.annoDB"
            # extract relevant annotation
            anno.genes = genes(anno.annoDB, columns=c("ENSEMBL", "SYMBOL", "GENEID"))
            # Create matching vectors for columns from input annotation
            # --------------------------------------------------------------------------
            selectID = as.character(anno.genes$GENEID)
            selectName = as.character(anno.genes$SYMBOL)
        }
    }
    # Checks if anno.genes should be used
    if (is.null(anno.annoDB) & !is.null(anno.genes)) {
        stopifnot(is(anno.genes, "GenomicRanges"))
        if (!is.null(anno.annoDB)) {
            msg = paste0("Parameter anno.annoDB and anno.genes are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            datasource = "anno.genes"
            # extract relevant annotation
            # check for duplicated annotation columns
            if (sum(colnames(mcols(anno.genes)) == match.geneID) > 1) {
                msg = paste0("The names of multiple columns of the annotation match the match.geneID parameter. Please use a unique column name for match.geneID in anno.genes.")
                stop(msg)
            }
            if (sum(colnames(mcols(anno.genes)) == match.geneName) > 1) {
                msg = paste0("The names of multiple columns of the annotation match the match.geneName parameter. Please use a unique column name for match.geneName in anno.genes.")
                stop(msg)
            }
            # check correct annotation columns
            inNames = c(match.geneID, match.geneName)
            annoColNames = colnames(mcols(anno.genes))
            presentNames = inNames[inNames %in% annoColNames]

            # extract columns from annotation based on available names
            present.cols = lapply(presentNames, function(x){
                match.col = mcols(anno.genes)[match(x, colnames(mcols(anno.genes)))][[1]]
                return(match.col)
            })
            names(present.cols) = presentNames

            # extract values from annotation for all present columns
            if (match.geneID %in% names(present.cols)) {
                selectID = present.cols[[match(match.geneID, names(present.cols))]]
            } else {
                msg = paste0("No meta column for ", match.geneID, " present. Creating custom ID.\n")
                if (!quiet) message(msg)
                selectID = paste0("CUSTOM", seq_along(anno.genes))
            }
            if (match.geneName %in% names(present.cols)) {
                selectName = present.cols[[match(match.geneName, names(present.cols))]]
            } else {
                msg = paste0("No meta column for ", match.geneName, " present.\n")
                if (!quiet) message(msg)
            }
        }
    }

    # check if gene to plot is in bds object -> must be present
    if (!(plot.geneID %in% getRanges(object)$geneID)) {
        msg0 = stop(paste0("Selected geneID: ", plot.geneID, " is not present in the meta columns of the given bds object. \n"))
        stop(msg0)
    }
    # check if gene to plot is in annotation object -> must be present
    if (!(plot.geneID %in% selectID)) {
        msg0 = stop(paste0("Selected geneID: ", plot.geneID, " is not present in the meta columns of the given gene annotation. \n"))
        stop(msg0)
    }

    transcript.regions.outlier.handle = match.arg(transcript.regions.outlier.handle, choices = c("first", "second", "both", "remove"))

    # --------------------------------------------------------------------------
    # MAIN
    # --------------------------------------------------------------------------

    # get ranges of binding sites for selected gene
    plot.bs = getRanges(object)
    plot.bs = plot.bs[plot.bs$geneID == plot.geneID]
    plot.bs = sort(plot.bs)
    plot.bs = resize(plot.bs, width = 1, fix = "center")

    # get gene ranges for selected gene
    plot.gene = anno.genes[anno.genes$gene_id == plot.geneID]
    plot.gene.name = plot.gene$gene_name

    # claculate relative position of binding sites in gene
    plot.gene.tiles = unlist(tile(plot.gene, n = plot.gene.n.tiles))
    plot.gene.tiles$position = 1:length(plot.gene.tiles)

    ols = findOverlaps(plot.bs, plot.gene.tiles)

    gene.strand = unique(as.character(strand(plot.gene)))
    if(gene.strand == "-") {
        position = plot.gene.n.tiles-(subjectHits(ols))
    } else {
        position = subjectHits(ols)
    }
    plot.bs$plot.position = position

    # --------------------------------------------------------------------------
    # PLOTS
    # --------------------------------------------------------------------------

    # make plots without transcript region annotation
    if (is.null(getRanges(object)$transcriptRegion)) {
        df.plot = mcols(plot.bs) %>% as.data.frame() %>%
            mutate(sig = (ifelse(bs.padj < alpha & abs(bs.log2FoldChange) > log2(lfc.cutoff), TRUE, FALSE))) %>%
            mutate(dir = factor(ifelse(bs.log2FoldChange > 0, "Up", "Down"), levels = c("Up", "Down"))) %>%
            mutate(sigDir = ifelse(sig == TRUE & dir == "Up", "Up", ifelse(sig == TRUE & dir == "Down", "Down", "Not"))) %>%
            mutate(sigDir = factor(sigDir, levels = c("Not", "Up", "Down")))

        p = ggplot(df.plot, aes(x = plot.position, y = bs.log2FoldChange, color = sigDir, fill = sigDir)) +
            geom_hline(yintercept = 0, color = "#4d4d4d") +
            geom_segment(aes(x=plot.position, xend=plot.position, y=0, yend=bs.log2FoldChange), size = 1, color = "#4d4d4d") +
            geom_point(size = 4, stroke = 1, shape = 21) +
            scale_fill_manual(values = c("Not" = "#999999", "Up" = "#68b1a5", "Down" = "#874C62")) +
            scale_color_manual(values = c("Not" = "#4d4d4d", "Up" = "#2b544d", "Down" = "#623747")) +
            theme_minimal() +
            guides(color = guide_legend(override.aes = list(size = 4))) +
            theme(legend.key.size = unit(1, 'cm'), legend.position = "top") +
            ylim(-max(abs(df.plot$bs.log2FoldChange)), max(abs(df.plot$bs.log2FoldChange))) +
            xlim(0,plot.gene.n.tiles) +
            labs(
                title = paste0(plot.gene.name, " (", plot.geneID, ")"),
                x = "Relative position on the gene",
                y = "Fold-change (log2)",
                fill = "Regulation",
                color = "Regulation"
            )
        return(p)
    }

    # make plots with transcript region annotation
    if (!is.null(getRanges(object)$transcriptRegion)) {

        expected.regions = c("INTRON", "CDS", "UTR3", "UTR5", "OTHER")
        all.regions.present = unique(plot.bs$transcriptRegion)

        # plot cases where no region duplication is present
        if (all(all.regions.present %in% expected.regions)) {
            df.plot = mcols(plot.bs) %>% as.data.frame() %>%
                mutate(sig = (ifelse(bs.padj < alpha & abs(bs.log2FoldChange) > log2(lfc.cutoff), TRUE, FALSE))) %>%
                mutate(dir = factor(ifelse(bs.log2FoldChange > 0, "Up", "Down"), levels = c("Up", "Down"))) %>%
                mutate(sigDir = ifelse(sig == TRUE & dir == "Up", "Up", ifelse(sig == TRUE & dir == "Down", "Down", "Not"))) %>%
                mutate(sigDir = factor(sigDir, levels = c("Not", "Up", "Down"))) %>%
                mutate(duplicated = "no")

            p = ggplot(df.plot, aes(x = plot.position, y = bs.log2FoldChange, color = sigDir, fill = sigDir, shape = transcriptRegion)) +
                geom_hline(yintercept = 0, color = "#4d4d4d") +
                geom_segment(aes(x=plot.position, xend=plot.position, y=0, yend=bs.log2FoldChange), size = 1, color = "#4d4d4d") +
                geom_point(size = 4, stroke = 1) +
                scale_shape_manual(values=c("INTRON" = 23, "CDS" = 22, "UTR3" = 25, "UTR5" = 24, "OTHER" = 21)) +
                scale_fill_manual(values = c("Not" = "#999999", "Up" = "#68b1a5", "Down" = "#874C62")) +
                scale_color_manual(values = c("Not" = "#4d4d4d", "Up" = "#2b544d", "Down" = "#623747")) +
                theme_minimal() +
                guides(color = guide_legend(override.aes = list(size = 4))) +
                theme(legend.key.size = unit(1, 'cm'), legend.position = "top") +
                ylim(-max(abs(df.plot$bs.log2FoldChange)), max(abs(df.plot$bs.log2FoldChange))) +
                xlim(0,plot.gene.n.tiles) +
                labs(
                    title = paste0(plot.gene.name, " (", plot.geneID, ")"),
                    x = "Relative position on the gene",
                    y = "Fold-change (log2)",
                    fill = "Regulation",
                    color = "Regulation",
                    shape = "Transcript Region"
                )
            return(p)
        }

        # plot cases where region duplication is present
        if (any(all.regions.present %in% expected.regions)) {
            msg0 = paste0("At least one additional transcript region type detected compared to the currently supported expected regions: ",
                          paste(expected.regions, collapse = ", "),
                          ". Using transcript.regions.outlier.handle=", transcript.regions.outlier.handle," to solve. \n")
            if (!quiet) {
                warning(msg0)
            }

            if (transcript.regions.outlier.handle == "remove") {
                df.plot = mcols(plot.bs) %>% as.data.frame() %>%
                    filter(transcriptRegion %in% expected.regions) %>%
                    mutate(sig = (ifelse(bs.padj < alpha & abs(bs.log2FoldChange) > log2(lfc.cutoff), TRUE, FALSE))) %>%
                    mutate(dir = factor(ifelse(bs.log2FoldChange > 0, "Up", "Down"), levels = c("Up", "Down"))) %>%
                    mutate(sigDir = ifelse(sig == TRUE & dir == "Up", "Up", ifelse(sig == TRUE & dir == "Down", "Down", "Not"))) %>%
                    mutate(sigDir = factor(sigDir, levels = c("Not", "Up", "Down"))) %>%
                    mutate(duplicated = "no")

                p = ggplot(df.plot, aes(x = plot.position, y = bs.log2FoldChange, color = sigDir, fill = sigDir, shape = transcriptRegion)) +
                    geom_hline(yintercept = 0, color = "#4d4d4d") +
                    geom_segment(aes(x=plot.position, xend=plot.position, y=0, yend=bs.log2FoldChange), size = 1, color = "#4d4d4d") +
                    geom_point(size = 4, stroke = 1) +
                    scale_shape_manual(values=c("INTRON" = 23, "CDS" = 22, "UTR3" = 25, "UTR5" = 24, "OTHER" = 21)) +
                    scale_fill_manual(values = c("Not" = "#999999", "Up" = "#68b1a5", "Down" = "#874C62")) +
                    scale_color_manual(values = c("Not" = "#4d4d4d", "Up" = "#2b544d", "Down" = "#623747")) +
                    theme_minimal() +
                    guides(color = guide_legend(override.aes = list(size = 4))) +
                    theme(legend.key.size = unit(1, 'cm'), legend.position = "top") +
                    ylim(-max(abs(df.plot$bs.log2FoldChange)), max(abs(df.plot$bs.log2FoldChange))) +
                    xlim(0,plot.gene.n.tiles) +
                    labs(
                        title = paste0(plot.gene.name, " (", plot.geneID, ")"),
                        x = "Relative position on the gene",
                        y = "Fold-change (log2)",
                        fill = "Regulation",
                        color = "Regulation",
                        shape = "Transcript Region",
                        caption = paste0("Not all binding sites might be shown due to: transcript.regions.outlier.handle=", transcript.regions.outlier.handle)
                    )
                return(p)
            }

            if (transcript.regions.outlier.handle == "first") {
                df.duplicated = mcols(plot.bs) %>% as.data.frame() %>% filter(grepl(",", transcriptRegion)) %>% mutate(duplicated = "yes") %>%
                    mutate(transcriptRegion = sapply(strsplit(transcriptRegion,","), `[`, 1))
                df.single = mcols(plot.bs) %>% as.data.frame() %>% filter(!grepl(",", transcriptRegion)) %>% mutate(duplicated = "no")
                df.plot = rbind.data.frame(df.duplicated, df.single) %>%
                    mutate(sig = (ifelse(bs.padj < alpha & abs(bs.log2FoldChange) > log2(lfc.cutoff), TRUE, FALSE))) %>%
                    mutate(dir = factor(ifelse(bs.log2FoldChange > 0, "Up", "Down"), levels = c("Up", "Down"))) %>%
                    mutate(sigDir = ifelse(sig == TRUE & dir == "Up", "Up", ifelse(sig == TRUE & dir == "Down", "Down", "Not"))) %>%
                    mutate(sigDir = factor(sigDir, levels = c("Not", "Up", "Down"))) %>%
                    mutate(transcriptRegion = stringr::str_trim(transcriptRegion, "left")) %>%
                    mutate(transcriptRegion = stringr::str_trim(transcriptRegion, "right"))

                p = ggplot(df.plot, aes(x = plot.position, y = bs.log2FoldChange, color = sigDir, fill = sigDir, shape = transcriptRegion)) +
                    geom_hline(yintercept = 0, color = "#4d4d4d") +
                    geom_segment(aes(x=plot.position, xend=plot.position, y=0, yend=bs.log2FoldChange), size = 1, color = "#4d4d4d") +
                    geom_point(size = 4, stroke = 1) +
                    scale_shape_manual(values=c("INTRON" = 23, "CDS" = 22, "UTR3" = 25, "UTR5" = 24, "OTHER" = 21)) +
                    scale_fill_manual(values = c("Not" = "#999999", "Up" = "#68b1a5", "Down" = "#874C62")) +
                    scale_color_manual(values = c("Not" = "#4d4d4d", "Up" = "#2b544d", "Down" = "#623747")) +
                    theme_minimal() +
                    guides(color = guide_legend(override.aes = list(size = 4))) +
                    theme(legend.key.size = unit(1, 'cm'), legend.position = "top") +
                    ylim(-max(abs(df.plot$bs.log2FoldChange)), max(abs(df.plot$bs.log2FoldChange))) +
                    xlim(0,plot.gene.n.tiles) +
                    labs(
                        title = paste0(plot.gene.name, " (", plot.geneID, ")"),
                        x = "Relative position on the gene",
                        y = "Fold-change (log2)",
                        fill = "Regulation",
                        color = "Regulation",
                        shape = "Transcript Region",
                        caption = paste0("Not all binding sites might be shown due to: transcript.regions.outlier.handle=", transcript.regions.outlier.handle)
                    )
                return(p)
            }

            if (transcript.regions.outlier.handle == "second") {
                df.duplicated = mcols(plot.bs) %>% as.data.frame() %>% filter(grepl(",", transcriptRegion)) %>% mutate(duplicated = "yes") %>%
                    mutate(transcriptRegion = sapply(strsplit(transcriptRegion,","), `[`, 2))
                df.single = mcols(plot.bs) %>% as.data.frame() %>% filter(!grepl(",", transcriptRegion)) %>% mutate(duplicated = "no")
                df.plot = rbind.data.frame(df.duplicated, df.single) %>%
                    mutate(sig = (ifelse(bs.padj < alpha & abs(bs.log2FoldChange) > log2(lfc.cutoff), TRUE, FALSE))) %>%
                    mutate(dir = factor(ifelse(bs.log2FoldChange > 0, "Up", "Down"), levels = c("Up", "Down"))) %>%
                    mutate(sigDir = ifelse(sig == TRUE & dir == "Up", "Up", ifelse(sig == TRUE & dir == "Down", "Down", "Not"))) %>%
                    mutate(sigDir = factor(sigDir, levels = c("Not", "Up", "Down"))) %>%
                    mutate(transcriptRegion = stringr::str_trim(transcriptRegion, "left")) %>%
                    mutate(transcriptRegion = stringr::str_trim(transcriptRegion, "right"))

                p = ggplot(df.plot, aes(x = plot.position, y = bs.log2FoldChange, color = sigDir, fill = sigDir, shape = transcriptRegion)) +
                    geom_hline(yintercept = 0, color = "#4d4d4d") +
                    geom_segment(aes(x=plot.position, xend=plot.position, y=0, yend=bs.log2FoldChange), size = 1, color = "#4d4d4d") +
                    geom_point(size = 4, stroke = 1) +
                    scale_shape_manual(values=c("INTRON" = 23, "CDS" = 22, "UTR3" = 25, "UTR5" = 24, "OTHER" = 21)) +
                    scale_fill_manual(values = c("Not" = "#999999", "Up" = "#68b1a5", "Down" = "#874C62")) +
                    scale_color_manual(values = c("Not" = "#4d4d4d", "Up" = "#2b544d", "Down" = "#623747")) +
                    theme_minimal() +
                    guides(color = guide_legend(override.aes = list(size = 4))) +
                    theme(legend.key.size = unit(1, 'cm'), legend.position = "top") +
                    ylim(-max(abs(df.plot$bs.log2FoldChange)), max(abs(df.plot$bs.log2FoldChange))) +
                    xlim(0,plot.gene.n.tiles) +
                    labs(
                        title = paste0(plot.gene.name, " (", plot.geneID, ")"),
                        x = "Relative position on the gene",
                        y = "Fold-change (log2)",
                        fill = "Regulation",
                        color = "Regulation",
                        shape = "Transcript Region",
                        caption = paste0("Not all binding sites might be shown due to: transcript.regions.outlier.handle=", transcript.regions.outlier.handle)
                    )
                return(p)
            }

            if (transcript.regions.outlier.handle == "both") {
                df.duplicated = mcols(plot.bs) %>% as.data.frame() %>% filter(grepl(",", transcriptRegion)) %>% mutate(duplicated = "yes") %>%
                    separate_longer_delim(transcriptRegion, delim = ",")
                df.single = mcols(plot.bs) %>% as.data.frame() %>% filter(!grepl(",", transcriptRegion)) %>% mutate(duplicated = "no")
                df.plot = rbind.data.frame(df.duplicated, df.single) %>%
                    mutate(sig = (ifelse(bs.padj < alpha & abs(bs.log2FoldChange) > log2(lfc.cutoff), TRUE, FALSE))) %>%
                    mutate(dir = factor(ifelse(bs.log2FoldChange > 0, "Up", "Down"), levels = c("Up", "Down"))) %>%
                    mutate(sigDir = ifelse(sig == TRUE & dir == "Up", "Up", ifelse(sig == TRUE & dir == "Down", "Down", "Not"))) %>%
                    mutate(sigDir = factor(sigDir, levels = c("Not", "Up", "Down"))) %>%
                    mutate(transcriptRegion = stringr::str_trim(transcriptRegion, "left")) %>%
                    mutate(transcriptRegion = stringr::str_trim(transcriptRegion, "right"))

                p = ggplot(df.plot, aes(x = plot.position, y = bs.log2FoldChange, color = sigDir, fill = sigDir, shape = transcriptRegion)) +
                    geom_hline(yintercept = 0, color = "#4d4d4d") +
                    geom_segment(aes(x=plot.position, xend=plot.position, y=0, yend=bs.log2FoldChange), size = 1, color = "#4d4d4d") +
                    geom_point(size = 4, stroke = 1) +
                    ggrepel::geom_label_repel(data = df.plot %>% filter(duplicated == "yes"), aes(label = transcriptRegion)) +
                    scale_shape_manual(values=c("INTRON" = 23, "CDS" = 22, "UTR3" = 25, "UTR5" = 24, "OTHER" = 21)) +
                    scale_fill_manual(values = c("Not" = "#999999", "Up" = "#68b1a5", "Down" = "#874C62")) +
                    scale_color_manual(values = c("Not" = "#4d4d4d", "Up" = "#2b544d", "Down" = "#623747")) +
                    theme_minimal() +
                    guides(color = guide_legend(override.aes = list(size = 4))) +
                    theme(legend.key.size = unit(1, 'cm'), legend.position = "top") +
                    ylim(-max(abs(df.plot$bs.log2FoldChange)), max(abs(df.plot$bs.log2FoldChange))) +
                    xlim(0,plot.gene.n.tiles) +
                    labs(
                        title = paste0(plot.gene.name, " (", plot.geneID, ")"),
                        x = "Relative position on the gene",
                        y = "Fold-change (log2)",
                        fill = "Regulation",
                        color = "Regulation",
                        shape = "Transcript Region",
                        size = "Duplicated Region",
                        caption = paste0("Some binding sites might not be visible due to overplotting: transcript.regions.outlier.handle=", transcript.regions.outlier.handle)
                    )
                return(p)
            }

        }
    }
}


