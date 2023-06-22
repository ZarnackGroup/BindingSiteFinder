

#' Plot the PureCLIP score distribution with global cutoff indicator
#'
#' A diagnostic function that plots the PureCLIP score distribution on a log2
#' scale. The function \code{\link{pureClipGlobalFilter}} is expected to be
#' executed prior to calling this plot function.
#'
#' @param object a \code{\link{BSFDataSet}} object
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{pureClipGlobalFilter}}
#'
#' @import ggplot2
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' # apply 5% filter
#' bds = pureClipGlobalFilter(object = bds, cutoff = 0.05)
#' pureClipGlobalFilterPlot(bds)
#'
#' @export
pureClipGlobalFilterPlot <- function(object) {

    # bind locally used variables
    x <- y <- quant <- name <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    if (is.null(object@params$pureClipGlobalFilter)) {
        msg0 = paste0("Global filter was not applied yet. Run BSFind() or pureClipGlobalFilter() to compute values first. \n")
        stop(msg0)
    }
    if (is.null(object@plotData$pureClipGlobalFilter$data) |
        is.null(object@plotData$pureClipGlobalFilter$cutoffPoint)) {
        msg1 = paste0("It seems like someting went wrong with your data. Please check the input and make sure to run pureClipGlobalFilter() or BSFind() before plotting. \n")
        stop(msg1)
    }

    optstr = object@params$pureClipGlobalFilter
    optstrNice = paste0("Cutoff=", optstr$cutoff, ", scoreMatchCol=", optstr$match.score)

    # get stored plotting data
    dt = data.frame(s = log2(object@plotData$pureClipGlobalFilter$data))
    cutpoint = log2(object@plotData$pureClipGlobalFilter$cutoffPoint)

    # calc density
    df = density(dt$s)
    df = data.frame(x = df$x, y = df$y)

    # set quantiles to show
    probs = seq(from = 0, to = 1, by = 0.1)
    quants = quantile(dt$s, prob = probs)
    df$quant = factor(findInterval(df$x, quants))

    # set cutoff to show
    cutoffNice = paste0(round(optstr$cutoff, 2) * 100, "%")
    dfCutoff = data.frame(pos = cutpoint, name = paste0(round(optstr$cutoff, 2) * 100, "%"))

    # Setup color pal
    colScaleFun = colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))

    p = ggplot(df, aes(x = x, y = y)) +
        geom_line() +
        geom_ribbon(aes(ymin = 0, ymax = y, fill = quant)) +
        scale_x_continuous(breaks = quants) +
        scale_fill_manual(values = colScaleFun(length(probs)+1)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6), legend.position = "none") +
        geom_vline(xintercept = cutpoint, linetype = "dashed", size = 1, color = "#b35900") +
        labs(
            title = "pureClipGlobalFilterPlot()",
            subtitle = optstrNice,
            x = "Score [log2]", y = "Density"
        ) +
        geom_label(data = dfCutoff, aes(x = pos, y = Inf, label = name), vjust = 1.2, color = "#b35900")
    return(p)
}


#' Plot the number of overlaps when assigning crosslink sites to genes
#'
#' A diagnostic function that plots the number of crosslink sites with their
#' respective overlapping rate. The function \code{\link{pureClipGeneWiseFilter}}
#' is expected to be executed prior to calling this plot function.
#'
#' @param object a \code{\link{BSFDataSet}} object
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{pureClipGeneWiseFilter}}
#'
#' @import ggplot2
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' # Load GRanges with genes
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#' # apply 5% gene-wise filter
#' bds = pureClipGeneWiseFilter(object = bds, anno.genes = gns, cutoff = 0.5,
#'  overlaps = "keepSingle")
#' duplicatedSitesPlot(bds)
#'
#' @export
duplicatedSitesPlot <- function(object) {

    # bind locally used variables
    `#N overlaps` <- Freq <- FreqNice <- nRes <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    if (is.null(object@params$geneWiseFilter)) {
        msg0 = paste0("Gene-wise filter was not applied yet. Run BSFind() or pureClipGeneWiseFilter() to compute values first. \n")
        stop(msg0)
    }
    if (is.null(object@plotData$geneWiseFilter$data) ) {
        msg1 = paste0("It seems like someting went wrong with your data. Please check the input and make sure to run pureClipGeneWiseFilter() or BSFind() before plotting. \n")
        stop(msg1)
    }

    optstr = object@params$geneWiseFilter
    optstrNice = paste0("Cutoff=", optstr$cutoff, ", overlaps=", optstr$overlaps)

    df = object@plotData$geneWiseFilter$data
    df = df %>% mutate(FreqNice = format(Freq, big.mark = ",", small.mark = "."))

    # Setup color pal
    cols = c("#576F72", rep("lightgrey", nrow(df)-1))

    # Make Plot
    # --------------------------------------------------------------------------
    # Standard plot
    p = ggplot(df, aes(x = `#N overlaps`, y = Freq, fill = `#N overlaps`, label = FreqNice)) +
        geom_col(color = "black") +
        scale_y_log10() +
        scale_fill_manual(values = cols) +
        theme_bw() +
        theme(legend.position = "none") +
        labs(
            title = "dupliatedSitesPlot()",
            subtitle = optstrNice,
            x = "#N overlaps",
            y = "#N (log10)"
        )

    # Add custom options based on overlappingLociOption
    if (optstr$overlaps == "keepSingle") {
        dfResolved = df %>% mutate(nRes = as.numeric(Freq) / as.numeric(`#N overlaps`))
        p = p +
            geom_col(data = dfResolved, aes(x = `#N overlaps`, y = nRes, fill = `#N overlaps`), fill = "#576F72", color = "black") +
            geom_text(data = df, aes(x = `#N overlaps`, y = Freq, label = FreqNice),
                      color = "#b35900", vjust = 1, hjust = 0.8, fontface = "bold", angle = 45)
    }
    if (optstr$overlaps == "removeAll") {
        dfResolved = df %>% mutate(nRes = ifelse(df$`#N overlaps` == 1, TRUE, FALSE))
        p = p +
            geom_col(data = subset(dfResolved, dfResolved$nRes == FALSE),
                     aes(x = `#N overlaps`, y = Freq, fill = `#N overlaps`),
                     fill = "#ECB390", color = "black") +
            geom_text(data = df, aes(x = `#N overlaps`, y = Freq, label = FreqNice),
                      color = "#b35900", vjust = 1, hjust = 0.5, fontface = "bold", angle = 45)
    }
    if (optstr$overlaps == "keepAll") {
        p = p +
            geom_col(data = df, aes(x = `#N overlaps`, y = Freq, fill = `#N overlaps`),
                     fill = "#576F72", color = "black") +
            geom_text(data = df, aes(x = `#N overlaps`, y = Freq, label = FreqNice),
                      color = "#b35900", vjust = 1, hjust = 0.5, fontface = "bold", angle = 45)
    }

    return(p)
}


#' Plot binding site merging diagnostics
#'
#' A diagnostic function that plots the number of regions to merge over the width
#' of these regions for each merging iteration calculated in
#' \code{\link{makeBindingSites}}. The function \code{\link{makeBindingSites}}
#' is expected to be executed prior to calling this plot function.
#'
#' @param object a \code{\link{BSFDataSet}} object
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{makeBindingSites}}
#'
#' @import ggplot2
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' bds = makeBindingSites(object = bds, bsSize = 9)
#' mergeCrosslinkDiagnosticsPlot(bds)
#'
#' @export
mergeCrosslinkDiagnosticsPlot <- function(object) {

    # bind locally used variables
    w <- s <- iteration <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))

    if (is.null(object@params$makeBindingSites)) {
        msg0 = paste0("Reproducibility filter was not applied yet. Run BSFind() or makeBindingSites() to compute values first. \n")
        stop(msg0)
    }
    if (is.null(object@plotData$makeBindingSites$mergeCsData) |
        is.null(object@plotData$makeBindingSites$data)) {
        msg1 = paste0("It seems like someting went wrong with your data. Please check the input and make sure to run makeBindingSites() or BSFind() before plotting. \n")
        stop(msg1)
    }

    optstr = object@params$makeBindingSites
    optstrNice = paste0("bsSize=", optstr$bsSize, ", minWidth=", optstr$minWidth, ", minCrosslinks=", optstr$minCrosslinks)

    dfPlot = object@plotData$makeBindingSites$mergeCsData
    dfPlot = dfPlot %>% mutate(iteration = as.factor(iteration), w = as.factor(w))

    # Make Plot
    # --------------------------------------------------------------------------
    # Standard plot
    cols = c("#716F81","#B97A95","#F6AE99","#F2E1C1","#BEAEE2","#F7DBF0","#CDF0EA","#F9F9F9",rep("grey",5))

    p = ggplot(dfPlot, aes(x = w, y = s, color = iteration, group = iteration)) +
        geom_point(size = 3) +
        geom_line(size = 1) +
        scale_color_manual(values = cols) +
        theme_bw() +
        labs(
            title = "mergeCrosslinkDiagnosticsPlot()",
            subtitle = optstrNice,
            x = "Region-to-fit width",
            y = "Region-to-fit #N",
            color = "Fitting iteration"
        ) +
        theme(legend.position = "bottom") +
        guides(colour = guide_legend(nrow = 1))

    return(p)
}

#' Plot binding site filter diagnostics
#'
#' A diagnostic function that plots the number of binding sites retained
#' after each filtering step calculated in \code{\link{makeBindingSites}}.
#' The function \code{\link{makeBindingSites}} is expected to be executed prior
#' to calling this plot function.
#'
#' @param object a \code{\link{BSFDataSet}} object
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{makeBindingSites}}
#'
#' @import ggplot2
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' bds = makeBindingSites(object = bds, bsSize = 9)
#' makeBsSummaryPlot(bds)
#'
#' @export
makeBsSummaryPlot <- function(object) {

    # bind locally used variables
    Option <- nRanges <- slice <- nNice <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))

    if (is.null(object@params$makeBindingSites)) {
        msg0 = paste0("Reproducibility filter was not applied yet. Run BSFind() or makeBindingSites() to compute values first. \n")
        stop(msg0)
    }
    if (is.null(object@plotData$makeBindingSites$mergeCsData) |
        is.null(object@plotData$makeBindingSites$data)) {
        msg1 = paste0("It seems like someting went wrong with your data. Please check the input and make sure to run makeBindingSites() or BSFind() before plotting. \n")
        stop(msg1)
    }

    optstr = object@params$makeBindingSites
    optstrNice = paste0("bsSize=", optstr$bsSize, ", minWidth=", optstr$minWidth, ", minCrosslinks=", optstr$minCrosslinks)

    dfPlot = object@plotData$makeBindingSites$data
    dfPlot = dfPlot %>%
        mutate(Option = factor(Option, levels = c(Option))) %>%
        slice(-1) %>%
        mutate(nNice = paste0(format(nRanges, big.mark = ",", decimal.mark = ".")))

    # Make Plot
    # --------------------------------------------------------------------------
    # Standard plot
    cols = c("#262A56","#8D7B68","#A4907C","#C8B6A6","#F1DEC9")

    p = ggplot(dfPlot, aes(x = Option, y = nRanges, fill = Option)) +
        geom_col() +
        scale_fill_manual(values = cols) +
        theme_bw() +
        theme(legend.position = "none") +
        labs(
            title = "makeBsSummaryPlot()",
            subtitle = optstrNice,
            x = "Options",
            y = "Ranges #N"
        ) +
        geom_text(data = dfPlot, aes(x = Option, y = nRanges, label = nNice),
                  color = "#b35900", vjust = 1, hjust = 0.8, fontface = "bold", angle = 45)

    return(p)
}


#' Plot to that shows the crosslink site distribution per replicate
#'
#' A diagnostic function that plots the number of crosslinks sites over the number
#' of crosslink in these sites and highlights the replicate specific
#' reproducibility cutoff that is derived from that distribution.
#' The function \code{\link{reproducibilityFilter}} is expected to be executed prior
#' to calling this plot function.
#'
#' @param object a \code{\link{BSFDataSet}} object
#' @param plotRange numeric; number of crosslinks per sites that should be shown
#' before summing them up
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{reproducibilityFilter}}
#'
#' @import ggplot2
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' bds = makeBindingSites(object = bds, bsSize = 9)
#' bds = reproducibilityFilter(bds)
#' reproducibilityFilterPlot(bds)
#'
#' @export
reproducibilityFilterPlot <- function(object, plotRange = 20) {

    # bind locally used variables
    value <- n <- condtion <- name <- condition <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))

    if (is.null(object@params$reproducibilityFilter)) {
        msg0 = paste0("Reproducibility filter was not applied yet. Run BSFind() or reproducibilityFilter() to compute values first. \n")
        stop(msg0)
    }
    if (is.null(object@plotData$reproducibilityFilterPlot$data) |
        is.null(object@plotData$reproducibilityFilterPlot$cutoffs) |
        is.null(object@plotData$reproducibilitySamplesPlot$data)) {
        msg1 = paste0("It seems like someting went wrong with your data. Please check the input and make sure to run reproducibilityFilter() or BSFind() before plotting. \n")
        stop(msg1)
    }

    optstr = object@params$reproducibilityFilter
    optstrNice = paste0("Cutoff=", optstr$cutoff, ", nReps=", optstr$nReps, ", minCrosslinks=", optstr$minCrosslinks)

    dfPlot = object@plotData$reproducibilityFilterPlot$data
    dfPlotCutoffs = object@plotData$reproducibilityFilterPlot$cutoffs

    # Apply plot range zoom in
    dfPlot$value[dfPlot$value > plotRange] = plotRange

    # Make Plot
    # --------------------------------------------------------------------------
    # Standard plot
    p = ggplot(dfPlot, aes(x = value, y = n)) +
        geom_col(color = "#2d5986", aes(fill = condition)) +
        facet_wrap(~name, scales = "free_x") +
        scale_fill_brewer(palette = "Set2") +
        theme_classic() +
        labs(
            title = "reproducibilityFilterPlot()",
            subtitle = optstrNice,
            x = "Crosslinks per site",
            y = "Number of sites",
            fill = "Condition"
        ) +
        geom_vline(data = dfPlotCutoffs, aes(xintercept = value, group = name), color = "#b35900", size = 1, linetype = "dashed") +
        geom_label(data = dfPlotCutoffs, aes(x = value, y = max(dfPlot$n), label = value), alpha = 0.5)
        # geom_label(data = dfPlotCutoffs, aes(x = value, y = 1, label = value), nudge_y = 2)

    return(p)
}


#' UpSet-plot to that shows how each replicate supports binding sites
#'
#' A diagnostic function that plots the set sizes for each replicate, indicating
#' how many binding site the specific replicate supports given its specific
#' threshold. The function \code{\link{reproducibilityFilter}} is expected to be
#' executed prior to calling this plot function.
#'
#' @param object a \code{\link{BSFDataSet}} object
#' @param nIntersections numeric; number of intersection to be shown
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{reproducibilityFilter}} \code{\link{reproducibilityFilterPlot}}
#'
#' @import ComplexHeatmap
#' @importFrom dplyr desc
#' @importFrom utils head tail
#' @importFrom grid gpar
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' bds = makeBindingSites(object = bds, bsSize = 9)
#' bds = reproducibilityFilter(bds)
#' reproducibilitySamplesPlot(bds)
#'
#' @export
reproducibilitySamplesPlot <- function(object, nIntersections = 20) {

    # bind locally used variables
    title <- rowTitle <- Freq <- `.` <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))

    if (is.null(object@params$reproducibilityFilter)) {
        msg0 = paste0("Reproducibility filter was not applied yet. Run BSFind() or reproducibilityFilter() to compute values first. \n")
        stop(msg0)
    }
    if (is.null(object@plotData$reproducibilitySamplesPlot$data)) {
        msg1 = paste0("It seems like someting went wrong with your data. Please check the input and make sure to run reproducibilityFilter() or BSFind() before plotting. \n")
        stop(msg1)
    }

    optstr = object@params$reproducibilityFilter
    optstrNice = paste0("Cutoff=", paste(optstr$cutoff, collapse = ","),
                        ", nReps=", paste(optstr$nReps, collapse = ","),
                        ", minCrosslinks=", paste(optstr$minCrosslinks, collapse = ","))

    # Get full combination matrix
    df = object@plotData$reproducibilitySamplesPlot$data
    # colnames(df) = colnames(df) %>% as.data.frame() %>% separate(col = ".", sep = "_", into = c(NA, "Sample")) %>% pull(Sample)
    m = make_comb_mat(df)

    rowTitle = "All intersections"

    # Reduce combination matrix for plotting
    if (length(comb_size(m)) > nIntersections) {
        plottingCutoff = comb_size(m) %>%
            as.data.frame() %>% rename("Freq" = ".") %>% arrange(desc(Freq)) %>%
            head(.,nIntersections) %>% tail(1) %>% pull(Freq)
        m = m[comb_size(m) >= plottingCutoff]
        rowTitle = paste0(nIntersections, " largest intersections")
    }

    title = paste0("reproducibilitySamplesPlot()\n", optstrNice)

    UpSet(m, column_title = title,
          row_title = rowTitle,
          comb_order = order(comb_size(m), decreasing = TRUE),
          top_annotation = upset_top_annotation(m, add_numbers = TRUE,
                                                numbers_gp = gpar(col = "#4C6793", fontsize = 8, fontface = "italic")),
          right_annotation = upset_right_annotation(m, add_numbers = TRUE, gp = gpar(fill = "#0B2447"),
                                                    numbers_gp = gpar(col = "#0B2447", fontsize = 8, fontface = "italic"),
                                                    numbers_rot = 25),
          comb_col = "#4C6793", bg_col = "white", pt_size = unit(.5, "cm"),
          border = T, lwd = 2, bg_pt_col = "grey",
          set_order = seq_along(set_size(m)))

}


#' UpSet-plot to that shows the gene type overlaps
#'
#' A diagnostic function that plots the gene types of binding sites on overlapping
#' loci genes. The function \code{\link{assignToGenes}} is expected to be
#' executed prior to calling this plot function.
#'
#' @param object a \code{\link{BSFDataSet}} object
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{assignToGenes}} \code{\link{targetGeneSpectrumPlot}}
#'
#' @import ComplexHeatmap
#' @importFrom grid gpar
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' # Load GRanges with genes
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#' bds = makeBindingSites(object = bds, bsSize = 9)
#' bds = assignToGenes(bds, anno.genes = gns)
#' geneOverlapsPlot(bds)
#'
#' @export
geneOverlapsPlot <- function(object) {

    # bind locally used variables
    title <- rowTitle <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))

    if (is.null(object@params$assignToGenes)) {
        msg0 = paste0("Assignment to anno.genes was not applied yet. Run BSFind() or assigneToGenes() to compute values first. \n")
        stop(msg0)
    }
    if (is.null(object@plotData$assignToGenes$dataOverlaps) |
        is.null(object@plotData$assignToGenes$dataSpectrum) ) {
        msg1 = paste0("It seems like someting went wrong with your data. Please check the input and make sure to run assignToGenes() or BSFind() before plotting. \n")
        stop(msg1)
    }

    optstr = object@params$assignToGenes
    optstrNice = paste0("Source=", optstr$source, ", Overlaps=", optstr$overlaps,
                        ifelse(!rlang::is_empty(optstr$rule), paste0(", rule=", paste(optstr$rule, collapse = ">")), ""))

    # Get full combination matrix
    df = object@plotData$assignToGenes$dataOverlaps
    m = make_comb_mat(df)

    rowTitle = "All intersections"
    title = paste0("geneOverlapsPlot()\n", optstrNice)

    UpSet(m, column_title = title,
          row_title = rowTitle,
          comb_order = order(comb_size(m), decreasing = TRUE),
          top_annotation = upset_top_annotation(m, add_numbers = TRUE,
                                                numbers_gp = gpar(col = "#4C6793", fontsize = 8, fontface = "italic")),
          right_annotation = upset_right_annotation(m, add_numbers = TRUE, gp = gpar(fill = "#0B2447"),
                                                    numbers_gp = gpar(col = "#0B2447", fontsize = 8, fontface = "italic"),
                                                    numbers_rot = 25),
          comb_col = "#4C6793", bg_col = "white", pt_size = unit(.5, "cm"),
          border = T, lwd = 2, bg_pt_col = "grey",
          set_order = seq_along(set_size(m)))
}

#' Bar-chart to show the hosting gene types of binding sites
#'
#' A diagnostic function that plots the gene type of the hosting gene for
#' each binding site. The function \code{\link{assignToGenes}} is expected to be
#' executed prior to calling this plot function.
#'
#' @param object a \code{\link{BSFDataSet}} object
#' @param showNGroups numeric; the number of different gene types to show
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{assignToGenes}} \code{\link{geneOverlapsPlot}}
#'
#' @import ggplot2
#' @importFrom dplyr lead desc mutate arrange
#' @importFrom utils head
#' @importFrom tibble add_row
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#' bds = makeBindingSites(object = bds, bsSize = 9)
#' bds = assignToGenes(bds, anno.genes = gns)
#' targetGeneSpectrumPlot(bds)
#'
#' @export
targetGeneSpectrumPlot <- function(object, showNGroups = 5) {

    # bind locally used variables
    GeneType <- Freq <- FreqNice <- nGenes <- nBs <- geneType <- nice <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))

    if (is.null(object@params$assignToGenes)) {
        msg0 = paste0("Assignment to anno.genes was not applied yet. Run BSFind() or assigneToGenes() to compute values first. \n")
        stop(msg0)
    }
    if (is.null(object@plotData$assignToGenes$dataOverlaps) |
        is.null(object@plotData$assignToGenes$dataSpectrum) ) {
        msg1 = paste0("It seems like someting went wrong with your data. Please check the input and make sure to run assignToGenes() or BSFind() before plotting. \n")
        stop(msg1)
    }

    optstr = object@params$assignToGenes
    optstrNice = paste0("Source=", optstr$source, ", Overlaps=", optstr$overlaps)

    df = object@plotData$assignToGenes$dataSpectrum

    # ----
    # Apply topN cutoff
    if (nrow(df) < showNGroups) {
        msg = paste0("Only ", nrow(df), " different gene types, displaying ", nrow(df), " gene types, not ", showNGroups, ". \n")
        showNGroups = nrow(df)
        message(msg)
    }
    if (nrow(df) > showNGroups) {
        otherFreq = df %>% arrange(desc(nGenes), desc(nBs)) %>% lead(showNGroups-1)
        dfOther = data.frame(geneType = "Other",
                             nBs = sum(otherFreq$nBs, na.rm = TRUE),
                             nGenes = sum(otherFreq$nGenes, na.rm = TRUE))
        df = df %>% as.data.frame() %>% arrange(desc(nGenes)) %>% head(showNGroups-1)
        df = rbind.data.frame(df, dfOther)
        df = df %>% mutate(nice = paste0(
            "#G:",
            format(nGenes, big.mark = ",", small.mark = "."),
            " (#BS:",
            format(nBs, big.mark = ",", small.mark = "."),
            ")")) %>%
            mutate(geneType = forcats::fct_reorder(geneType, nGenes))
    }
    if (nrow(df) == showNGroups) {
        df = df %>% arrange(desc(nGenes), desc(nBs)) %>%
            mutate(nice = paste0(
                "#G:",
                format(nGenes, big.mark = ",", small.mark = "."),
                " (#BS:",
                format(nBs, big.mark = ",", small.mark = "."),
                ")")) %>%
            mutate(geneType = forcats::fct_reorder(geneType, nGenes))

    }
    # ----
    # Setup color pal
    cols = RColorBrewer::brewer.pal(n = showNGroups+2, name = "Greys")
    cols = cols[1:showNGroups+1]

    # Make Plot
    # --------------------------------------------------------------------------
    p = ggplot(df, aes(x = geneType, y = nGenes, fill = geneType, label = nice)) +
        geom_col(color = "black") +
        scale_y_log10() +
        coord_flip() +
        scale_fill_manual(values = cols) +
        theme_bw() +
        theme(legend.position = "none") +
        geom_text(color = "#b35900", hjust = 1, fontface = "bold") +
        labs(
            title = "targetGeneSpectrumPlot()",
            subtitle = optstrNice,
            x = "",
            y = "#N (log10)"
        )

    return(p)
}



#' UpSet-plot to that shows the transcript region overlaps
#'
#' A diagnostic function that plots the transcript regions of binding sites
#' on overlapping loci. The function \code{\link{assignToTranscriptRegions}} is expected to be
#' executed prior to calling this plot function.
#'
#' @param object a \code{\link{BSFDataSet}} object
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{assignToTranscriptRegions}} \code{\link{transcriptRegionSpectrumPlot}}
#'
#' @import ComplexHeatmap
#' @importFrom grid gpar
#' @importFrom dplyr desc
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[2])
#' bds = makeBindingSites(object = bds, bsSize = 9)
#' bds = assignToGenes(bds, anno.genes = gns)
#' bds = assignToTranscriptRegions(object = bds, anno.transcriptRegionList = regions)
#' transcriptRegionOverlapsPlot(bds)
#'
#' @export
transcriptRegionOverlapsPlot <- function(object) {

    # bind locally used variables
    title <- rowTitle <- Freq <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))

    if (is.null(object@params$assignToTranscriptRegions)) {
        msg0 = paste0("Reproducibility filter was not applied yet. Run BSFind() or assignToTranscriptRegions() to compute values first. \n")
        stop(msg0)
    }
    if (is.null(object@plotData$assignToTranscriptRegions$dataOverlaps) |
        is.null(object@plotData$assignToTranscriptRegions$dataSpectrum)) {
        msg1 = paste0("It seems like someting went wrong with your data. Please check the input and make sure to run assignToTranscriptRegions() BSFind() before plotting. \n")
        stop(msg1)
    }
    optstr = object@params$assignToTranscriptRegions
    optstrNice = paste0("Source=", optstr$source, ", overlaps=", optstr$overlaps,
                        ifelse(!rlang::is_empty(optstr$rule), paste0(", rule=", paste(optstr$rule, collapse = ">")), ""))

    # Get full combination matrix
    df = object@plotData$assignToTranscriptRegions$dataOverlaps
    m = make_comb_mat(df, mode = "distinct")

    rowTitle = "All intersections"
    title = paste0("transcriptRegionOverlapsPlot()\n", optstrNice)

    UpSet(m, column_title = title,
          row_title = rowTitle,
          comb_order = order(comb_size(m), decreasing = TRUE),
          top_annotation = upset_top_annotation(m, add_numbers = TRUE,
                                                numbers_gp = gpar(col = "#4C6793", fontsize = 8, fontface = "italic")),
          right_annotation = upset_right_annotation(m, add_numbers = TRUE, gp = gpar(fill = "#0B2447"),
                                                    numbers_gp = gpar(col = "#0B2447", fontsize = 8, fontface = "italic"),
                                                    numbers_rot = 25),
          comb_col = "#4C6793", bg_col = "white", pt_size = unit(.5, "cm"),
          border = T, lwd = 2, bg_pt_col = "grey",
          set_order = seq_along(set_size(m)))
}


#' Bar-chart to show the hosting transcript regions of binding sites
#'
#' A diagnostic function that plots the transcript regions of the hosting gene for
#' each binding site. The function \code{\link{assignToTranscriptRegions}} is expected to be
#' executed prior to calling this plot function.
#'
#' @param object a \code{\link{BSFDataSet}} object
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{assignToTranscriptRegions}} \code{\link{transcriptRegionOverlapsPlot}}
#'
#' @import ggplot2
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[2])
#' bds = makeBindingSites(object = bds, bsSize = 9)
#' bds = assignToGenes(bds, anno.genes = gns)
#' bds = assignToTranscriptRegions(object = bds, anno.transcriptRegionList = regions)
#' transcriptRegionSpectrumPlot(bds)
#'
#' @export
transcriptRegionSpectrumPlot <- function(object) {

    # bind locally used variables
    TranscriptRegion <- FreqNice <- Freq <-  NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))

    if (is.null(object@params$assignToTranscriptRegions)) {
        msg0 = paste0("Reproducibility filter was not applied yet. Run BSFind() or assignToTranscriptRegions() to compute values first. \n")
        stop(msg0)
    }
    if (is.null(object@plotData$assignToTranscriptRegions$dataOverlaps) |
        is.null(object@plotData$assignToTranscriptRegions$dataSpectrum)) {
        msg1 = paste0("It seems like someting went wrong with your data. Please check the input and make sure to run assignToTranscriptRegions() or BSFind() before plotting. \n")
        stop(msg1)
    }

    optstr = object@params$assignToTranscriptRegions
    optstrNice = paste0("Source=", ifelse(!is.null(optstr$anno.annoDB), "anno.annoDB", "anno.transcriptRegionList"), ", overlaps=", optstr$overlaps,
                        ifelse(! rlang::is_empty(optstr$rule), paste0(", rule=", paste(optstr$rule, collapse = ">")), ""))

    # Get full combination matrix
    df = object@plotData$assignToTranscriptRegions$dataSpectrum
    df = df %>% as.data.frame() %>% arrange(desc(Freq)) %>%
        mutate(FreqNice = format(Freq, big.mark = ",", small.mark = ".")) %>%
        mutate(TranscriptRegion = forcats::fct_reorder(TranscriptRegion, Freq)) %>%
        mutate(TranscriptRegion = tolower(TranscriptRegion))

    # Nice formatting of names
    df$TranscriptRegion = unlist(lapply(df$TranscriptRegion, function(x){
        ifelse(x == "utr3", "3'UTR",
               ifelse(x == "utr5", "5'UTR",
                      ifelse(x == "cds", "CDS", .capitalize(x))))
    }))

    # set color scheme
    cols = c("#A7D2CB", "#F2D388", "#C98474", "#874C62", "#576F72", "#B4CDE6", "#7D6E83", "#D0B8A8")
    cols = cols[1:nrow(df)]

    p = ggplot(df, aes(x = TranscriptRegion, y = Freq, fill = TranscriptRegion, label = FreqNice)) +
        geom_col(color = "darkgrey", size = 0.5) +
        theme_bw() +
        theme(legend.position = "none") +
        geom_text(color = "#b35900", hjust = 0.5, fontface = "bold") +
        scale_fill_manual(values = cols) +
        # scale_fill_brewer(palette = "Set2") +
        coord_flip(clip = "on", expand = TRUE) +
        labs(title = "transcriptRegionSpectrumPlot()",
             subtitle = optstrNice,
             x = "Transcript regions",
             y = "#N Binding sites")
    return(p)
}


#' Plot the PureCLIP score distribution after re-assignment
#'
#' A diagnostic function that plots the PureCLIP score distribution on a log2
#' scale after the re-assignment on binding site level.
#' The function \code{\link{annotateWithScore}} is expected to be executed
#' prior to calling this plot function.
#'
#' @param object a \code{\link{BSFDataSet}} object
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{annotateWithScore}}
#'
#' @import ggplot2
#' @importFrom stats density
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' bds1 = makeBindingSites(object = bds, bsSize = 9)
#' bds1 = annotateWithScore(bds1, match.ranges = getRanges(bds))
#' globalScorePlot(bds1)
#'
#' @export
globalScorePlot <- function(object) {

    # bind locally used variables
    x <- y <- quant <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    if (is.null(object@params$annotateWithScore)) {
        msg0 = paste0("Score transfer was not applied yet. Run BSFind() or annotateWithScore() to compute values first. \n")
        stop(msg0)
    }
    if (is.null(object@plotData$annotateWithScore$data) ) {
        msg1 = paste0("It seems like someting went wrong with your data. Please check the input and make sure to run annotateWithScore() or BSFind() before plotting. \n")
        stop(msg1)
    }

    optstr = object@params$annotateWithScore
    optstrNice = paste0("match.score=", optstr$match.score, ", match.option=", optstr$match.option)

    # get stored plotting data
    dt = data.frame(s = log2(object@plotData$annotateWithScore$data))

    # calc density
    df = density(dt$s)
    df = data.frame(x = df$x, y = df$y)

    # set quantiles to show
    probs = seq(from = 0, to = 1, by = 0.1)
    quants = quantile(dt$s, prob = probs)
    df$quant = factor(findInterval(df$x, quants))

    # Setup color pal
    colScaleFun = colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))

    p = ggplot(df, aes(x = x, y = y)) +
        geom_line() +
        geom_ribbon(aes(ymin = 0, ymax = y, fill = quant)) +
        scale_x_continuous(breaks = quants) +
        scale_fill_manual(values = colScaleFun(length(probs)+1)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
        labs(
            title = "globalScorePlot()",
            subtitle = optstrNice,
            x = "Score [log2]", y = "Density"
        )
    return(p)
}


#' Plot the signal-to-flank score for varying gene-wise filter and binding
#' site width
#'
#' A diagnostic function that plots the the signal-to-flank score as a mean
#' for each binding site width and gene-wise filter as indicated when executing
#' \code{\link{estimateBsWidth}}. Additionally a mean of means visualizes the
#' overall trend and a red line indicates the suggested optimal binding site
#' width. The function \code{\link{estimateBsWidth}} is expected to be
#' executed prior to calling this plot function.
#'
#' @param object a \code{\link{BSFDataSet}} object
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{estimateBsWidth}}
#'
#' @import ggplot2
#' @importFrom stats sd
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#' bds = estimateBsWidth(bds, anno.genes = gns, est.maxBsWidth = 19,
#'  geneResolution = "coarse", bsResolution = "coarse", est.subsetChromosome = "chr22")
#' estimateBsWidthPlot(bds)
#' @export
estimateBsWidthPlot <- function(object) {

    # bind locally used variables
    bsSize <- ms <- geneWiseFilter <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))

    if (is.null(object@params$estimateBsWidth) | is.null(object@params$bsSize) |
        is.null(object@params$geneFilter)) {
        msg0 = paste0("Binding site width was not estimated. Run estimateBsWidth() or BSFind() first. \n")
        stop(msg0)
    }
    if (is.null(object@plotData$estimateBsWidth$data)) {
        msg1 = paste0("It seems like someting went wrong with your data. Please check the input and make sure to run estimateBsWidth() or BSFind() before plotting. \n")
    }

    optstr = object@params$estimateBsWidth
    optstrNice = paste0("bsResolution=", optstr$bsResolution,
                        ", geneFilter_resolution=", optstr$geneResolution,
                        ", subsetChromosome=", optstr$est.subsetChromosome)

    est.bsSize = object@params$bsSize
    est.GeneFilter = object@params$geneFilter
    df = object@plotData$estimateBsWidth$data

    dfPlot = df %>%
        mutate(bsSize = as.factor(bsSize)) %>%
        mutate(geneWiseFilter = paste0((geneWiseFilter)*100, "%")) %>%
        mutate(geneWiseFilter = factor(geneWiseFilter, levels = c(paste0(seq(0,95, by = 5), "%"))))

    dfMean = dfPlot %>%
        group_by(bsSize) %>%
        summarise(ms = mean(score), sd = sd(score)) %>%
        mutate(geneWiseFilter = "mean")


    dfAnno = dfPlot %>%
        mutate(bsSize = as.numeric(as.character(bsSize))) %>%
        group_by(geneWiseFilter) %>%
        filter(bsSize == max(bsSize)) %>%
        mutate(bsSize = as.factor(bsSize))

    cols = c("#f8ede3", "#f4e8dd", "#f1e3d7", "#edddd1", "#e9d8cb", "#e6d3c5",
             "#e2cebf", "#dec9b9", "#dbc3b3", "#d7bead", "#d4b9a6", "#d0b4a0",
             "#ccaf9a", "#c9a994", "#c5a48e", "#c19f88", "#be9a82", "#ba957c",
             "#b68f76", "#b38a70", "#af856a")

    p = ggplot(dfMean, aes(x = bsSize, y = ms, group = geneWiseFilter)) +
        geom_line(data = dfPlot, aes(x = bsSize, y = score, group = geneWiseFilter, color = geneWiseFilter), size = 0.8) +
        geom_point(data = dfPlot, aes(x = bsSize, y = score, group = geneWiseFilter, color = geneWiseFilter), size = 2) +
        scale_colour_manual(values = cols) +
        geom_errorbar(aes(ymin = ms-sd, ymax = ms+sd), color = "#7D6E83", width = 0.2, size = 1) +
        geom_line(color = "#7D6E83", size = 1) +
        geom_point(size = 4, shape = 21, color = "#7D6E83", fill = "#7D6E83", stroke = 1) +
        theme_bw() +
        theme(legend.position = "bottom") +
        guides(color=guide_legend(nrow=2,byrow=TRUE)) +
        labs(
            title = "estimateBsWidthPlot()",
            subtitle = optstrNice,
            y = "Signal-to-flank score",
            x = "Binding site width",
            color = "Gene-wise filter"
        ) +
        geom_vline(xintercept = dfMean$bsSize[dfMean$bsSize == est.bsSize], color = "#b35900", linetype = "dashed", size = 1) +
        geom_label(data = dfMean[dfMean$bsSize == est.bsSize,], aes(x = Inf, y = Inf),
                   label = paste0("Size=", est.bsSize, "/ filter=", est.GeneFilter), vjust = 1, size = 3, hjust = 1, alpha = 0.3)
    return(p)
}


#' Step-wise flowchart plot
#'
#' An overview plot that shows all workflow functions that were executed on the
#' current object, with all input and output binding site numbers and major options
#' that were used. The function can be called at any time in the analysis. Most
#' optimal usage is after a full run of the wrapper function \code{\link{BSFind}}.
#'
#' @param object a \code{\link{BSFDataSet}} object
#'
#' @return a plot of type \code{\link{ggplot}}
#'
#' @seealso \code{\link{BSFind}}
#'
#' @import ggplot2
#' @importFrom dplyr n lead slice left_join
#' @importFrom utils tail
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[2])
#' bds = BSFind(bds, anno.genes = gns, anno.transcriptRegionList = regions,
#'  est.subsetChromosome = "chr22")
#' processingStepsFlowChart(bds)
#'
#' @export
processingStepsFlowChart <- function(object) {

    # bind locally used variables
    xmin <- ymin <- xmax <- ymax <- Type <- NULL
    y <- x <- Step <- from <- to <- s_e <- Options <- `#N In` <- `#N Out` <- NULL
    y1 <- y2 <- id <- x1 <- x2 <- opt <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))

    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    # get results
    res = object@results
    res = format(res, big.mark = ",", decimal.mark = ".")
    colnames(res) = c("Step", "Type", "#N In", "#N Out", "%", "Options")

    flowNode = res %>%
        mutate(x = 0, y = n():1) %>%
        mutate(y = y *2) %>%
        mutate(xmin = x - 0.45,
               xmax = x + 0.45,
               ymin = y - 0.65,
               ymax = y + 0.65) %>%
        mutate(Type = as.factor(Type)) %>%
        as.data.frame()

    df = res %>%
        mutate(from = Step, to = lead(Step)) %>%
        slice(-n())

    flowEdges = df %>% dplyr::select(from, to) %>%
        mutate(id = row_number()) %>%
        pivot_longer(cols = c("from", "to"),
                     names_to = "s_e",
                     values_to = "Step") %>%
        left_join(flowNode, by = "Step") %>%
        dplyr::select(-c(Step, Type, y, xmin, xmax)) %>%
        mutate(y = ifelse(s_e == "from", ymin, ymax)) %>%
        dplyr::select(-c(ymin, ymax, Options))

    flowArrowIn = flowNode %>%
        mutate(x2 = xmin -0.01, x1 = xmin -0.25, y2 = ymax -0.45, y1 = ymax +0.4) %>%
        mutate(`#N In` = paste0("#N=", `#N In`),`#N Out` = paste0("#N=", `#N Out`))

    flowArrowOut = flowNode %>%
        mutate(x1 = xmin -0.01, x2 = xmin -0.25, y1 = ymax -0.7, y2 = ymax -1.5) %>%
        mutate(`#N In` = paste0("#N=", `#N In`),`#N Out` = paste0("#N=", `#N Out`))

    flowOutLast = flowArrowIn %>% tail(1) %>%
        mutate(y1 = y1-2, y2 = y2-2)

    suppressWarnings({
        flowOption = flowNode %>%
            separate(Options, sep = ",", into = c("opt")) %>% # throws a warning, that additional columns are dropped
            mutate(x2 = xmax+0.25, x1 = xmax +0.25, y2 = ymax -1, y1 = ymax -0.65)

        flowOptionArrow = flowNode %>%
            separate(Options, sep = ",", into = c("opt")) %>% # throws a warning, that additional columns are dropped
            mutate(x2 = xmax+0.01, x1 = xmax +0.25, y2 = ymax -0.85, y1 = ymax -0.65)
    })

    if (!is.null(object@params$estimateBsWidth)) {
        cols = c("#F2D7D9", "#D3CEDF", "#ECE5C7", "#9CB4CC")
    } else {
        cols = c("#F2D7D9", "#D3CEDF", "#9CB4CC")
    }

    p = ggplot() +
        geom_rect(data = flowNode, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = Type), color = "#4d4d4d") +
        geom_text(data = flowNode, aes(x = x, y = y, label = Step), color = "#585c45") +
        geom_path(data = flowEdges,
                  mapping = aes(x = x, y = y, group = id),
                  colour = "#585c45",
                  arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
        geom_curve(data = flowArrowIn, aes(x = x1, y = y1, xend = x2, yend = y2),
                   arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size = 0.4,
                   color = "gray20", curvature = -0.2) +
        geom_curve(data = flowArrowOut, aes(x = x1, y = y1, xend = x2, yend = y2),
                   arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size = 0.4,
                   color = "gray20", curvature = -0.2) +
        geom_label(data = flowArrowIn, aes(x = x1, y = y1, label = `#N In`), size = 4, hjust = 1, color = "#4d4d4d", fill = "#EDEDED") +
        geom_label(data = flowOutLast, aes(x = x1, y = y1, label = `#N Out`), size = 4, hjust = 1, color = "#4d4d4d", fill = "#EDEDED") +
        xlim(-1.3,1.3) +
        scale_fill_manual(values = cols) +
        theme_void() +
        labs(title = "BindingSiteFinder Flowchart",
             subtitle = "Use the function processingStepsTable(), with option=extended to get a table with the full overview of all options that \nwere used. ",
             fill = "") +
        theme(legend.position = "bottom",
              plot.title = element_text(hjust = 0.5, color = "#748DA6", size = 20, face = "bold"),
              plot.subtitle = element_text(color = "#748DA6", size = 8, face = "italic")) +
        geom_label(data = flowOption, aes(x = x1, y = y1, label = opt), size = 3, hjust = 0, color = "#4d4d4d", fill = "#EDEDED", fontface = "italic") +
        geom_curve(data = flowOptionArrow, aes(x = x1, y = y1, xend = x2, yend = y2),
                   arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size = 0.4,
                   color = "gray20", curvature = -0.2)

    return(p)
}


#' Plot crosslink events coverage over range
#'
#' A diagnostic plot function that allows to check the coverage of crosslink
#' events over different merged regions. The coverage is shown as mean over all
#' replicates and conditions, with a standard deviation corridor.
#'
#' If \code{object} is a single BSFDataObject a single coverage plot will be
#' drawn, whereas if it is a list of BSFDataObjects, then faceting is used to
#' make a plot for each list element.
#'
#' @param object a BSFDataSet, or a list of BSFDataSet
#' @param width a numeric value that defines the plotting ranges
#' @param name plot title
#' @param ... further arguments passed to ggplot
#'
#' @return a plot of type \code{ggplot2} displaying the crosslink coverage over
#' the ranges of the given \code{\link{BSFDataSet}}
#'
#' @seealso \code{\link{BSFDataSet}}
#'
#' @import ggplot2
#' @importFrom matrixStats colSds
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' # plotting a single object
#' bds <- makeBindingSites(object = bds, bsSize = 9, minWidth = 2,
#' minCrosslinks = 2, minClSites = 1)
#' rangeCoveragePlot(bds, width = 20)
#'
#' # plotting multiple objects
#' bds1 <- makeBindingSites(object = bds, bsSize = 3, minWidth = 2,
#' minCrosslinks = 2, minClSites = 1, sub.chr = "chr22")
#' bds2 <- makeBindingSites(object = bds, bsSize = 9, minWidth = 2,
#' minCrosslinks = 2, minClSites = 1, sub.chr = "chr22")
#' l = list(`1. bsSize = 3` = bds1, `2. bsSize = 9` = bds2)
#' rangeCoveragePlot(l, width = 20)
#'
#' @export
rangeCoveragePlot <-
    function(object, width, name = "Coverage Plot", ...) {
        # bind locally used variables
        position <- sd <- xmin <- xmax <- ymin <- ymax <- NULL

        # make plot if only a single object is given
        if (!is(object, "list")) {
            stopifnot(is(object, "BSFDataSet"))

            if (!is.numeric(width)) {
                stop("width has to be numeric")
            }

            bsSize = unique(width(getRanges(object)))
            rng = getRanges(object)
            rngResize = resize(rng, width = 1, fix = "center") + width
            newObj = setRanges(object, rngResize)

            df = coverageOverRanges(
                newObj, returnOptions = "merge_ranges_keep_positions",
                silent = TRUE)
            df = as.data.frame(df)

            df = data.frame(
                mean = colMeans(df),
                sd = colSds(as.matrix(df)),
                position = c(rev(-seq_len(width)), 0, seq_len(width))
            )
            p = ggplot(...) +
                annotate(
                    "rect",
                    xmin = (-floor(bsSize / 2)),
                    xmax = (floor(bsSize / 2)),
                    ymin = -Inf,
                    ymax = Inf,
                    alpha = 0.2,
                    color = "grey"
                ) +
                geom_ribbon(
                    data = df,
                    aes(
                        x = position,
                        ymin = mean - sd,
                        ymax = mean + sd
                    ),
                    fill = "#8cb3d9"
                ) +
                geom_line(
                    data = df,
                    aes(x = position, y = mean),
                    size = 1,
                    color = "#2d5986"
                ) +
                theme_classic() +
                ggtitle(name) +
                xlab("Position (nt)") +
                ylab("Crosslink events (mean of replicates)")

        }
        # make plot if a list of BSFDataSet objects is given
        if (is(object, "list")) {
            stopifnot(vapply(object, function(x) {
                is(x, "BSFDataSet")},
                FUN.VALUE = logical(1))
            )

            objectNames = names(object)
            df = lapply(seq_along(object), function(x) {
                rngResize = resize(getRanges(object[[x]]),
                                   width = 1,
                                   fix = "center") + width
                newObject = setRanges(object[[x]], rngResize)

                df = coverageOverRanges(
                    newObject, returnOptions = "merge_ranges_keep_positions",
                    silent = TRUE)
                df = as.data.frame(df)

                df = data.frame(
                    mean = colMeans(df),
                    sd = colSds(as.matrix(df)),
                    position = c(rev(-seq_len(width)), 0, seq_len(width)),
                    name = objectNames[[x]]
                )

                return(df)
            })
            df = do.call(rbind, df)
            bsSize = vapply(object, function(x) {
                unique(width(getRanges(x)))
            }, FUN.VALUE = integer(1))

            dfFrame = data.frame(
                name = names(bsSize),
                xmin = (-floor(bsSize / 2)),
                xmax = (floor(bsSize / 2)),
                ymin = -Inf,
                ymax = Inf
            )

            p = ggplot(...) +
                geom_rect(
                    data = dfFrame,
                    aes(
                        xmin = xmin,
                        xmax = xmax,
                        ymin = ymin,
                        ymax = ymax
                    ),
                    alpha = 0.2,
                    color = "grey"
                ) +
                geom_ribbon(
                    data = df,
                    aes(
                        x = position,
                        ymin = mean - sd,
                        ymax = mean + sd
                    ),
                    fill = "#8cb3d9"
                ) +
                geom_line(
                    data = df,
                    aes(x = position, y = mean),
                    size = 1,
                    color = "#2d5986"
                ) +
                theme_classic() +
                facet_wrap( ~ name) +
                xlab("Position (nt)") +
                ylab("Crosslink events (mean of replicates)")
        }
        return(p)
    }


#' Plot summarized results of the different binding site merging and filtering
#' steps
#'
#' Bar charts produced for the different filter steps in the binding site
#' merging routine. Depending on the selected option (\code{select}) all or
#' only a user defined filter can be shown.
#'
#' If \code{object} is a single BSFDataObject a single coverage plot will be
#' drawn, whereas if it is a list of BSFDataObjects, then faceting is used to
#' make a plot for each list element.
#'
#' @param object a BSFDataObject, with the makeBindingSites function already run
#' @param select one of "all", "filter", "inputRanges",
#' "minCLSites", "mergeCrosslinkSites", "minCrosslinks",
#' "centerIsClSite" or "centerIsSummit".
#' Defines which parameter is selected for plotting.
#' @param ... further arguments passed to ggplot
#'
#' @return a plot of type ggplot after the \code{\link{makeBindingSites}}
#' function was run
#'
#' @import ggplot2
#'
#' @seealso \code{\link{makeBindingSites}}
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' # plotting a single object
#' bds0 <- makeBindingSites(object = bds, bsSize = 9, minWidth = 2,
#' minCrosslinks = 2, minClSites = 1)
#' mergeSummaryPlot(bds0)
#'
#' # plotting mulitple obejcts
#' bds1 <- makeBindingSites(object = bds, bsSize = 9, minWidth = 2,
#' minCrosslinks = 2, minClSites = 1, sub.chr = "chr22")
#' bds2 <- makeBindingSites(object = bds, bsSize = 9, minWidth = 2,
#' minCrosslinks = 2, minClSites = 3, sub.chr = "chr22")
#' l = list(`1. bsSize = 3` = bds1, `2. bsSize = 9` = bds2)
#' rangeCoveragePlot(l, width = 20)
#'
#' @export
mergeSummaryPlot <- function(object,
                             select = c(
                                 "all",
                                 "filter",
                                 "inputRanges",
                                 "minClSites",
                                 "mergeCrosslinkSites",
                                 "minCrosslinks",
                                 "centerIsClSite",
                                 "centerIsSummit"
                             ),
                             ...) {
    # bind locally used variables
    Option <- nRanges <- name <- NULL

    # make plot if only a single summary is given
    if (!is(object, "list")) {
        stopifnot(is(object, "BSFDataSet"))
        smry = getSummary(object)
        # filter not possible if only a single summary is given
        select = "all"
        df = smry
        df$Option = factor(df$Option, levels = c(df$Option))
        p = ggplot(df, aes(
            x = Option,
            y = nRanges,
            fill = Option
        ), ...) +
            geom_col() +
            scale_fill_brewer(palette = "Dark2") +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 90),
                  legend.position = "none") +
            ggtitle("Processing steps")

    }
    # make plot if a list of summaries is give
    if (is(object, "list")) {
        select = match.arg(
            select,
            choices = c(
                "all",
                "filter",
                "inputRanges",
                "minClSites",
                "mergeCrosslinkSites",
                "minCrosslinks",
                "centerIsClSite",
                "centerIsSummit"
            )
        )

        stopifnot(all(vapply(object, function(o) {
            (is(o, "BSFDataSet"))},
            FUN.VALUE = logical(1))))

        smry = lapply(object, getSummary)
        df = do.call("rbind", smry)
        df$name = rep(names(smry), each = 6)
        df$Option = factor(df$Option, levels = c(smry[[1]][[1]]))
        if (select == "all") {
            p = ggplot(df, aes(
                x = Option,
                y = nRanges,
                fill = Option
            ), ...) +
                geom_col() +
                scale_fill_brewer(palette = "Dark2") +
                theme_classic() +
                theme(axis.text.x = element_text(angle = 90),
                      legend.position = "none") +
                facet_wrap( ~ name) +
                ggtitle("Processing steps")
        }
        if (select == "filter") {
            df = df[df$Option != c("inputRanges", "mergeCrosslinkSites"), ]
            p = ggplot(df, aes(
                x = Option,
                y = nRanges,
                fill = Option
            ), ...) +
                geom_col() +
                scale_fill_brewer(palette = "Dark2") +
                theme_classic() +
                theme(axis.text.x = element_text(angle = 90),
                      legend.position = "none") +
                facet_wrap( ~ name) +
                ggtitle("Filtering options")
        }
        if (select == "inputRanges" |
            select == "mergeCrosslinkSites" |
            select == "minCrosslinks" |
            select == "centerIsClSite" |
            select == "centerIsSummit" | select == "minClSites") {
            df = df[df$Option == c(select), ]
            p = ggplot(df, aes(
                x = name,
                y = nRanges,
                fill = name
            ), ...length()) +
                geom_col() +
                scale_fill_brewer(palette = "Dark2") +
                theme_classic() +
                theme(axis.text.x = element_text(angle = 90),
                      legend.position = "none") +
                ggtitle(select)
        }

    }

    return(p)
}

#' Plot to that shows how many replicates support each binding site
#'
#' Plotting function for settings specified in
#' \code{\link{reproducibilityFilter}}.
#'
#' @param object a BSFDataSet object
#' @param cutoff a vector of length = 1, or of length = levels(meta$conditions)
#' with a single number (between 0-1) indicating the quantile cutoff
#' @param min.crosslinks numeric of length = 1, defines the lower boundary for
#' the minimum number of crosslinks a binding site has to be supported by all
#' replicates, regardless of the replicate specific quantile threshold
#' @param max.range maximum number of crosslinks per sites that should be shown
#' @param ... further arguments passed to ggplot
#'
#' @return a plot of type \code{ggplot2} showing the per replicate
#' reproducibility cutoffs based on a given quantile threshold
#'
#' @seealso \code{\link{reproducibilityFilter}}
#'
#' @import ggplot2
#' @importFrom ggforce geom_mark_rect
#'
#' @examples
#'
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' # merge binding sites
#' bds <- makeBindingSites(object = bds, bsSize = 9, minWidth = 2,
#' minCrosslinks = 2, minClSites = 1)
#'
#' # use the same cutoff for both conditions
#' suppressWarnings(reproducibiliyCutoffPlot(bds, max.range = 20, cutoff = c(0.05, 0.05)))
#'
#' # use different cutoffs for each condition
#' suppressWarnings(reproducibiliyCutoffPlot(bds, max.range = 20, cutoff = c(0.1, 0.05)))
#'
#'
#' @export
reproducibiliyCutoffPlot <-
    function(object,
             cutoff = 0.05,
             min.crosslinks = 1,
             max.range = 20,
             ...) {

        # bind locally used variables
        value <- name <- condition <- applyTo <- NULL

        # deprecation notice
        .Deprecated("reproducibilityFilterPlot")

        stopifnot(is(object, "BSFDataSet"))

        cond = getMeta(object)$condition
        df = coverageOverRanges(
            object, returnOptions = "merge_positions_keep_replicates",
            silent = TRUE)
        df = as.data.frame(mcols(df))

        if (length(cutoff) == 1) {
            if(length(levels(cond)) > 1) {
                stop("Only one cutoff is given for multiple conditions.")
            }
            # calculate sample specific thresholds
            qSel = .selectQuantilesSingleCondtion(
                covDf = df,
                userCond = cond,
                userNreps = 1,
                userCutoff = cutoff
            )
            # apply minimal crosslink threshold
            qSel$value = ifelse(qSel$value < min.crosslinks,
                                qSel$value + min.crosslinks, qSel$value)

            # prepare objects for plotting
            df[df > max.range] = max.range
            df = df %>% pivot_longer(everything())
            df$condition = vapply(strsplit(df$name, "_"), `[`, 2,
                                  FUN.VALUE = character(1))

            qSel$lower = ifelse(qSel$value == min.crosslinks, TRUE, FALSE)

            p = ggplot(df,
                       aes(x = value, group = name)) +
                scale_fill_brewer(palette = "Set1") +
                scale_color_brewer(palette = "Set1") +
                geom_bar(fill = "#2d5986", color = "#2d5986") +
                facet_wrap( ~name, scales = "free_x") +
                geom_vline(data = qSel,
                           aes(xintercept = value, group = name),
                           color = "darkgrey", size = 1.5) +
                theme_classic() +
                geom_label(data = qSel, aes(x = value, y = 1, label = value)) +
                labs(
                    title = paste0(paste(paste0("Cutoff: ", unique(qSel$per), " ",
                                                unique(qSel$sel)), collapse = "; "),
                                   "; min.crosslinks = ", min.crosslinks),
                    x = "Crosslinks per site",
                    y = "Number of sites",
                    fill = "condition"
                )

        }
        if (length(cutoff) > 1) {
            if (length(levels(cond)) == 1) {
                stop("Multiple cutoffs are given but only one condition exists")
            }
            # calculate sample specific thresholds
            qSel = .selectQuantilesMultipleConditions(
                covDf = df,
                userCond = cond,
                userNreps = 1,
                userCutoff = cutoff
            )
            # apply minimal crosslink threshold
            qSel$value = ifelse(qSel$value < min.crosslinks,
                                qSel$value + min.crosslinks, qSel$value)

            # prepare objects for plotting
            df[df > max.range] = max.range
            df = df %>% pivot_longer(everything())
            df$condition = vapply(strsplit(df$name, "_"), `[`, 2,
                                  FUN.VALUE = character(1))
            idx = match(df$condition, unique(cond))
            df$cutoff = cutoff[idx]
            rectDf = data.frame(
                condition = cond,
                name = unique(df$name),
                value = 1
            )
            qSel$lower = ifelse(qSel$value == min.crosslinks, TRUE, FALSE)

            p = ggplot(df,
                       aes(x = value, group = name)) +
                geom_rect(data = rectDf,
                          aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                              fill = condition), alpha = 0.3) +
                scale_fill_brewer(palette = "Set1") +
                # scale_color_brewer(palette = "Set1") +
                geom_bar(fill = "#2d5986", color = "#2d5986") +
                facet_wrap( ~name, scales = "free_x") +
                geom_vline(data = qSel,
                           aes(xintercept = value, group = name,
                               color = applyTo), size = 1.5) +
                theme_classic() +
                geom_label(data = qSel, aes(x = value, y = 1, label = value)) +
                labs(
                    title = paste0(paste(paste0("Cutoff: ", unique(qSel$per), " ",
                                                unique(qSel$sel)), collapse = "; "),
                                   "; min.crosslinks = ", min.crosslinks),
                    x = "Crosslinks per site",
                    y = "Number of sites",
                    fill = "condition"
                ) +
                scale_color_discrete(guide="none")
        }
        return(p)
    }


#' Plot that shows the binding site support ratio
#'
#' Function that shows a ratio to determine how well a given binding
#' site with is supported by the crosslink coverage of the data.
#' Ratios are computed using the \code{\link{supportRatio}} function.
#'
#' The higher the ratio, the more does the given binding site width captures
#' the enrichment of crosslinks compared the the local surrounding. A ratio
#' equal to 1 would mean no enrichemnt at all.
#'
#' @param object a BSFDataSet object
#' @param bsWidths a numeric vector indicating the different binding site
#' width to compute the ratio for
#' @param bsFlank optional; a numeric vector of the same length as
#' \code{bsWidth} used to specify the width of the flanking regions
#' @param ... further arguments passed to \code{makeBindingSites}
#'
#' @return an object of class \code{ggplot2}
#'
#' @import ggplot2
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' suppressWarnings(supportRatioPlot(bds, bsWidths = c(3,7),
#' minWidth = 1, minClSites = 1, minCrosslinks = 2))
#'
#' @export
supportRatioPlot <- function(object, bsWidths, bsFlank = NA, ...){
    # deprecation notice
    .Deprecated("estimateBsWidthPlot")

    df = supportRatio(object, bsWidths, bsFlank, ...)

    p = ggplot(df, aes(x = bsWidths, y = supportRatio, group = 1)) +
        geom_point() +
        geom_line() +
        labs(
            title = "Binding site widths signal ratios",
            x = "Binding site widths",
            y = "Signal to flank ratio"
        ) +
        theme_bw()
    p
}

