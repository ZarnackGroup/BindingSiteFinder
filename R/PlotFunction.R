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
#' @importFrom matrixStats colSds
#' @import ggplot2 GenomicRanges
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
                    fill = "lightblue"
                ) +
                geom_line(
                    data = df,
                    aes(x = position, y = mean),
                    size = 1,
                    color = "blue"
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
                    fill = "lightblue"
                ) +
                geom_line(
                    data = df,
                    aes(x = position, y = mean),
                    size = 1,
                    color = "blue"
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
#' @seealso \code{\link{makeBindingSites}}
#'
#' @import ggplot2 GenomicRanges
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
#' @importFrom tidyr pivot_longer
#' @import ggplot2 GenomicRanges
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
#' # use a single quantile cutoff
#' reproducibiliyCutoffPlot(bds, max.range = 20, cutoff = c(0.05))
#'
#' # use condition specific quantile cutoffs
#' reproducibiliyCutoffPlot(bds, max.range = 20, cutoff = c(0.1, 0.05))
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

        stopifnot(is(object, "BSFDataSet"))

        cond = getMeta(object)$condition
        df = coverageOverRanges(
            object, returnOptions = "merge_positions_keep_replicates",
            silent = TRUE)
        df = as.data.frame(mcols(df))

        if (length(cutoff) == 1) {
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

            rectDf = data.frame(
                condition = cond,
                name = unique(df$name),
                value = 1
            )

            p = ggplot(df, aes(
                x = value,
                group = name,
                fill = condition
            ), ...) +
                geom_rect(
                    data = rectDf,
                    aes(
                        xmin = -Inf,
                        xmax = Inf,
                        ymin = -Inf,
                        ymax = Inf,
                        color = condition
                    ),
                    alpha = 0.3
                ) +
                scale_fill_brewer(palette = "Dark2") +
                scale_color_brewer(palette = "Dark2") +
                geom_bar(fill = "lightblue", color = "blue") +
                facet_wrap( ~ name) +
                geom_vline(
                    data = qSel,
                    aes(xintercept = value, group = name),
                    color = "darkgrey",
                    size = 1.5
                ) +
                theme_classic() +
                ggtitle(paste(paste0(
                    "Cutoff: ", unique(qSel$per), " ",
                    unique(qSel$sel)
                ), collapse = "; ")) +
                xlab("Crosslinks per site") +
                ylab("Number of sites")
        }

        if (length(cutoff) > 1) {
            if (length(levels(cond)) == 1) {
                stop("multiple cutoffs are given but only one condition exists")
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

            p = ggplot(df, aes(
                x = value,
                group = name,
                fill = condition
            ), ...) +
                geom_rect(
                    data = rectDf,
                    aes(
                        xmin = -Inf,
                        xmax = Inf,
                        ymin = -Inf,
                        ymax = Inf,
                        color = condition
                    ),
                    alpha = 0.3
                ) +
                scale_fill_brewer(palette = "Dark2") +
                scale_color_brewer(palette = "Dark2") +
                geom_bar(fill = "lightblue", color = "blue") +
                facet_wrap( ~name) +
                geom_vline(
                    data = qSel,
                    aes(
                        xintercept = value,
                        group = name,
                        color = applyTo
                    ),
                    size = 1.5
                ) +
                theme_classic() +
                ggtitle(paste(vapply(levels(cond), function(x) {
                    paste0("Cutoff: ",
                           unique(qSel$per[qSel$applyTo == x]), " ", x)
                }, FUN.VALUE = character(1)),
                collapse = "; ")) +
                xlab("Crosslinks per site") +
                ylab("Number of sites")
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
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' supportRatioPlot(bds, bsWidths = c(3,7),
#' minWidth = 1, minClSites = 1, minCrosslinks = 2)
#'
#' @export
supportRatioPlot <- function(object, bsWidths, bsFlank = NA, ...){
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
