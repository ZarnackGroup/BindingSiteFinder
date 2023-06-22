#' Plot signal coverage of selected ranges
#'
#' Function plots the coverage of the CLIP data in the signal slot and plots it
#' as coverage. The plot is centered around a given binding site, which can be
#' selected by an index.
#'
#' @param object a \link{BSFDataSet} object
#' @param plotIdx integer, index of the range to plot
#' @param flankPos numeric, number of nucleotides by which the plotting frame is
#' symmetrically extended
#' @param shiftPos numeric, nucleotide positions by which the current plotting
#' range should be shifted
#' @param mergeReplicates logical, if replicates should be merge per
#' condition (TRUE) or if every replicates should be shown separately (FALSE)
#' @param autoscale logical, if y-axis should be scaled to the maximum for all
#' replicates (TRUE), or not (FALSE)
#' @param highlight logical, if the central range should be highlighted (TRUE),
#' or not (FALSE)
#' @param showCentralRange logical, if the central range should be shown (TRUE),
#' or not (FALSE)
#' @param customRange \code{GenomicRanges}, a custom range object to be shown
#' underneath the coverage tracks
#' @param customRange.name character, the name of the customRange track
#' @param customAnnotation \code{GenomicRanges} or \code{TxDb}, a custom annotation
#' for eg. gene, exons, etc. to be shown underneath the coverage tracks
#' @param customAnnotation.name character, the name of the customAnnotation track
#' @param title character, set plotting title
#' @param colorPalette vector, hex colors used for the conditions
#'
#' @return an object of class \code{GVIZ}
#'
#' @importFrom GenomicRanges findOverlaps
#' @importFrom dplyr everything group_by summarize pull mutate rename
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#'
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' bindingSiteCoveragePlot(bds, plotIdx = 3, flankPos = 10)
#'
#' @export
bindingSiteCoveragePlot <- function(object, plotIdx, flankPos, shiftPos = NULL,
                                    mergeReplicates = FALSE, autoscale = FALSE,
                                    highlight = TRUE, showCentralRange = TRUE,
                                    customRange = NULL, customRange.name = "custom",
                                    customAnnotation = NULL, customAnnotation.name = "anno",
                                    title = NULL, colorPalette = NULL
) {

    # test parameter input validity
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is(mergeReplicates, "logical"))
    stopifnot(is(autoscale, "logical"))
    stopifnot(is(highlight, "logical"))
    stopifnot(is(showCentralRange, "logical"))

    if (!all(is.null(customRange) | is(customRange, "GenomicRanges"))) {
        stop("customRange not of the right type. Must be GenomicRanges or NULL")
    }

    if (!all(is.null(customAnnotation) |
             is(customAnnotation, "GenomicRanges") |
             is(customAnnotation, "TxDb"))) {
        stop("customAnnotation not of the right type. Must be GenomicRanges, TxDb or NULL")
    }

    # bind locally used variables
    `.` <- name <- value <- condition <- NULL

    # replicates should be merged
    if(isTRUE(mergeReplicates)) {
        object = collapseReplicates(object, collapseAll = FALSE)
    }

    # get initial objects for plotting
    rng = getRanges(object)
    sgn = getSignal(object)
    mta = getMeta(object)
    rngToPlot = rng[plotIdx]
    rngToPlotExtended = rng[plotIdx] + flankPos

    # shift window if needed
    if(!is.null(shiftPos)) {
        rngToPlotExtended = shift(rngToPlotExtended, shiftPos)
    }

    # find overlapping additional ranges
    idx = findOverlaps(rng, rngToPlotExtended)
    rngSub = rng[queryHits(idx)]

    if (is.null(names(rngSub))) {
        rngSub$symbol = paste0("BS:", seq_along(rngSub))
    } else {
        rngSub$symbol = names(rngSub)
    }

    # setting frame for plot
    currChr = unique(as.character(seqnames(rngToPlot)))
    currStart = start(rngToPlotExtended)
    currEnd = end(rngToPlotExtended)
    currStrand = as.character(strand(rngToPlotExtended))

    centerStart = start(rngToPlot)
    centerWidth = width(rngToPlot)-1

    # determine if window is too large for histogram boarder
    if (width(rngToPlotExtended) < 300) {
        plotRangeSmall = TRUE
    } else {
        plotRangeSmall = FALSE
    }

    # make tracks to plot
    axisTrack = Gviz::GenomeAxisTrack(fontcolor = "black")

    rngTrackOthers = Gviz::GeneRegionTrack(rngSub, chromosome = currChr,
                                     start = currStart, end = currEnd,
                                     name = "BS Other", fontsize = 6,
                                     exonAnnotation = "symbol", fontcolor.item="black",
                                     fill = "salmon", col = "black", col.line = "black",
                                     col.frame = "black", col.axis="black",
                                     col.border.title=NULL, frame = TRUE,
                                     col.title = "black", background.title="white")

    # compute coverage over all replicates
    sgnCov = suppressMessages(suppressWarnings(coverageOverRanges(setRanges(object, rngToPlotExtended),
                                                                  returnOptions = "merge_ranges_keep_positions")))

    # reverse signal if it is on the minus strand
    if (as.character(strand(rngToPlot)) == "-") {
        sgnCov = t(apply(sgnCov, 1, rev))
    }

    sgnRng = GenomicRanges::tile(rngToPlotExtended, width = 1) %>%
        unlist(., use.names = F)
    values(sgnRng) = t(sgnCov)

    if (isTRUE(autoscale)) {
        max = mcols(sgnRng) %>% as.data.frame() %>% max()
        max = rep(max, nrow(mta))
    } else {
        max = mcols(sgnRng) %>% as.data.frame() %>%
            tidyr::pivot_longer(dplyr::everything()) %>% dplyr::group_by(name) %>%
            dplyr::summarize(max = max(value)) %>% dplyr::pull(max)
    }

    # set color palette
    if (is.null(colorPalette)) {
        colorPalette = c("#A7D2CB", "#F2D388", "#C98474", "#874C62")
    }

    if (isTRUE(mergeReplicates)) {
        mta = data.frame(condition = as.factor(levels(mta$condition)))
    }

    if (length(levels(mta$condition)) <= length(colorPalette)) {
        replicatColors = lapply(seq_along(levels(mta$condition)), function(x){
            idx = which(levels(mta$condition)[x] == mta$condition)
            col = rep(colorPalette[x], length(idx))
            return(col)
        })
        names(replicatColors) = levels(mta$condition)
        if (!isTRUE(mergeReplicates)) {
            replicatColors = unlist(replicatColors, use.names = TRUE) %>% as.data.frame() %>%
                tibble::rownames_to_column("condition") %>%
                dplyr::mutate(condition = substring(condition, 1, nchar(condition)-1)) %>%
                dplyr::rename('col' = '.')
        }
        if (isTRUE(mergeReplicates)) {
            replicatColors = unlist(replicatColors, use.names = TRUE) %>% as.data.frame() %>%
                tibble::rownames_to_column("condition") %>%
                dplyr::rename('col' = '.')
        }
        mta$color = replicatColors$col[match(mta$condition, replicatColors$condition)]
    }
    if (length(levels(mta$condition)) == 1) {
        mta$color = "#A7D2CB"
    }
    if (length(levels(mta$condition)) > length(colorPalette)) {
        warning("The color palette supports only up to 4 different conditions at the moment. Setting all colors to grey.")
        mta$color = "#808080"
    }

    sgnTracks = lapply(seq_along(mcols(sgnRng)), function(x){
        currSample = sgnRng[,x]
        histOutlineColor = ifelse(plotRangeSmall, "black", mta$color[x])
        dTrack = Gviz::DataTrack(currSample, type = 'hist', ylim = c(0,max[x]),
                           col.histogram=histOutlineColor, fill.histogram = mta$color[x],
                           name = names(mcols(currSample)),
                           col = "black", col.line = "black",
                           col.frame = "black", col.axis="black",
                           col.border.title=NULL, frame = FALSE,
                           col.title = "black", background.title="white")
        return(dTrack)
    })

    # all main tracks are created at this point
    if (!isTRUE(showCentralRange)) {
        # remove central range from plotting if wanted
        allTracks = c(axisTrack, sgnTracks)
    } else {
        # all main tracks are added -> this is the default look of the plot
        allTracks = c(axisTrack, sgnTracks, rngTrackOthers)
    }

    # add additional custom range to be shown under the binding sites
    if (!is.null(customRange)) {
        # find overlapping custom ranges
        idx = findOverlaps(customRange, rngToPlotExtended)
        customRangeSub = customRange[queryHits(idx)]
        rngCustom = Gviz::GeneRegionTrack(customRangeSub, chromosome = currChr,
                                    start = currStart, end = currEnd,
                                    name = customRange.name, fontsize = 6,
                                    fill = "#B7C4CF", col = "black", col.line = "black",
                                    col.frame = "black", col.axis="black",
                                    col.border.title=NULL, frame = TRUE,
                                    col.title = "black", background.title="white")
        allTracks = c(allTracks, rngCustom)
    }
    # add additional gene annotation to be shown underneath
    if (!is.null(customAnnotation)) {
        rngGene = Gviz::GeneRegionTrack(customAnnotation, chromosome = currChr,
                                  start = currStart, end = currEnd,
                                  name = customAnnotation.name, fontsize = 6,
                                  fill = "#D7C0AE", col = "black", col.line = "black",
                                  col.frame = "black", col.axis="black",
                                  col.border.title=NULL, frame = TRUE,
                                  col.title = "black", background.title="white")
        allTracks = c(allTracks, rngGene)
    }
    # add highlighting option to the central range
    if(isTRUE(highlight)) {
        allTracks = Gviz::HighlightTrack(trackList = allTracks,
                                   start = centerStart, width = centerWidth, chromosome = currChr)
    }

    # manage plot title
    if (is.null(title)) {
        if (is.null(names(rngToPlot))) {
            title = paste0("BS:", as.character(plotIdx), "/ ",
                           currChr, ": ", currStart, "-", currEnd, " ", currStrand)
        } else {
            title = paste0("BS:", as.character(names(rngToPlot)), "/ ",
                           currChr, ": ", currStart, "-", currEnd, " ", currStrand)
        }
    }


    Gviz::plotTracks(trackList = allTracks,
               from = currStart, to = currEnd, main = title, cex.main = 1, col = "black")
}


