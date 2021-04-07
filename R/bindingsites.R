#' Define equally sized binding sites from peak calling results and iCLIP crosslink events.
#'
#' This function performs the merging of single nucleotide crosslink sites into
#' binding sites of a user defined width (\code{bsSize}). Depending on the desired
#' output width crosslink sites with a distance closer than \code{bsSize} -1 are concatenated.
#' Initially all input regions are concatenated and then imperatively merged and extended.
#' Concatenated regions smaller than \code{minWidth} are removed prior to the merge
#' and extension routine. This prevents outlier crosslink pileup, eg. mapping artifacts
#' to be integrated into the final binding sites. All remaining regions are further
#' processed and regions larger than the desired output width are imperatively split up
#' by setting always the position with the highest number of crosslinks as center.
#' Regions smaller than the desired width are symmetrically extended.
#' Resulting binding sites are then filtered by the defined constraints.
#'
#' The \code{bsSize} argument defines the final output width of the merged binding sites.
#' I has to be an ood number, to ensure that a binding site has a distinct center.
#'
#' The \code{minWidth} parameter is used to describe the minimum width a ranges has
#' to be after the initial concatenation step. For example: Consider bsSize = 9 and
#' minWidht = 2. Then all initial crosslink sites that are closer to each other than
#' 8 nucleotides (bsSize -1) will be concatenated. Any of these ranges with less than
#' 2 nucleotides of width will be removed.
#'
#' The argument \code{minCrosslinks} defines how many single nucleotide resolition
#' crosslink events need to overlap a binding site to pass. A default of 2 is very inclusive
#' resulting in only the removal of the most extreme cases.
#'
#' The \code{minClSites} argument defines how many positions of the binding site
#' must have been covered by the original crosslink site input. If the input was
#' based on the single nucleotide crosslink positions computed by PureCLIP than
#' this filter checks for the number of positions originally identified by PureCLIP
#' in the computed binding sites. The default of \code{minClSites} = 1 essentially
#' deactivates this filter.
#'
#' The options \code{centerIsClSite} and \code{centerIsSummit} ensure that the center
#' of each binding site is covered by an initial crosslink site and represents the
#' summit of crosslink events in the binding site, respectively.
#'
#' The option \code{sub.chr} allows to run the binding site merging on a smaller subset
#' (eg. "chr1") for improoved computational speed when testing the effect of various
#' binding site width and filtering options.
#'
#' @param object a BSFDataSet object (see \code{\link{BSFDataSet}})
#' @param bsSize an odd integer value specifiying the size of the output binding sites
#' @param minWidth the minimum size of regions that are subjected to the iterative
#' merging routine, after the initial region concatenation.
#' @param minCrosslinks the minimal number of crosslink events that have to overlap a
#' final binding site
#' @param minClSites the minimal number of crosslink sites that have to overlap a
#' final binding site
#' @param centerIsClSite logical, whether the center of a final binding site must
#' be covered by an intial crosslink site
#' @param centerIsSummit logical, whether the center of a final binding site must
#' exhibit the highest number of crosslink events
#' @param sub.chr chromosome identifier (eg, chr1, chr2) used for subsetting the
#' BSFDataSet object. This option can be used for testing different parameter options
#'
#' @return an object of type BSFDataSet with modified ranges
#'
#' @import GenomicRanges methods
#'
#' @examples
#'
#' # load data
#' csFile <- system.file("extdata", "PureCLIP_crosslink_sites_example.bed",
#'  package="BindingSiteFinder")
#' cs = rtracklayer::import(con = csFile, format = "BED")
#' clipFiles <- system.file("extdata", package="BindingSiteFinder")
#'
#' # two experimental conditions
#' meta = data.frame(condition = factor(c("WT", "WT", "KD", "KD"), levels = c("KD", "WT")),
#' clPlus = list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE),
#' clMinus = list.files(clipFiles, pattern = "minus.bw$", full.names = TRUE))
#' bds = BSFDataSet(ranges = cs, meta = meta)
#'
#' # standard options, no subsetting
#' bds <- makeBindingSites(object = bds, bsSize = 9, minWidth = 2,
#' minCrosslinks = 2, minClSites = 1)
#'
#' # standard options, with subsetting
#' bds <- makeBindingSites(object = bds, bsSize = 9, minWidth = 2,
#' minCrosslinks = 2, minClSites = 1, sub.chr = "chr1")

#' @export
makeBindingSites <- function(object,
                             bsSize,
                             minWidth = 2,
                             minCrosslinks = 2,
                             minClSites = 1,
                             centerIsClSite = TRUE,
                             centerIsSummit = TRUE,
                             sub.chr = NA) {
    stopifnot(is(object, "BSFDataSet"))

    # check integrity of input values
    if (missing(bsSize)) {
        stop("bsSize not set. Please provide a desired binding site width.")
    }
    if (c(bsSize %% 2) == 0) {
        stop("bsSize is even. An odd number is required to have a distinct binding site center.")
    }
    if (!is.numeric(bsSize)) {
        stop("bsSize must be of type numeric")
    }
    if (!is.numeric(minWidth)) {
        stop("minWidth must be of type numeric")
    }
    if (!is.numeric(minClSites)) {
        stop("minClSites must be of type numeric")
    }
    if (!is.numeric(minCrosslinks)) {
        stop("minCrosslinks must be of type numeric")
    }

    #---------------------------------------------------------------------------
    # subsetting options
    if (!is.na(sub.chr)) {
        if (!is.character(sub.chr)) {
            stop("sub.chr must be of type character")
        }
        # subset the signal
        objectSub = .subsetByChr(object = object, chr = sub.chr)
        sgn = getSignal(objectSub)
        rngS0 = getRanges(objectSub)
    }
    # no subsetting
    if (is.na(sub.chr)) {
        sgn = getSignal(object)
        rngS0 = getRanges(object)
    }
    if (length(rngS0) == 0) {
        stop("0 ranges as input.")
    }

    #---------------------------------------------------------------------------
    # prepare data for merging
    sgnMerge = .collapesReplicates(sgn)

    # execute filter and merging routines
    rngS1 = .mergeCrosslinkSites(
        rng = rngS0,
        sgn = sgnMerge,
        bsSize = bsSize,
        minWidth = minWidth
    )
    stopifnot(length(rngS1) > 0)

    rngS2 = .filterMinCrosslinks(rng = rngS1,
                                 sgn = sgnMerge,
                                 minCrosslinks = minCrosslinks)
    stopifnot(length(rngS2) > 0)

    rngS3 = .filterMinClSites(rng = rngS2,
                              rng0 = rngS0,
                              minClSites = minClSites)
    stopifnot(length(rngS3) > 0)

    if (isTRUE(centerIsClSite)) {
        rngS4 = .filterCenterClSite(rng = rngS3, rng0 = rngS0)
    } else {
        rngS4 = rngS3
    }
    stopifnot(length(rngS4) > 0)

    if (isTRUE(centerIsSummit)) {
        rngS5 = .filterCenterSummit(rng = rngS4, sgn = sgnMerge)
    } else {
        rngS5 = rngS4
    }
    stopifnot(length(rngS5) > 0)

    # summarize number of ranges after each step
    reportDf = data.frame(
        Option = c(
            "inputRanges",
            "mergeCrosslinkSites",
            "minCrosslinks",
            "minClSites",
            "centerIsClSite",
            "centerIsSummit"
        ),
        nRanges = c(
            length(rngS0),
            length(rngS1),
            length(rngS2),
            length(rngS3),
            ifelse(isTRUE(centerIsClSite),
                   length(rngS4), NA),
            ifelse(isTRUE(centerIsSummit),
                   length(rngS5), NA)
        )
    )

    # update BSFDataSet with new ranges information
    objectNew = setRanges(object, rngS5)
    objectNew = setSignal(objectNew, sgn)
    objectNew = setSummary(objectNew, reportDf)
    ClipDS = objectNew
    return(ClipDS)
}


################################################################################
# unexported functions

.mergeCrosslinkSites <- function(rng,
                                 sgn,
                                 bsSize,
                                 minWidth) {
    # summarize signal over all replicates for mergeing
    sgnMergePlus = sgn$signalPlus
    sgnMergeMinus = sgn$signalMinus

    rngS1 = GenomeInfoDb::keepStandardChromosomes(rng, pruning.mode = "coarse")

    ### Merge peaks for given bs size
    rngS2 = reduce(rngS1, min.gapwidth = bsSize - 1)

    ### Remove peaks smaller than min width
    rngS3 = rngS2[width(rngS2) >= minWidth]
    names(rngS3) = seq_along(rngS3)

    ### Center detection and extension
    rngCenterPlus <- GRanges()
    rngCenterMinus <- GRanges()
    rngToProcessPlus <- subset(rngS3, strand == "+")
    rngToProcessMinus <- subset(rngS3, strand == "-")

    Counter = 0
    while (TRUE) {
        # quit if no more regions to check
        if (length(rngToProcessMinus) == 0 &
            length(rngToProcessPlus) == 0) {
            break
        } else {
            if (length(rngToProcessPlus) != 0) {
                # get max xlink position of each peak
                peaksMaxPosPlus = as.matrix(sgnMergePlus[rngToProcessPlus])
                peaksMaxPosPlus[is.na(peaksMaxPosPlus)] = -Inf
                peaksMaxPosPlus = max.col(peaksMaxPosPlus, ties.method = "first")

                # make new peaks centered arround max position
                currentPeaksPlus = rngToProcessPlus
                start(currentPeaksPlus) = start(currentPeaksPlus) + peaksMaxPosPlus -1
                end(currentPeaksPlus) = start(currentPeaksPlus)
                currentPeaksPlus = currentPeaksPlus + ((bsSize - 1) / 2)
                # store peaks
                rngCenterPlus = c(rngCenterPlus, currentPeaksPlus)
                # remove peak regions from rest of possible regions
                currentPeaksPlus = as(currentPeaksPlus + ((bsSize - 1) / 2), "GRangesList")

                # update peak regions that are left for processing
                rngToProcessPlus = unlist(psetdiff(rngToProcessPlus, currentPeaksPlus))
            }
            if (length(rngToProcessMinus) != 0) {
                peaksMaxPosMinus = as.matrix(sgnMergeMinus[rngToProcessMinus])
                peaksMaxPosMinus[is.na(peaksMaxPosMinus)] = -Inf
                peaksMaxPosMinus = max.col(peaksMaxPosMinus, ties.method = "last")

                currentPeaksMinus = rngToProcessMinus
                start(currentPeaksMinus) = start(currentPeaksMinus) + peaksMaxPosMinus -1
                end(currentPeaksMinus) = start(currentPeaksMinus)
                currentPeaksMinus = currentPeaksMinus + ((bsSize - 1) / 2)

                rngCenterMinus = c(rngCenterMinus, currentPeaksMinus)

                currentPeaksMinus = as(currentPeaksMinus + ((bsSize - 1) /
                                                                2), "GRangesList")

                rngToProcessMinus = unlist(psetdiff(rngToProcessMinus, currentPeaksMinus))
            }
            Counter = Counter + 1
        }
    }
    rngS4 = c(rngCenterPlus, rngCenterMinus)
    rngS4 = GenomeInfoDb::sortSeqlevels(rngS4)
    rngS4 = sort(rngS4)
    return(rngS4)
}

.filterMinCrosslinks <- function(rng,
                                 sgn,
                                 minCrosslinks) {
    sgnMergePlus = sgn$signalPlus
    sgnMergeMinus = sgn$signalMinus

    # check if both strands exists
    if ("+" %in% unique(strand(rng))) {
        # split by strand plus
        rngCurrPlus = rng[strand(rng) == "+"]
        rngCurrPlusMat = as.matrix(sgnMergePlus[rngCurrPlus])
        rngCurrPlus = rngCurrPlus[apply((rngCurrPlusMat > 0), 1, sum) > minCrosslinks]
    }
    if (!"+" %in% unique(strand(rng))) {
        rngCurrPlus = NULL
    }

    if ("-" %in% unique(strand(rng))) {
        # split by strand minus
        rngCurrMinus = rng[strand(rng) == "-"]
        rngCurrMinusMat = as.matrix(sgnMergeMinus[rngCurrMinus])
        rngCurrMinus = rngCurrMinus[apply((rngCurrMinusMat > 0), 1, sum) > minCrosslinks]
    }
    if (!"-" %in% unique(strand(rng))) {
        rngCurrMinus = NULL
    }

    # combine sort return
    rngCurr = c(rngCurrPlus, rngCurrMinus)
    rngCurr = GenomeInfoDb::sortSeqlevels(rngCurr)
    rngCurr = sort(rngCurr)
    return(rngCurr)
}

#' @importFrom S4Vectors queryHits
.filterCenterClSite <- function(rng, rng0) {
    rngCurr = rng - (unique(width(rng)) - 1) / 2
    rngCurr = rng[queryHits(findOverlaps(rngCurr, rng0))]
    # combine sort return
    rngCurr = GenomeInfoDb::sortSeqlevels(rngCurr)
    rngCurr = sort(rngCurr)
    return(rngCurr)
}

.filterCenterSummit <- function(rng, sgn) {
    sgnMergePlus = sgn$signalPlus
    sgnMergeMinus = sgn$signalMinus

    if ("+" %in% unique(strand(rng))) {
        rngPlus = rng[strand(rng) == "+"]
        rngCurrPlusMat = as.matrix(sgnMergePlus[rngPlus])
        rngCurrPlusCount = apply(rngCurrPlusMat, 1, max)
        rngCurrPlus = rngPlus[rngCurrPlusCount == rngCurrPlusMat[, ((unique(width(rng)) -1) / 2) + 1]]
    }
    if (!"+" %in% unique(strand(rng))) {
        rngCurrPlus = NULL
    }

    if ("-" %in% unique(strand(rng))) {
        rngMinus = rng[strand(rng) == "-"]
        rngCurrMinusMat = as.matrix(sgnMergeMinus[rngMinus])
        rngCurrMinusCount = apply(rngCurrMinusMat, 1, max)
        rngCurrMinus = rngMinus[rngCurrMinusCount == rngCurrMinusMat[, ((unique(width(rng)) - 1) / 2) + 1]]
    }
    if (!"-" %in% unique(strand(rng))) {
        rngCurrMinus = NULL
    }

    # combine sort return
    rngCurr = c(rngCurrPlus, rngCurrMinus)
    rngCurr = GenomeInfoDb::sortSeqlevels(rngCurr)
    rngCurr = sort(rngCurr)
    return(rngCurr)
}

.filterMinClSites <- function(rng, rng0, minClSites) {
    overlaps = findOverlaps(rng, rng0)
    freq = table(queryHits(overlaps))
    idx = as.numeric(names(freq[freq >= minClSites]))
    rngCurr = rng[idx]

    # combine sort return
    rngCurr = GenomeInfoDb::sortSeqlevels(rngCurr)
    rngCurr = sort(rngCurr)
    return(rngCurr)
}
