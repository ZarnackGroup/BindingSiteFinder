#' Replicate reproducibility filter function
#'
#' For each replicate the number of binding sites with a certain number of crosslinks is
#' calculated. A quantile based threshold (\code{cutoff}) is applied to each replicate.
#' This indicates how many of the merged binding sites are supported by crosslinks
#' from the respective replicate. Next, one can specify how many replicates need to
#' pass the defined threshold for a binding site to be considered reproducible.
#'
#' If \code{cutoff} is a single number then the indicated cutoff will be applied to
#' all replicates. If it is a vector then each element in the vector is applied to all
#' replicates of the respective condition. The order is hereby given by the levels
#' of the condition column of the meta data (see \code{\link{BSFDataSet}},
#' \code{\link{getMeta}}). If the condition specific filter is applied, a meta column
#' is added to the GRanges of the \code{BSFDataSet} object, indicating the support for each condition.
#'
#' If \code{n.reps} is a single number then this number is used as treshold for all
#' binding sites. If it is a vector then it is applied to the replicates of the respective
#' condition (like in \code{cutoff}). This allows the application of different thresholds
#' for experiments of different experimental conditions. If the condition specific filter
#' is applied, a meta column is added to the GRanges of the \code{BSFDataSet} object,
#' indicating the support for each condition.
#'
#' @param object a BSFDataSet object
#' @param cutoff a vector of length = 1, or of length = levels(getMeta(object)$conditions) with
#' a single number (between 0-1) indicating the quantile cutoff
#' @param n.reps a vector of length = 1, or of length = levels(getMeta(object)$conditions) indicating
#' how many replicates need to meet their threshold for a binding site to be called reproducible.
#' @param min.crosslinks numeric of length = 1, defines the lower boundary for the minimum
#' number of crosslinks a binding site has to be supported by all replicates, regardless
#' of the replicate specific quantile threshold
#' @param returnType one of "BSFDataSet" or "data.frame". "BSFDataSet" is the default and
#' "matrix" can be used for easy plotting.
#'
#' @return an object of type BSFDataSet
#'
#' @import tidyr GenomicRanges
#'
#' @examples
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
#' # merge binding sites
#' bds <- makeBindingSites(object = bds, bsSize = 9, minWidth = 2,
#' minCrosslinks = 2, minClSites = 1)
#'
#' # use default return with single threshold
#' s = reproducibilityFilter(bds, cutoff = c(0.05), n.reps = c(3))
#'
#' # use default return with condition specific threshold
#' s = reproducibilityFilter(bds, cutoff = c(0.1, 0.05), n.reps = c(1, 2))
#'
#' # use data.frame return type for plotting
#' s = reproducibilityFilter(bds, cutoff = c(0.1, 0.05), n.reps = c(1, 2),
#' returnType = "data.frame")
#' library(ComplexHeatmap)
#' m = make_comb_mat(s)
#' UpSet(m)
#'
#' @export
reproducibilityFilter <- function(object,
                                  cutoff = 0.05,
                                  n.reps = 2,
                                  min.crosslinks = 1,
                                  returnType = c("BSFDataSet", "data.frame")) {
    stopifnot(is(object, "BSFDataSet"))

    cond = getMeta(object)$condition
    df = coverageOverRanges(object, returnType = "data.frame")

    if (length(cutoff) == 1) {
        # calculate sample specific thresholds
        qSel = .selectQuantilesSingleCondtion(
            covDf = df,
            userCond = cond,
            userNreps = n.reps,
            userCutoff = cutoff
        )
        # apply minimal crosslink threshold
        qSel$value = ifelse(qSel$value < min.crosslinks,
                            qSel$value + min.crosslinks,
                            qSel$value)
        matchIdx = match(qSel$name, colnames(df))

        # calculate replicate support based on quantile cutoff
        s = apply(df, 1, function(x) {
            ifelse(x > qSel$value[matchIdx], 1, 0)
        }) %>%
            t %>% as.data.frame() #TODO use this for plotting the upset plot
        support = rowSums(s) > n.reps
        # store results in output
        newRanges = getRanges(object)
        newRanges = newRanges[support]
        newObject = setRanges(object, newRanges)
    }

    if (length(cutoff) > 1) {
        if (length(levels(cond)) == 1) {
            stop("multiple cutoffs are given but only one condition exists")
        }
        # calculate sample specific thresholds
        qSel = .selectQuantilesMultipleConditions(
            covDf = df,
            userCond = cond,
            userNreps = n.reps,
            userCutoff = cutoff
        )
        # apply minimal crosslink threshold
        qSel$value = ifelse(qSel$value < min.crosslinks,
                            qSel$value + min.crosslinks,
                            qSel$value)
        matchIdx = match(qSel$name, colnames(df))

        # calculate reproducibility per condition
        s = apply(df, 1, function(x) {
            ifelse(x > qSel$value[matchIdx], 1, 0)
        }) %>%
            t %>% as.data.frame() #TODO use this for plotting the upset plot
        sSplit = sapply(levels(cond), function(x) {
            s %>% dplyr::select(contains(x)) %>% rowSums()
        }) %>% as.data.frame()
        idx = match(colnames(sSplit), qSel$applyTo)
        support = apply(sSplit, 1, function(x) {
            x >= qSel$n.reps[idx]
        }) %>%
            t %>% as.data.frame()

        supportAll = apply(support, 1, any)

        newRanges = getRanges(object)
        mcols(newRanges) = support
        newRanges = newRanges[supportAll]
        newObject = setRanges(object, newRanges)
    }

    returnType = match.arg(returnType, choices = c("BSFDataSet", "data.frame"))

    if (returnType == "BSFDataSet") {
        retObj = newObject
    }
    if (returnType == "data.frame") {
        retObj = s
    }

    return(retObj)
}


#' Coverage function for BSFDataSet objects
#'
#' The crosslink coverage is computed for all ranges in the the given \code{BSFDataSet}
#' object (see \code{\link{BSFDataSet}} for details). The coverage of the crosslinks
#' stored in the singnal slot of the \code{BSFDataSet} is computed over the ranges
#' in the \code{BSFDataSet} object.
#'
#' If \code{returnType} is set to GRanges, the coverge information is stored in the
#' metadata slot.
#'
#' @param object a BSFDataSet object
#' @param merge logical, if coverage should be merged per replicate or reported
#' for each nucleotide in the range individually
#' @param returnType one of "GRanges", "matrix" or "data.frame"
#'
#' @return an object of class specified in \code{returnType}
#' @import GenomicRanges
#'
#' @examples
#' # load data
#' csFile <- system.file("extdata", "PureCLIP_crosslink_sites_example.bed",
#'  package="BindingSiteFinder")
#' cs = rtracklayer::import(con = csFile, format = "BED")
#' clipFiles <- system.file("extdata", package="BindingSiteFinder")
#'
#' # one experimental condition
#' meta = data.frame(condition = c("WT", "WT", "WT", "WT"),
#' clPlus = list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE),
#' clMinus = list.files(clipFiles, pattern = "minus.bw$", full.names = TRUE))
#' bds = BSFDataSet(ranges = cs, meta = meta)
#'
#' rng = coverageOverRanges(bds)
#'
#' @export
coverageOverRanges <- function(object,
                               merge = TRUE,
                               returnType = c("GRanges", "matrix", "data.frame")) {
    stopifnot(is(object, "BSFDataSet"))
    # split by strand
    rng = getRanges(object)
    rngPlus = rng[strand(rng) == "+"]
    rngMinus = rng[strand(rng) == "-"]
    # prepare signal
    sgn = getSignal(object)
    # signal coverage is reported for each position in the range of the peak
    if (!isTRUE(merge)) {
        # manage return type
        # only return type data.frame is possible with this option
        returnType = match.arg(returnType,
                               choices = c("GRanges", "matrix", "data.frame"))
        if (returnType != "data.frame") {
            warning("Only return type 'data.frame' possible with non-merged output.")
        }
        returnType = "data.frame"

        if (length(rngPlus) > 0) {
            matPlus = lapply(sgn$signalPlus, function(x) {
                as.matrix(x[rngPlus])
            })
            covPlus = do.call(rbind, lapply(matPlus, colSums))
        }
        if (length(rngPlus) == 0) {
            covPlus = 0
        }
        if (length(rngMinus) > 0) {
            matMinus = lapply(sgn$signalMinus, function(x) {
                as.matrix(x[rngMinus])
            })
            covMinus = do.call(rbind, lapply(matMinus, colSums))
            # flip orientation of minus strand coverage
            covMinus = covMinus %>% as.data.frame() %>% rev() %>% as.matrix()
        }
        if (length(rngMinus) == 0) {
            covMinus = 0
        }
        covDf = covPlus + covMinus
        retObj = as.data.frame(covDf)
    }
    # signal is merged over all positions in the range
    if (isTRUE(merge)) {
        mcols(rngPlus) = as.matrix(
            do.call(cbind, lapply(sgn$signalPlus, function(x) {
                sum(x[rngPlus])
            })))
        mcols(rngMinus) = as.matrix(
            do.call(cbind, lapply(sgn$signalMinus, function(x) {
                sum(x[rngMinus])
            })))
        # sort ranges
        rngCov = c(rngPlus, rngMinus)
        rngCov = GenomeInfoDb::sortSeqlevels(rngCov)
        rngCov = sort(rngCov)
        # manage return type
        returnType = match.arg(returnType,
                               choices = c("GRanges", "matrix", "data.frame"))
        if (returnType == "GRanges") {
            retObj = rngCov
        }
        if (returnType == "matrix") {
            retObj = as.matrix(mcols(rngCov))
        }
        if (returnType == "data.frame") {
            retObj = as.data.frame(mcols(rngCov))
        }
    }
    return(retObj)
}


#' Annotation function for BSFDataSet object
#'
#' This function can be used to annotate a \code{BSFDataSet} object with merged binding
#' sites with scores from the initial ranges (eg. PureCLIP scores).
#'
#' @param object a BSFDataSet object
#' @param scoreRanges a GRanges object, with numeric column named 'score'
#'
#' @return an object of class BSFDataSet with updated meta columns of the ranges
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits
#'
#' @examples
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
#' # merge binding sites
#' bds <- makeBindingSites(object = bds, bsSize = 9, minWidth = 2,
#' minCrosslinks = 2, minClSites = 1)
#'
#' # annotate with original pureCLIP score
#' bdsRe = annotateWithScore(bds, cs)
#'
#' @export
annotateWithScore <- function(object,
                              scoreRanges) {

    # bind locally used variables
    qHits <- NULL

    stopifnot(is(object, "BSFDataSet"))

    msg = tryCatch(getSummary(object), error=function(e) e, warning=function(w) w)
    if (is(msg, "warning")) {
        stop("Function makeBindingSites() was not run.")
    }

    stopifnot(is(scoreRanges, "GRanges"))

    if (!any(colnames(mcols(scoreRanges)) == "score")) {
        stop("Ranges do not have a meta column named score.")
    }
    if (!is.numeric(scoreRanges$score)) {
        stop("Score column must have numeric values.")
    }

    rng = getRanges(object)
    ol = findOverlaps(rng, scoreRanges)

    if (length(ol) == 0) {
        stop("Ranges of 'object' and 'scoreRanges' do not overlap.")
    }

    matchDF = data.frame(qHits = queryHits(ol),
                         sHits = subjectHits(ol),
                         score = scoreRanges$score[subjectHits(ol)])

    scores = dplyr::group_by(matchDF, qHits) %>%
        dplyr::summarize(pSum = sum(score),
                         pMax = max(score),
                         pMean = mean(score), .groups = "drop") %>%
        as.data.frame()

    mcols(rng)$scoreSum = scores$pSum
    mcols(rng)$scoreMean = scores$pMean
    mcols(rng)$scoreMax = scores$pMax

    newObject = setRanges(object, rng)
    return(newObject)
}




#' Support ratio function for BSFDataSet objects
#'
#' Functions that computes a ratio to determine how well a given binding site with
#' is supported by the crosslink coverage of the data. For a given \code{BSFDataSet}
#' object binding sites are computed for each width indicated in the \code{bsWidths}
#' vector (using the \code{\link{coverageOverRanges}} function). These coverages
#' are compared to the coverage of regions flanking the binding sites. If not
#' indicated in \code{bsFlank} these regions are of the same width as the binding
#' sites.
#'
#' Testing the width of 3nt for example, would result in a coverage within all
#' 3nt wide binding sites (c1) and a coverage computed on the adjacent 3nt
#' flanking the binding sites up- and downstream (f1, f2). Based on these numbers
#' the ratio is computed by: c1/(1/2(f1+f2)).
#'
#' The median over all ratios is reported as representative value.
#'
#' @param object a BSFDataSet object
#' @param bsWidths a numeric vector indicating the different binding site
#' width to compute the ratio for
#' @param bsFlank optional; a numeric vector of the same length as \code{bsWidth}
#' used to specify the width of the flanking regions
#' @param ... further arguments passed to \code{makeBindingSites}
#'
#' @return an object of class \code{data.frame}
#'
#' @examples
#' # load data
#' csFile <- system.file("extdata", "PureCLIP_crosslink_sites_example.bed",
#'  package="BindingSiteFinder")
#' cs = rtracklayer::import(con = csFile, format = "BED")
#' clipFiles <- system.file("extdata", package="BindingSiteFinder")
#'
#' # one experimental condition
#' meta = data.frame(condition = c("WT", "WT", "WT", "WT"),
#' clPlus = list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE),
#' clMinus = list.files(clipFiles, pattern = "minus.bw$", full.names = TRUE))
#' bds = BSFDataSet(ranges = cs, meta = meta)
#'
#' supportRatio(bds, bsWidths = c(3,7))
#'
#' @export
supportRatio <- function(object, bsWidths, bsFlank = NA, ...) {
    stopifnot(is(object, "BSFDataSet"))
    if (!any(is.numeric(bsWidths))) {
        stop("bsWidth needs to be numeric. ")
    }
    if (any(round(bsWidths) != bsWidths)) {
        stop("bsWidth is not an integer. ")
    }
    # check bsFlank vector if set
    if (!any(is.na(bsFlank))) {
        if (!any(is.numeric(bsFlank))) {
            stop("bsFlank needs to be numeric. ")
        }
        if (any(round(bsFlank) != bsFlank)) {
            stop("bsFlank is not an integer. ")
        }
        if (length(bsWidths) != length(bsFlank)) {
            stop("bsWidths and bsFlank needs to be vectors of the same length. ")
        }
    }
    # use same width for flanking regions as the binding sites are
    if (any(is.na(bsFlank))) {
        bsFlank = bsWidths
    }
    # calculate ratio
    objList = lapply(bsWidths, function(x){
        makeBindingSites(object = object, bsSize = x, ...)
    })
    objScore = lapply(seq_along(bsWidths), function(x){
        .computeSupportRatio(object = objList[[x]], flankSize = bsFlank[x])
    })
    # return results
    resDf = data.frame(bsWidths = factor(bsWidths),
                       supportRatio = unlist(objScore))
    return(resDf)
}


