#' Replicate reproducibility filter function
#'
#' For each replicate the number of binding sites with a certain number of
#' crosslinks is calculated. A quantile based threshold (\code{cutoff}) is
#' applied to each replicate. This indicates how many of the merged binding
#' sites are supported by crosslinks from the respective replicate. Next, one
#' can specify how many replicates need to pass the defined threshold for a
#' binding site to be considered reproducible.
#'
#' If \code{cutoff} is a single number then the indicated cutoff will be
#' applied to all replicates. If it is a vector then each element in the vector
#' is applied to all replicates of the respective condition. The order is
#' hereby given by the levels of the condition column of the meta data
#' (see \code{\link{BSFDataSet}},\code{\link{getMeta}}). If the condition
#' specific filter is applied, a meta column is added to the GRanges of the
#' \code{BSFDataSet} object, indicating the support for each condition.
#'
#' If \code{n.reps} is a single number then this number is used as treshold for
#' all binding sites. If it is a vector then it is applied to the replicates of
#' the respective condition (like in \code{cutoff}). This allows the
#' application of different thresholds for experiments of different
#' experimental conditions. If the condition specific filter is applied, a meta
#' column is added to the GRanges of the \code{BSFDataSet} object,
#' indicating the support for each condition.
#'
#' @param object a BSFDataSet object
#' @param cutoff a vector of length = 1, or of length =
#' levels(getMeta(object)$conditions) with a single number (between 0-1)
#' indicating the quantile cutoff
#' @param n.reps a vector of length = 1, or of length = l
#' evels(getMeta(object)$conditions) indicating how many replicates need to
#' meet their threshold for a binding site to be called reproducible.
#' @param min.crosslinks numeric of length = 1, defines the lower boundary for
#' the minimum number of crosslinks a binding site has to be supported by all
#' replicates, regardless of the replicate specific quantile threshold
#' @param returnType one of "BSFDataSet" or "data.frame". "BSFDataSet" is the
#' default and "matrix" can be used for easy plotting.
#'
#' @return an object of type BSFDataSet
#'
#' @import tidyr GenomicRanges
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' # merge binding sites
#' bds <- makeBindingSites(object = bds, bsSize = 9, minWidth = 2,
#' minCrosslinks = 2, minClSites = 1)
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
                                  cutoff,
                                  n.reps,
                                  min.crosslinks = 1,
                                  returnType = c("BSFDataSet", "data.frame")) {
    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))

    if (length(cutoff) != length(n.reps)) {
        stop("Number of values for 'cutoff' does not match the number of values
             for 'n.reps'. ")
    }

    metaData = getMeta(object)
    numberOfConditions = length(levels(metaData$condition))
    if (numberOfConditions != length(reproducibilityCutoff) && !is.null(reproducibilityCutoff)) {
        msg = paste0("Reproducibility filter cutoff does not match the number of conditions. You specified: ",
                     length(reproducibilityCutoff), ", but there are: ", numberOfConditions, "\n")
        msg2 = paste0("The specified cutoff (",  reproducibilityCutoff, ") ",
                      "is applied to all conditions (",
                      paste(as.character(levels(metaData$condition)), collapse = ",") ,") \n")
        warning(paste0(msg, msg2))
        defaultReproducibilityCutoff = rep(reproducibilityCutoff, numberOfConditions)
    }
    if (is.null(reproducibilityCutoff)) {
        msg = paste0("Reproducibility cutoff not defined. Defaults to 0.05 for each condition. \n")
        message(msg)
        defaultReproducibilityCutoff = rep(0.05, numberOfConditions)
    } else {
        defaultReproducibilityCutoff = reproducibilityCutoff
    }
    # Manage parameter n.reps
    if (numberOfConditions != length(n.reps) && !is.null(n.reps)) {
        msg = paste0("Parameter n.reps does not match the number of conditions. You specified: ",
                     length(n.reps), ", but there are: ", numberOfConditions, "\n")
        msg2 = paste0("n.reps defaults to N-1 for each condition. \n")
        warning(paste0(msg, msg2))

        n.conditions = table(metaData$condition) %>% as.data.frame()
        defaultNreps = n.conditions$Freq -1
    }
    if (is.null(n.reps)) {
        msg = paste0("Parameter n.reps not defined. Defaults to N-1 for each condition. \n")
        message(paste0(msg))
        n.conditions = table(metaData$condition) %>% as.data.frame()
        defaultNreps = n.conditions$Freq -1
    } else {
        defaultNreps = n.reps
    }

    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    cond = getMeta(object)$condition
    df = as.data.frame(mcols(coverageOverRanges(
        object, returnOptions = "merge_positions_keep_replicates",
        silent = TRUE)))

    # Manage single cutoff for single condition
    if (length(cutoff) == 1) {
        if(length(levels(cond)) > 1) {
            stop("Only one cutoff is given for multiple conditions.")
        }
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
        support = rowSums(s) >= n.reps
        # store results in output
        newRanges = getRanges(object)
        newRanges = newRanges[support]
        newObject = setRanges(object, newRanges)
    }

    # Manage multiple cutoffs for multipe conditions
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
        sSplit = vapply(levels(cond), function(x) {
            s %>% dplyr::select(contains(x)) %>% rowSums()
        }, FUN.VALUE = numeric(nrow(s))) %>% as.data.frame()
        # idx = match(colnames(sSplit), qSel$applyTo)
        idx = match(colnames(sSplit), qSel$sel)
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

    # Manage return options
    returnType = match.arg(returnType, choices = c("BSFDataSet", "data.frame"))
    if (returnType == "BSFDataSet") {
        retObj = newObject
    }
    if (returnType == "data.frame") {
        retObj = s
    }

    return(retObj)
}


#' Annotation function for BSFDataSet object
#'
#' This function can be used to annotate a \code{BSFDataSet} object with
#' merged binding sites with scores from the initial ranges
#' (eg. PureCLIP scores).
#'
#' @param object a BSFDataSet object
#' @param scoreRanges a GRanges object, with numeric column named 'score'
#'
#' @return an object of class BSFDataSet with updated meta columns of the ranges
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits
#'
#' @examples
#' if (.Platform$OS.type != "windows") {
#'     # load data
#'     csFile <- system.file("extdata", "PureCLIP_crosslink_sites_examples.bed",
#'                         package="BindingSiteFinder")
#'     cs = rtracklayer::import(con = csFile, format = "BED",
#'     extraCols=c("additionalScores" = "character"))
#'     cs$additionalScores = NULL
#'     clipFiles <- system.file("extdata", package="BindingSiteFinder")
#'     # two experimental conditions
#'     meta = data.frame(
#'     id = c(1,2,3,4),
#'     condition = factor(c("WT", "WT", "KD", "KD"),
#'     levels = c("KD", "WT")),
#'     clPlus = list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE),
#'     clMinus = list.files(clipFiles, pattern = "minus.bw$",
#'      full.names = TRUE))
#'     bds = BSFDataSetFromBigWig(ranges = cs, meta = meta, silent = TRUE)
#'
#'     # merge binding sites
#'     bds <- makeBindingSites(object = bds, bsSize = 9, minWidth = 2,
#'     minCrosslinks = 2, minClSites = 1)
#'
#'     # annotate with original pureCLIP score
#'     bdsRe = annotateWithScore(bds, cs)
#' }
#' @export
annotateWithScore <- function(object,
                              scoreRanges) {

    # bind locally used variables
    qHits <- NULL

    stopifnot(is(object, "BSFDataSet"))

    msg = tryCatch(getSummary(object), error=function(e) e,
                   warning=function(w) w)
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
#' Functions that computes a ratio to determine how well a given binding site
#' with is supported by the crosslink coverage of the data. For a given
#' \code{BSFDataSet} object binding sites are computed for each width indicated
#' in the \code{bsWidths} vector (using the \code{\link{coverageOverRanges}}
#' function). These coverages are compared to the coverage of regions flanking
#' the binding sites. If not indicated in \code{bsFlank} these regions are of
#' the same width as the binding sites.
#'
#' Testing the width of 3nt for example, would result in a coverage within all
#' 3nt wide binding sites (c1) and a coverage computed on the adjacent 3nt
#' flanking the binding sites up- and downstream (f1, f2). Based on these
#'  numbers the ratio is computed by: c1/(1/2(f1+f2)).
#'
#' The median over all ratios is reported as representative value.
#'
#' @param object a BSFDataSet object
#' @param bsWidths a numeric vector indicating the different binding site
#' width to compute the ratio for
#' @param bsFlank optional; a numeric vector of the same length as
#' \code{bsWidth} used to specify the width of the flanking regions
#' @param sub.chr chromosome identifier (eg, chr1, chr2) used for subsetting the
#' BSFDataSet object.
#' @param approximate logical; if binding sites should be approxiamted by their
#' center position instead of running \code{makeBindingSites} on every iteration
#' @param ... further arguments passed to \code{makeBindingSites}
#'
#' @return an object of class \code{data.frame}
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' supportRatio(bds, bsWidths = c(3,7))
#'
#' @export
supportRatio <- function(object, bsWidths, bsFlank = NA, sub.chr, approximate = FALSE, ...) {
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
            stop("bsWidths and bsFlank needs to be vectors
                 of the same length. ")
        }
    }
    # use same width for flanking regions as the binding sites are
    if (any(is.na(bsFlank))) {
        bsFlank = bsWidths
    }
    # compute binding sites
    print("make bs")
    if (isTRUE(approximate)) {
        objList = lapply(bsWidths, function(x){
            .approximateBindingSites(object = object, bsSize = x, sub.chr = sub.chr)
        })
    } else {
        objList = lapply(bsWidths, function(x){
            makeBindingSites(object = object, bsSize = x, sub.chr = sub.chr, ...)
        })
    }
    # calculate ratio
    print("calc ratio")
    objScore = lapply(seq_along(bsWidths), function(x){
        .computeSupportRatio(object = objList[[x]], flankSize = bsFlank[x])
    })
    # return results
    resDf = data.frame(bsWidths = factor(bsWidths),
                       supportRatio = unlist(objScore))
    return(resDf)
}


