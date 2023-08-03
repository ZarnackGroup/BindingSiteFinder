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
#' If \code{nReps} is a single number then this number is used as treshold for
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
#' @param nReps a vector of length = 1, or of length = l
#' evels(getMeta(object)$conditions) indicating how many replicates need to
#' meet their threshold for a binding site to be called reproducible.
#' @param minCrosslinks numeric of length = 1, defines the lower boundary for
#' the minimum number of crosslinks a binding site has to be supported by all
#' replicates, regardless of the replicate specific quantile threshold
#' @param returnType one of "BSFDataSet" or "data.frame". "BSFDataSet" is the
#' default and "matrix" can be used for easy plotting.
#' @param n.reps deprecated -> use nReps instead
#' @param min.crosslinks deprecated -> use minCrosslinks instead
#' @param quiet logical; whether to print messages
#'
#' @return an object of type BSFDataSet
#'
#' @import tidyr GenomicRanges lifecycle
#' @importFrom dplyr count select
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' # merge binding sites
#' bds <- makeBindingSites(object = bds, bsSize = 9)
#'
#' # use default return with condition specific threshold
#' bds = reproducibilityFilter(bds, cutoff = 0.1, nReps = 1)
#'
#' @export
reproducibilityFilter <- function(object,
                                  cutoff = NULL,
                                  nReps = NULL,
                                  minCrosslinks = 1,
                                  returnType = c("BSFDataSet", "data.frame"),
                                  n.reps = lifecycle::deprecated(),
                                  min.crosslinks = lifecycle::deprecated(),
                                  quiet = FALSE) {
    # initialize local vairables
    name <- value <- NULL

    # Argument deprecation warnings
    # --------------------------------------------------------------------------

    if (lifecycle::is_present(n.reps)) {
        lifecycle::deprecate_warn("2.0.0", "reproducibilityFilter(n.reps = )", "reproducibilityFilter(nReps = )")
        nReps = n.reps
    }

    if (lifecycle::is_present(min.crosslinks)) {
        lifecycle::deprecate_warn("2.0.0", "reproducibilityFilter(min.crosslinks = )", "reproducibilityFilter(minCrosslinks = )")
        minCrosslinks = min.crosslinks
    }

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is.logical(quiet))

    if (length(cutoff) != length(nReps)) {
        stop("Number of values for 'cutoff' does not match the number of values
             for 'nReps'. ")
    }

    metaData = getMeta(object)

    numberOfConditions = length(levels(metaData$condition))
    if (numberOfConditions != length(cutoff) && !is.null(cutoff)) {
        msg = paste0("Reproducibility filter cutoff does not match the number of conditions. You specified: ",
                     length(cutoff), ", but there are: ", numberOfConditions, "\n")
        msg2 = paste0("The specified cutoff (",  cutoff, ") ",
                      "is applied to all conditions (",
                      paste(as.character(levels(metaData$condition)), collapse = ",") ,") \n")
        warning(paste0(msg, msg2))
        defaultCutoff = rep(cutoff, numberOfConditions)
    }
    if (is.null(cutoff)) {
        msg = paste0("Reproducibility cutoff not defined. Defaults to 0.05 for each condition. \n")
        if(!quiet) message(msg)
        defaultCutoff = rep(0.05, numberOfConditions)
    }
    if (numberOfConditions == length(cutoff) && !is.null(cutoff)) {
        defaultCutoff = cutoff
    }
    # Manage parameter nReps
    if (numberOfConditions != length(nReps) && !is.null(nReps)) {
        msg = paste0("Parameter nReps does not match the number of conditions. You specified: ",
                     length(nReps), ", but there are: ", numberOfConditions, "\n")
        msg2 = paste0("nReps defaults to N-1 for each condition. \n")
        warning(paste0(msg, msg2))

        n.conditions = table(metaData$condition) %>% as.data.frame()
        defaultNreps = n.conditions$Freq -1
    }
    if (is.null(nReps)) {
        msg = paste0("Parameter nReps not defined. Defaults to N-1 for each condition. \n")
        if(!quiet) message(paste0(msg))
        n.conditions = table(metaData$condition) %>% as.data.frame()
        defaultNreps = n.conditions$Freq -1
    }
    if (numberOfConditions == length(nReps) && !is.null(nReps)) {
        defaultNreps = nReps
    }

    # ---
    # Store function parameters in list
    optstr = list(cutoff = defaultCutoff, nReps = defaultNreps, minCrosslinks = minCrosslinks)
    object@params$reproducibilityFilter = optstr

    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    rngInitial = getRanges(object)
    cond = metaData$condition
    # get number of crosslinks per binding site and replicate
    df = as.data.frame(mcols(coverageOverRanges(
        object, returnOptions = "merge_positions_keep_replicates",
        silent = TRUE)))

    # Manage single cutoff for single condition
    if (length(defaultCutoff) == 1) {
        if(length(levels(cond)) > 1) {
            stop("Only one cutoff is given for multiple conditions.")
        }
        # calculate sample specific thresholds
        qSel = .selectQuantilesSingleCondtion(
            covDf = df,
            userCond = cond,
            userNreps = defaultNreps,
            userCutoff = defaultCutoff
        )
        # apply minimal crosslink threshold
        qSel$value = ifelse(qSel$value < minCrosslinks,
                            minCrosslinks,
                            qSel$value)
        matchIdx = match(qSel$name, colnames(df))

        # ---
        # Store data for diagnostic plot in list
        dfPlot = df %>%
            pivot_longer(everything()) %>%
            group_by(name, value) %>% dplyr::count() %>%
            separate(name, into = c(NA, "condition"), sep = "_", remove = FALSE)
        object@plotData$reproducibilityFilterPlot$data = dfPlot
        object@plotData$reproducibilityFilterPlot$cutoffs = qSel

        # calculate replicate support based on quantile cutoff
        s = apply(df, 1, function(x) {
            ifelse(x > qSel$value[matchIdx], 1, 0)
        }) %>%
            t() %>% as.data.frame()

        # ---
        # Store results for plotting
        object@plotData$reproducibilitySamplesPlot$data = s

        # Filter binding sites ranges by replicate support
        support = rowSums(s) >= defaultNreps
        newRanges = getRanges(object)
        newRanges = newRanges[support]
        newObject = setRanges(object, newRanges)
    }

    # Manage multiple cutoffs for multiple conditions
    if (length(defaultCutoff) > 1) {
        if (length(levels(cond)) == 1) {
            stop("multiple cutoffs are given but only one condition exists")
        }
        # calculate sample specific thresholds
        qSel = .selectQuantilesMultipleConditions(
            covDf = df,
            userCond = cond,
            userNreps = defaultNreps,
            userCutoff = defaultCutoff
        )
        # apply minimal crosslink threshold
        qSel$value = ifelse(qSel$value < minCrosslinks,
                            minCrosslinks,
                            qSel$value)
        matchIdx = match(qSel$name, colnames(df))

        # ---
        # Store data for diagnostic plot in list
        dfPlot = df %>%
            pivot_longer(everything()) %>%
            group_by(name, value) %>% dplyr::count() %>%
            separate(name, into = c(NA, "condition"), sep = "_", remove = FALSE)
        object@plotData$reproducibilityFilterPlot$data = dfPlot
        object@plotData$reproducibilityFilterPlot$cutoffs = qSel

        # calculate reproducibility per condition
        s = apply(df, 1, function(x) {
            ifelse(x > qSel$value[matchIdx], 1, 0)
        }) %>%
            t %>% as.data.frame()

        # ---
        # Store results for plotting
        object@plotData$reproducibilitySamplesPlot$data = s

        # Filter binding sites ranges by replicate support
        sSplit = vapply(levels(cond), function(x) {
            s %>% dplyr::select(contains(x)) %>% rowSums()
        }, FUN.VALUE = numeric(nrow(s))) %>% as.data.frame()
        idx = match(colnames(sSplit), qSel$sel)
        support = apply(sSplit, 1, function(x) {
            x >= qSel$defaultNreps[idx]
        }) %>%
            t %>% as.data.frame()

        supportAll = apply(support, 1, any)

        newRanges = getRanges(object)
        # mcols(newRanges) = support
        mcols(newRanges) = cbind.data.frame(mcols(newRanges), support)
        newRanges = newRanges[supportAll]
        newObject = setRanges(object, newRanges, quiet = quiet)
    }

    # ---
    # Store for results
    resultLine = data.frame(
        funName = "reproducibilityFilter()", class = "binding sites",
        nIn = length(rngInitial), nOut = length(newRanges),
        per = paste0(round(length(newRanges)/ length(rngInitial), digits = 2)*100,"%"),
        options = paste0("Cutoff=", paste(optstr$cutoff, collapse = ", "),
                         ", nReps=", paste(optstr$nReps, collapse = ", "),
                         ", minCrosslinks=", paste(optstr$minCrosslinks, collapse = ", "))
    )
    newObject@results = rbind(newObject@results, resultLine)

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
#' @param match.ranges a GRanges object, with numeric column for the score to match
#' @param match.score character; meta column name of the crosslink site
#' \code{\link{GenomicRanges}} object that holds the score to match
#' @param match.option character; option how score should be matched
#' @param quiet logical; whether to print messages
#' @param scoreRanges deprecated -> use match.ranges instead
#' @param MatchColScore deprecated -> use match.score instead
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
annotateWithScore <- function(object, # bindingSiteFinder
                              match.ranges = NULL,
                              match.score = "score",
                              match.option = c("max", "sum", "mean"),
                              scoreRanges = lifecycle::deprecated(),
                              MatchColScore = lifecycle::deprecated(),
                              quiet = FALSE
) {

    # bind locally used variables
    qHits <- NULL

    # Argument deprecation warnings
    # --------------------------------------------------------------------------

    if (lifecycle::is_present(scoreRanges)) {
        lifecycle::deprecate_warn("2.0.0", "annotateWithScore(scoreRanges = )", "annotateWithScore(match.ranges = )")
        match.ranges = scoreRanges
    }

    if (lifecycle::is_present(MatchColScore)) {
        lifecycle::deprecate_warn("2.0.0", "annotateWithScore(MatchColScore = )", "annotateWithScore(match.score = )")
        match.score = MatchColScore
    }

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is.logical(quiet))

    # check input ranges
    if(is.null(match.ranges)) {
        msg = paste0("No ranges to match the score from. Please provide a GRanges object with a matching score column in the meta data.")
        stop(msg)
    }
    stopifnot(is(match.ranges, "GRanges"))

    # handle options
    match.option = match.arg(match.option, choices = c("max", "sum", "mean"))

    # match score column
    scoreColNames = colnames(mcols(match.ranges))
    if (!(match.score %in% scoreColNames)) {
        msg = paste0("Matching columns (", match.score,
                     ") is not present in the provided match.ranges. \n")
        stop(msg)
    }

    # ---
    # Store function parameters in list
    optstr = list(match.score = match.score, match.option = match.option)
    object@params$annotateWithScore = optstr

    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    rngInitial = getRanges(object)
    rng = getRanges(object)
    ol = findOverlaps(rng, match.ranges)

    if (length(ol) == 0) {
        stop("Ranges of 'object' and 'match.ranges' do not overlap.")
    }

    matchDF = data.frame(qHits = queryHits(ol),
                         sHits = subjectHits(ol),
                         score = match.ranges$score[subjectHits(ol)])

    if (match.option == "max") {
        score = dplyr::group_by(matchDF, qHits) %>%
            dplyr::summarize(score = max(score), .groups = "drop") %>%
            as.data.frame()
    }
    if (match.option == "sum") {
        score = dplyr::group_by(matchDF, qHits) %>%
            dplyr::summarize(score = sum(score), .groups = "drop") %>%
            as.data.frame()
    }
    if (match.option == "mean") {
        score = dplyr::group_by(matchDF, qHits) %>%
            dplyr::summarize(score = mean(score), .groups = "drop") %>%
            as.data.frame()
    }

    # ---
    # Store results for plotting
    object@plotData$annotateWithScore$data = score$score

    # ---
    # Store for results
    resultLine = data.frame(
        funName = "annotateWithScore()", class = "binding sites",
        nIn = length(rngInitial), nOut = length(rng),
        per = paste0(round(length(rng)/ length(rngInitial), digits = 2)*100,"%"),
        options = paste0("MatchOption=", optstr$match.option, ", match.score=", match.score)
    )
    object@results = rbind(object@results, resultLine)

    mcols(rng)$score = score$score
    newObject = setRanges(object, rng, quiet = quiet)
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
#' @param ... further arguments passed to \code{makeBindingSites}
#'
#' @return an object of class \code{data.frame}
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' suppressWarnings(supportRatio(bds, bsWidths = c(3,7)))
#'
#' @export
supportRatio <- function(object, bsWidths, bsFlank = NA, sub.chr = NA, ...) {
    # deprecation notice
    .Deprecated("estimateBsWidth")

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
    objList = lapply(bsWidths, function(x){
        makeBindingSites(object = object, bsSize = x, sub.chr = sub.chr, ...)
    })
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



#' Calculate signal-to-flank score
#'
#' This function calculates the signal-to-flank ratio for all present binding
#' sites.
#'
#' Each input range is treated as a binding site. For a particular binding site
#' all overlapping crosslinks are summed up and divided by the normalized sum of
#' the crosslinks in the two adjecent regions of the same size. This is done
#' for all bining sites and the ratio is reported as a score.
#'
#' @param object a BSFDataSet object
#' @param flank character; how the flanking region shoule be set. Options are
#' 'bs', 'manual'
#' @param flank.size numeric; if flank='manual' provide the desired flanking size
#' @param quiet logical; whether to print messages
#'
#' @return an object of class \code{\link{BSFDataSet}} with signal-to-flank ratios
#' added to the meta column of the ranges
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' bds = makeBindingSites(bds, bsSize = 5)
#' bds = calculateSignalToFlankScore(bds)
#'
#' @export
calculateSignalToFlankScore <- function(
        object,
        flank = c("bs", "manual"),
        flank.size = NULL,
        quiet = FALSE
){
    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is.logical(quiet))

    flank = match.arg(flank, choices = c("bs", "manual"))

    # ---
    # Store function parameters in list
    optstr = list(flank = flank)
    object@params$calculateSignalToFlankScore = optstr

    # PREPARE TEST RANGES + SIGNAL
    # --------------------------------------------------------------------------
    # collapse signal from replicates
    sgnMerge = .collapseSamples(getSignal(object))
    # prepare ranges
    rng = getRanges(object)

    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    # handle all plus ranges
    cRangePlus = subset(rng, strand == "+")
    if (length(cRangePlus) > 0) {
        bsSumPlus = sum(sgnMerge$signalPlus[cRangePlus])
        # make flanking range
        if (flank == "bs") {
            extendedRangePlus = cRangePlus + width(cRangePlus)
        }
        if (flank == "manual") {
            if (is.null(flank.size)) {
                msg = paste0("If 'flank' is set to be 'manual', then provide a valide flank.size as numeric value reflecting in nucleotides, or use option 'bs'.\n")
            } else {
                extendedRangePlus = cRangePlus + flank.size
            }
        }
        # get sum over flanking range
        exSumPlus = sum(sgnMerge$signalPlus[extendedRangePlus])
        mcols(extendedRangePlus)$signalToFlankRatio = (bsSumPlus / exSumPlus)
    } else {
        msg = paste0("No ranges on '+' strand. \n")
        if(!quiet) warning(msg)
        extendedRangePlus = cRangePlus
    }
    # handle all minus ranges
    cRangeMinus = subset(rng, strand == "-")
    if (length(cRangeMinus) > 0) {
        bsSumMinus = sum(sgnMerge$signalMinus[cRangeMinus])
        # make flanking range
        if (flank == "bs") {
            extendedRangeMinus = cRangeMinus + width(cRangeMinus)
        }
        if (flank == "manual") {
            if (is.null(flank.size)) {
                msg = paste0("If 'flank' is set to be 'manual', then provide a valide flank.size as numeric value reflecting in nucleotides, or use option 'bs'.\n")
            } else {
                extendedRangeMinus = cRangeMinus + flank.size
            }
        }
        # get sum over flanking range
        extendedRangeMinus = cRangeMinus + cRangeMinus$bsSize
        exSumMinus = sum(sgnMerge$signalMinus[extendedRangeMinus])
        mcols(extendedRangeMinus)$signalToFlankRatio = (bsSumMinus / exSumMinus)
    } else {
        msg = paste0("No ranges on '-' strand. \n")
        if(!quiet) warning(msg)
        extendedRangeMinus = cRangeMinus
    }

    rngNew = c(extendedRangePlus, extendedRangeMinus)
    rngNew = .sortRanges(rngNew)
    object = setRanges(object, rngNew)

    # ---
    # Store for results
    rng = getRanges(object)
    resultLine = data.frame(
        funName = "calculateSignalToFlankScore()", class = "estimate",
        nIn = length(rng), nOut = length(rng),
        per = paste0(round(length(rng)/ length(rng), digits = 2)*100,"%"),
        options = paste0("flank=", optstr$flank)
    )
    object@results = rbind(object@results, resultLine)

    return(object)
}
