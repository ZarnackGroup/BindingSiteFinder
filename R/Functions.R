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
#' @param cutoff numeric; percentage cutoff to be used for the
#' reproducibility quantile filtering
#' @param nReps numeric; number of replicates that must meet the cutoff
#' defined in \code{cutoff} for a binding site to be called reproducible.
#' Defaults to N-1.
#' @param minCrosslinks numeric; minimal number of crosslinks a binding
#' site needs to have to be called reproducible. Acts as a lower boundary for
#' \code{cutoff}. Defaults to 1.
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
        object, returnOptions = "merge_positions_keep_replicates")))

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

    idx = match(rng$bsID, rngNew$bsID)
    mcols(rng) = mcols(rngNew[idx])

    object = setRanges(object, rng)

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



#' Combine multiple \code{\link{BSFDataSet}} objects
#'
#' This function combines all \code{\link{BSFDataSet}} objects from the input
#' list into a single \code{\link{BSFDataSet}} object.
#'
#' Meta-data tables are added to each other by performing a row-wise bind,
#' basically adding all meta data tables underneath each other.
#'
#' The default way of signal combination is merging all signal lists on the level
#' of the indivdual samples. One can also force a re-load of the signal list
#' component by using \code{force.reload=TRUE}.
#' The signal can be combined by
#'
#' The ranges are combined by adding both granges objects together. With option
#' \code{overlaps.fix} one can decide if partially overlapping ranges should be
#' combined into a single range or not. If this option is FALSE one is likely to
#' have overlapping binding sites after the merge. If this option is TRUE, then
#' the combined coverage is used to guide the new center point for these cases.
#'
#' The \code{combine.bsSize} option allows one to set a unique bsSize for all
#' objects that should be combined. Although it is recommended to combine only
#' objects with the same bsSize this option can be used to ensure that the
#' merged result has the same bsSize for all ranges.
#'
#' @param list list; a list of objects from class \code{\link{BSFDataSet}} that
#' should be combined
#' @param overlaps.fix logical; if partially overlapping binding sites should be
#' re-centered
#' @param combine.bsSize numeric; the binding site size that the merged sites
#' should have. Default=NULL, then bsSize is taken from the input objects
#' in \code{list}.
#' @param combine.name character; meta table name of the combined object.
#' Default=NULL; then name is set to 'combined'
#' @param force.reload logical; whether the signal should be derived from the
#' merge of the input objects given in \code{list} or if the signal should be
#' re-loaded from the path given in the meta data.
#' @param quiet logical; whether to print messages or not
#' @param veryQuiet logical; whether to print status messages or not
#'
#' @return an object of class \code{\link{BSFDataSet}} with ranges, signal and
#' meta data resulting from the merge of the input objects.
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' # make binding sites
#' bds = makeBindingSites(bds, bsSize = 7)
#'
#' # split ranges in two groups
#' allRanges = getRanges(bds)
#' set.seed(1234)
#' idx = sample(1:length(allRanges), size = length(allRanges)/2, replace = FALSE)
#' r1 = allRanges[idx]
#' r2 = allRanges[-idx]
#'
#' # splite meta data
#' allMeta = getMeta(bds)
#' m1 = allMeta[1:2,]
#' m2 = allMeta[3:4,]
#'
#' # create new objects
#' bds1 = setRanges(bds, r1)
#' bds2 = setRanges(bds, r2)
#' bds1 = setMeta(bds1, m1)
#' bds2 = setMeta(bds2, m2)
#' bds1 = setName(bds1, "test1")
#' bds2 = setName(bds2, "test2")
#'
#' # merge two objects with '+' operator
#' c1 = bds1 + bds2
#'
#' # merge two objects from list
#' list = list(bds1, bds2)
#' c1 = combineBSF(list = list, overlaps.fix = TRUE,
#'  combine.bsSize = NULL, combine.name = NULL, quiet = TRUE)
#'
#' @export
combineBSF <- function(list, # list of class BSFDataSet
                       overlaps.fix = TRUE,
                       combine.bsSize = NULL,
                       combine.name = NULL,
                       force.reload = FALSE,
                       quiet = TRUE,
                       veryQuiet = FALSE
){
    # initialize local variables
    this.bsSize <- resizedHit <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(all(unlist(lapply(list, is, "BSFDataSet"))))
    if (!is.null(combine.bsSize)){
        stopifnot(is(combine.bsSize, "numeric"))
    }

    # check dataset names
    list = lapply(seq_along(list), function(x){
        this.dataset = list[[x]]
        this.name = getName(list[[x]])
        if (length(this.name) == 0) {
            # No name was given to this dataset
            # -> make name and set
            new.name = paste0("bds", x)
            this.dataset = setName(list[[x]], new.name)

            # Inform user
            msg0 = paste0("No name present for dataset No: ", x, ".\n")
            msg1 = paste0("Name is set to: ", new.name, ".\n")
            msg2 = paste0("Prevent this message by setting a name using 'setNames()'.\n")
            if (!quiet) warning(c(msg0, msg1, msg2))
        }
        return(this.dataset)
    })

    # Combine meta
    # --------------------------------------------------------------------------
    comb.meta = do.call(rbind, lapply(list, getMeta))

    # check for duplicated path
    if (any(duplicated(comb.meta$clPlus))) {
        msg = paste0("Duplicated path for plus strand found. Please check input.\n")
        if (!quiet) warning(c(msg))
    }
    if (any(duplicated(comb.meta$clMinus))) {
        msg = paste0("Duplicated path for minus strand found. Please check input.\n")
        if (!quiet) warning(c(msg))
    }
    # set ids
    comb.meta$id = 1:nrow(comb.meta)

    # set name
    comb.meta$datasets = comb.meta$name
    if (is.null(combine.name)) {
        comb.meta$name = "combined"
    } else {
        comb.meta$name = combine.name
    }

    # Combine signal
    # --------------------------------------------------------------------------
    # combine signal if user wants it or not
    if (!isTRUE(force.reload)) {
        # combine plus strand
        signal.plus = lapply(list, function(x){
            this.signal = getSignal(x)
            this.signal.plus = this.signal$signalPlus
            return(this.signal.plus)
        })
        signal.plus = unlist(signal.plus)
        names(signal.plus) = paste0(comb.meta$id, "_", comb.meta$condition)
        # combine minus strand
        signal.minus = lapply(list, function(x){
            this.signal = getSignal(x)
            this.signal.minus = this.signal$signalMinus
            return(this.signal.minus)
        })
        signal.minus = unlist(signal.minus)
        names(signal.minus) = paste0(comb.meta$id, "_", comb.meta$condition)

        comb.signal = list(signalPlus = signal.plus,
                           signalMinus = signal.minus)
    }

    # Combine ranges
    # --------------------------------------------------------------------------
    comb.ranges = lapply(seq_along(list), function(x){
        this.range = getRanges(list[[x]])
        mcols(this.range)$dataset = getName(list[[x]])
        return(this.range)
    })
    comb.ranges = do.call(c, comb.ranges)

    # fix overlaps or not
    if (!isTRUE(overlaps.fix)) {
        msg = paste0("Overlaps are not resoloved. This could lead to overlapping binding site.\n")
        if (!quiet) warning(msg)
    }

    # check for identical bsSize
    if (!length(unique(width(comb.ranges))) == 1) {
        # not all ranges have the same size
        msg0 = paste0("Not all ranges are of the same size. \n")
        msg1 = paste0("Found ranges of bsSize = c(", paste(unique(width(comb.ranges)), collapse = ","), ").\n")
        # do not fix ranges
        if (!isTRUE(overlaps.fix)) {
            if (!quiet) warning(c(msg0, msg1))
        } else {
            # ranges should be fixed
            if (is.null(combine.bsSize)) {
                msg2 = paste0("Use 'combine.bsSize' to set a reference size for all ranges.\n")
                stop (c(msg0, msg1, msg2))
            } else {
                msg2 = paste0("Manual size of ", combine.bsSize, "nt is used for merging.\n")
                this.bsSize = combine.bsSize
            }
            if (!quiet) warning(c(msg0, msg1, msg2))
        }
    } else {
        # all ranges have equal size
        if (unique(width(comb.ranges)) == 1) {
            # ranges to combine are of 1 nt range
            msg0 = paste0("Range to combine have width of 1nt.\n")
            msg1 = paste0("Make sure to run binding site definition first; see `BSFind()`.\n")
            msg2 = paste0("Alternativly set `overlaps.fix=FALSE` to simply add ranges.\n")
            if (!quiet) warning(c(msg0,msg1,msg2))
        }
        if (!is.null(combine.bsSize)) {
            # set combined bsSize manually
            this.bsSize = combine.bsSize
        } else {
            # set combined bsSize from data
            this.bsSize = unique(width(comb.ranges))
        }
    }

    # ---
    # Store function parameters in list
    optstr = list(overlaps.fix = overlaps.fix, combine.bsSize = this.bsSize,
                  combine.name = combine.name, force.reload = force.reload)

    # Resize ranges overlaps
    # --------------------------------------------------------------------------
    # don't fix overlaps, keep as it is
    if (!isTRUE(overlaps.fix)) {
        # apply new bsSize for all ranges if needed
        if (!is.null(this.bsSize)) {
            # turn ranges into midpoint representations
            comb.ranges = resize(comb.ranges, fix = "center", width = this.bsSize)
        }
        # create combined BSFDataSet
        if (!veryQuiet) message("Creating merged object...")
        # use reload option or not
        if (isTRUE(force.reload)){
            resize.bds = BSFDataSetFromBigWig(ranges = comb.ranges, meta = comb.meta, silent = quiet)
        } else {
            resize.bds = BSFDataSet(ranges = comb.ranges, meta = comb.meta, signal = comb.signal, silent = quiet)
        }
    } else {
        # turn ranges into midpoint representations
        resize.ranges = resize(comb.ranges, fix = "center", width = 1)
        # create combined BSFDataSet
        if (!veryQuiet) message("Creating merged object...")
        # use reload option or not
        if (isTRUE(force.reload)) {
            resize.bds = BSFDataSetFromBigWig(ranges = resize.ranges, meta = comb.meta, silent = quiet)
        } else {
            resize.bds = BSFDataSet(ranges = resize.ranges, meta = comb.meta, signal = comb.signal, silent = quiet)
        }
        # resolve overlaps by merging
        if (!veryQuiet) message("Fixing partial overlaps...")
        resize.bds = makeBindingSites(object = resize.bds, bsSize = this.bsSize,
                                      minWidth = 0, minCrosslinks = 0, minClSites = 0,
                                      centerIsSummit = FALSE, centerIsClSite = FALSE,
                                      quiet = TRUE)
        # remove stats from makeBiningSites
        resize.bds@results = data.frame()
        resize.bds@params = list()
        resize.bds@plotData = list()
    }

    # Add meta information from original ranges
    # --------------------------------------------------------------------------
    # match resized binding sites with all input binding sites
    if (!veryQuiet) message("Merging meta info...")

    resized.ranges = getRanges(resize.bds)
    ols = findOverlaps(resized.ranges, comb.ranges)

    matched.ranges = resized.ranges[queryHits(ols)]
    matched.ranges$combHit = subjectHits(ols)
    matched.ranges$resizedHit = queryHits(ols)

    mcols(matched.ranges) = cbind(combHit = matched.ranges$combHit,
                                  resizedHit = matched.ranges$resizedHit,
                                  mcols(comb.ranges[matched.ranges$combHit]))

    # group by matching
    grouped.meta = mcols(matched.ranges) %>%
        as.data.frame(, row.names = NULL) %>%
        group_by(resizedHit)

    # match meta data from overlapping binding sites
    total.meta = .matchMetaData(grouped.meta)

    # attach meta data to ranges
    fix.ranges = resized.ranges
    mcols(fix.ranges) = cbind(mcols(fix.ranges), total.meta)

    # manage other columns
    mcols(fix.ranges)$bsID = paste0("BS", seq_along(fix.ranges))
    if ("bsSize" %in% colnames(mcols(resized.ranges))) {
        mcols(fix.ranges)$bsSize = max(resized.ranges$bsSize)
    }
    mcols(fix.ranges)$resizedHit = NULL

    # set final object
    comb.bds = setRanges(resize.bds, fix.ranges)

    # ---
    # Store for results
    resultLine = data.frame(
        funName = "combineBSF()", class = "transform",
        nIn = length(comb.ranges), nOut = length(fix.ranges),
        per = paste0(round(length(fix.ranges)/ length(comb.ranges), digits = 2)*100,"%"),
        options = paste0("overlaps.fix=", optstr$overlaps.fix, ", combine.bsSize=", optstr$combine.bsSize)
    )
    comb.bds@results = rbind(comb.bds@results, resultLine)
    comb.bds@params$combineBSF = optstr

    return(comb.bds)
}
