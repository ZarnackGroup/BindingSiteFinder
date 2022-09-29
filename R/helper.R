.sortRanges <- function(rng) {
    rngSort = GenomeInfoDb::sortSeqlevels(rng)
    rngSort = sort(rngSort)
    return(rngSort)
}

# not exported functions that are for internal use only
.collapesReplicates <- function(signal) {
    p = signal$signalPlus
    m = signal$signalMinus

    pSum = 0
    for (i in seq_along(p)) {
        pSum = pSum + p[[i]]
    }
    names(pSum) = names(p[[1]])
    mSum = 0
    for (i in seq_along(m)) {
        mSum = mSum + m[[i]]
    }
    names(mSum) = names(m[[1]])

    mergedSignal = list(signalPlus = (pSum), signalMinus = (mSum))
    return(mergedSignal)
}

#' @importFrom stats quantile
.selectQuantilesMultipleConditions <-
    function(covDf, userCond, userNreps, userCutoff) {
        # bind locally used variables
        per <- applyTo <- NULL

        # construct general cutoff matrix in 1% steps
        q = as.data.frame(apply(covDf, 2, function(x) {
            quantile(x, probs = seq(0, 1, by = 0.01))
        }))
        q$cut = seq(0, 1, by = 0.01)
        q$per = rownames(q)

        # select that part of q that was chosen by user
        qSel = q[q$cut %in% userCutoff, ]
        applyDf = data.frame(levels(userCond), userCutoff)
        # applyDf = data.frame(levels(userCond), userCutoff, userCond)

        idx = match(qSel$cut, applyDf$userCutoff)
        qSel$applyTo = applyDf$levels.userCond.[idx]
        qSel = qSel %>% pivot_longer(-c(cut, per, applyTo))
        qSel$sel = vapply(strsplit(qSel$name, "_"), `[`, 2,
                          FUN.VALUE = character(1))
        if (length(unique(userCutoff)) > 1) {
            qSel = qSel[qSel$applyTo == qSel$sel, ]
        }
        # add n.reps support to df
        nRepsDf = data.frame(n.reps = userNreps, applyTo = levels(userCond))
        idx = match(qSel$sel, nRepsDf$applyTo)
        qSel$n.reps = nRepsDf$n.reps[idx]
        return(qSel)
    }

#' @importFrom stats quantile
.selectQuantilesSingleCondtion <-
    function(covDf, userCond, userNreps, userCutoff) {
        # bind locally used variables
        per <- NULL

        # construct general cutoff matrix in 1% steps
        q = as.data.frame(apply(covDf, 2, function(x) {
            quantile(x, probs = seq(0, 1, by = 0.01))
        }))
        q$cut = seq(0, 1, by = 0.01)
        q$per = rownames(q)

        # select that part of q that was chosen by user
        qSel = q[q$cut %in% userCutoff, ]
        qSel = qSel %>% pivot_longer(-c(cut, per))
        qSel$sel = vapply(strsplit(qSel$name, "_"), `[`, 2,
                          FUN.VALUE = character(1))

        # add n.reps support to df
        nRepsDf = data.frame(n.reps = userNreps, applyTo = levels(userCond))
        idx = match(qSel$sel, nRepsDf$applyTo)
        qSel$n.reps = nRepsDf$n.reps[idx]
        return(qSel)
    }

.subsetByChr <- function(object, chr) {
    # subset ranges
    rng = getRanges(object)
    rngSub = rng[seqnames(rng) == chr,]

    # subset the signal
    sgn = getSignal(object)
    sgnSub = lapply(sgn, function(selStrand) {
        lapply(selStrand, function(chrList) {
            chrList[names(chrList) == chr]
        })
    })

    objectNew = setRanges(object, rngSub)
    objectNew = setSignal(objectNew, sgnSub)
    return(objectNew)
}

#' @importFrom stats median
.computeSupportRatio <- function(object, flankSize){
    stopifnot(is(object, "BSFDataSet"))
    if(! length(flankSize) == 1) {
        stop("flankSize must be a single integer value. ")
    }
    if (! is.numeric(flankSize)) {
        stop("flankSize needs to be numeric. ")
    }
    if (round(flankSize) != flankSize) {
        stop("flankSize is not an integer. ")
    }

    # c0 = rowSums(coverageOverRanges(object, merge = TRUE,
    #                                 returnType = "data.frame"))

    c0 = rowSums(as.data.frame(mcols(coverageOverRanges(
        object, returnOptions = "merge_positions_keep_replicates",
        silent = TRUE))))

    objMod = object
    objMod = setRanges(objMod, getRanges(object) + flankSize)
    # c1 = rowSums(coverageOverRanges(objMod, merge = TRUE,
    #                                 returnType = "data.frame"))
    c1 = rowSums(as.data.frame(mcols(coverageOverRanges(
        objMod, returnOptions = "merge_positions_keep_replicates",
        silent = TRUE))))

    # use 0.1 as pseudocount for the score
    score = c0 / (((c1 - c0) + 0.1) / 2)
    # report median over all binding sites
    score = median(score)
    return(score)
}


#' Helper function that fixes seqnames errors in ranges and signal
#'
#' Ranges and signal list must be in the same structure as it is used by
#' the \link{BSFDataSet}. The function removes all ranges and signal entries
#' where the seqlevels are not part of both signal and ranges.
#'
#' If only ranges are problematic, keepStandardChromosomes might solve the
#' problem too.
#'
#' @param ranges GRanges object as it is used to build a \link{BSFDataSet}
#' @param signal SignalList object as it is used to build a \link{BSFDataSet}
#'
#' @return GRanges list for ranges and SignalList for signal
#'
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomeInfoDb dropSeqlevels
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' rng = getRanges(bds)
#' sgn = getSignal(bds)
#'
#' r = as.data.frame(granges(rng))
#' d = data.frame(seqnames = "chr6", start = 1, end = 1, width = 1, strand = "+")
#' rngNew = makeGRangesFromDataFrame(rbind(d, r), keep.extra.columns = TRUE)
#'
#' fixed = forceEqualNames(rngNew, sgn)
#'
#' bdsNew = setRanges(bds, fixed$ranges)
#' bdsNew = setSignal(bdsNew, fixed$signal)
#'
#' @export
forceEqualNames <- function(ranges, signal) {
    # check which chromosomes do not fit
    rngChrs = sort(as.character(unique(seqnames(ranges))))
    check = lapply(signal, function(currStrand){
        lapply(currStrand, function(currSample){
            currChrs = sort(names(currSample))
            currChrs
        })
    })
    l = append(list(), c(check$signalPlus, check$signalMinus,
                         list('rngChr' = rngChrs)))
    chrsToUse = Reduce(intersect, l)
    # fix ranges
    rngChrsNew = rngChrs[match(chrsToUse, rngChrs)]
    ranges = subset(ranges, match(seqnames(ranges), rngChrsNew))
    # fix siganl
    signal = lapply(signal, function(currStrand) {
        lapply(currStrand, function(currSample) {
            currSample[match(chrsToUse, names(currSample))]
        })
    })
    # clean seqlevels
    ranges = GenomeInfoDb::dropSeqlevels(ranges,
                                         value = seqlevels(ranges)[
                                             !match(seqlevels(ranges),
                                                    unique(seqnames(ranges)),
                                                    nomatch = 0) > 0])

    fixed = list(ranges = ranges, signal = signal)
    return(fixed)
}
