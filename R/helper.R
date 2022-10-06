.sortRanges <- function(rng) {
    rngSort = GenomeInfoDb::sortSeqlevels(rng)
    rngSort = sort(rngSort)
    return(rngSort)
}

#' Collapse signal from replicates
#'
#' Collapses all replicates merges all samples from a \link{BSFDataSet} object
#' into a single signal stream, only split by minus and plus strand.
#'
#' @param object a \code{BSFDataSet} object
#' @param collapseAll TRUE/FALSE, if all samples should be collapsed (TRUE), or
#' if samples should be kept separate by condition (FALSE)
#' @return object of type \code{\link{BSFDataSet}} with updated signal
#'
#' @seealso \code{\link{BSFDataSet}}
#'
#' @examples
#'
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' bdsNew = collapseReplicates(bds)
#'
#' @export
collapseReplicates <- function(object, collapseAll = FALSE) {
    # stop if object is not of class BSF
    stopifnot(is(object, "BSFDataSet"))
    # stop if collapseAll parameter is not logical
    if (!is.logical(collapseAll)) {
        stop("Option collapseAll not logical, set to TRUE/FALSE")
    }

    sgn = getSignal(object)
    mta = getMeta(object)

    # collapse per condition in meta table
    if (!isTRUE(collapseAll)) {
        # get conditions to split by
        cond = levels(mta$condition)
        # handle plus strand
        plus = sgn$signalPlus
        idx = lapply(cond, function(currCond){
            grep(currCond, names(plus))
        })
        plusSum = lapply(seq_along(idx), function(x){
            currSampleIdx = idx[[x]]
            p = sgn$signalPlus[currSampleIdx]
            pSum = 0
            for (i in seq_along(p)) {
                pSum = pSum + p[[i]]
            }
            names(pSum) = names(p[[1]])
            return(pSum)
        })
        names(plusSum) = cond
        # handle minus strand
        minus = sgn$signalMinus
        idx = lapply(cond, function(currCond){
            grep(currCond, names(minus))
        })
        minusSum = lapply(seq_along(idx), function(x){
            currSampleIdx = idx[[x]]
            m = sgn$signalMinus[currSampleIdx]
            mSum = 0
            for (i in seq_along(m)) {
                mSum = mSum + m[[i]]
            }
            names(mSum) = names(m[[1]])
            return(mSum)
        })
        names(minusSum) = cond
        mrgSgn = list(signalPlus = (plusSum), signalMinus = (minusSum))
        newObj = setSignal(object, mrgSgn)
    }

    # collapse all replicates regardless of the conditions
    if (isTRUE(collapseAll)) {

        p = sgn$signalPlus
        m = sgn$signalMinus

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

        # set condition name
        lp = list(pSum)
        names(lp) = "All"
        lm = list(mSum)
        names(lm) = "All"

        mrgSgn = list(signalPlus = (lp), signalMinus = (lm))
        newObj = setSignal(object, mrgSgn)
    }

    return(newObj)
}

# for internal use only
.collapseSamples <- function(signal) {
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

.checkForDropSeqlevels <- function(ranges, signal, dropSeqlevels) {
    msg = NULL
    # check input signal
    if (isTRUE(dropSeqlevels)) {
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
        rangesNew = subset(ranges, match(seqnames(ranges), rngChrsNew))
        # fix signal
        signalNew = lapply(signal, function(currStrand) {
            lapply(currStrand, function(currSample) {
                currSample[match(chrsToUse, names(currSample))]
            })
        })
        # logg which chromosome was removed
        removedChr = setdiff(unique(as.character(unlist(l))), chrsToUse)

        # check if signal was modified
        if (!identical(signal, signalNew)) {
            msg = paste0("Fixed signal input, removing chr: ",
                         paste(removedChr, collapse = " "))
        }
        if (!identical(ranges, rangesNew)) {
            msg = paste0("Fixed ranges input, removing chr: ",
                         paste(removedChr, collapse = " "))
        }
    }
    if (!isTRUE(dropSeqlevels)) {
        rngChrs = sort(as.character(unique(seqnames(ranges))))
        check = lapply(signal, function(currStrand){
            lapply(currStrand, function(currSample){
                currChrs = sort(names(currSample))
                identical(currChrs, rngChrs)
            })
        })
        if(!all(unlist(check))) {
            check = lapply(signal, function(currStrand){
                lapply(currStrand, function(currSample){
                    currChrs = sort(names(currSample))
                    currChrs
                })
            })
            l = append(list(), c(check$signalPlus, check$signalMinus,
                                 list('rngChr' = rngChrs)))
            chrsToUse = Reduce(intersect, l)
            # logg which chromosome was removed
            removedChr = setdiff(unique(as.character(unlist(l))), chrsToUse)
            # not all identical
            msg = paste0("dropSeqlevels is FALSE and chromosome found in ",
                           "ranges and singal do not match on chr: ",
                         removedChr)
            warning(msg)
        }
        rangesNew = ranges
        signalNew = signal
    }
    fixed = list(ranges = rangesNew, signal = signalNew, msg = msg)
    return(fixed)
}


