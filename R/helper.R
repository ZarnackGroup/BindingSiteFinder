.sortRanges <- function(rng) {
    rngSort = GenomeInfoDb::sortSeqlevels(rng)
    rngSort = sort(rngSort)
    return(rngSort)
}

.collapseSamples <- function(signal) {
    pSum <- Reduce(`+`, signal$signalPlus)
    names(pSum) <- names(signal$signalPlus[[1]])
    mSum <- Reduce(`+`, signal$signalMinus)
    names(mSum) <- names(signal$signalMinus[[1]])

    mergedSignal <- list(signalPlus = pSum, signalMinus = mSum)
    return(mergedSignal)
}

.subsetByChr <- function(object, chr, quiet = FALSE) {
    # subset ranges
    rng = getRanges(object)
    # rngSub = rng[seqnames(rng) %in% chr]
    rngSub = subset(rng, match(seqnames(rng), chr))

    # subset the signal
    sgn = getSignal(object)
    sgnSub = lapply(sgn, function(selStrand) {
        lapply(selStrand, function(chrList) {
            chrList[names(chrList) %in% chr]
        })
    })

    objectNew = setRanges(object, rngSub, quiet = quiet)
    objectNew = setSignal(objectNew, sgnSub, quiet = quiet)
    return(objectNew)
}

.subsetByChr_old <- function(object, chr, quiet = FALSE) {
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

    objectNew = setRanges(object, rngSub, quiet = quiet)
    objectNew = setSignal(objectNew, sgnSub, quiet = quiet)
    return(objectNew)
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

.reduceSignalToFrame <- function(object, frame, quiet = FALSE) {
    rng = getRanges(object)
    # extend ranges by desired frame on which signal should be kept
    rng = rng + frame
    newObj = setRanges(object, rng, quiet = quiet)
    # drop signal outside of protective frame
    newObj = newObj[seq_along(rng), drop=TRUE]
    # put original sized ranges back in place
    rng = rng-frame
    newObj = setRanges(newObj, rng, quiet = quiet)

    return(newObj)
}

.capitalize <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
          sep="", collapse=" ")
}

