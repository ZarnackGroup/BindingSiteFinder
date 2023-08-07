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
        nRepsDf = data.frame(defaultNreps = userNreps, applyTo = levels(userCond))
        idx = match(qSel$sel, nRepsDf$applyTo)
        qSel$defaultNreps = nRepsDf$defaultNreps[idx]
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
        nRepsDf = data.frame(defaultNreps = userNreps, applyTo = levels(userCond))
        idx = match(qSel$sel, nRepsDf$applyTo)
        qSel$defaultNreps = nRepsDf$defaultNreps[idx]
        qSel$value = round(qSel$value, digits = 0)
        return(qSel)
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

    c0 = rowSums(as.data.frame(mcols(coverageOverRanges(
        object, returnOptions = "merge_positions_keep_replicates",
        quiet = TRUE))))

    objMod = object
    objMod = suppressMessages(setRanges(objMod, getRanges(object) + flankSize))
    c1 = rowSums(as.data.frame(mcols(coverageOverRanges(
        objMod, returnOptions = "merge_positions_keep_replicates",
        quiet = TRUE))))

    # use 0.1 as pseudocount for the score
    score = c0 / (((c1 - c0) + 0.1) / 2)
    # report median over all binding sites
    score = median(score)
    return(score)
}

.resolveGeneOverlapsWithRule <- function(rng, ols, rule, selectID, selectName, selectType) {
    # initialize local variables
    bsIndex <- geneType <- choice <- NULL

    # takes a Hits object from the binding site and gene annotation matching
    # and a rule by which binding sites on multiple genes are matched
    # and the initial ranges
    # returns the initial binding site ranges with matched gene info

    # check if input is complete
    # -> note that selectType cannot be missing since then this option is not
    # available !!!
    if (is.null(selectID)) {
        # the gene_id is missing
        olsInfo = data.frame(geneIndex = queryHits(ols),
                             bsIndex = subjectHits(ols),
                             geneName = selectName[queryHits(ols)],
                             geneType = selectType[queryHits(ols)]) %>%
            group_by(bsIndex) %>%
            arrange(bsIndex) %>%
            mutate(choice = rule$idx[match(geneType, rule$gene_type)]) %>%
            arrange(choice, .by_group = TRUE) %>%
            slice_head(n = 1)
        # matching index
        idx = match(rng$currIdx, olsInfo$bsIndex)
        # handle intergenic cases
        rngIntergenic = rng[which(is.na(idx))]
        rngIntergenic$geneName = NA
        rngIntergenic$geneType = "Intergenic"
        # handle normal cases
        rngRest = rng[countOverlaps(rng, rngIntergenic) == 0]
        rngRest$geneName[idx[!is.na(idx)]] = olsInfo$geneName
        rngRest$geneType[idx[!is.na(idx)]] = olsInfo$geneType
    }
    if (is.null(selectName)) {
        # the gene_name is missing
        olsInfo = data.frame(geneIndex = queryHits(ols),
                             bsIndex = subjectHits(ols),
                             geneID = selectID[queryHits(ols)],
                             geneType = selectType[queryHits(ols)]) %>%
            group_by(bsIndex) %>%
            arrange(bsIndex) %>%
            mutate(choice = rule$idx[match(geneType, rule$gene_type)]) %>%
            arrange(choice, .by_group = TRUE) %>%
            slice_head(n = 1)
        # matching index
        idx = match(rng$currIdx, olsInfo$bsIndex)
        # handle intergenic cases
        rngIntergenic = rng[which(is.na(idx))]
        rngIntergenic$geneID = NA
        rngIntergenic$geneType = "Intergenic"
        # handle normal cases
        rngRest = rng[countOverlaps(rng, rngIntergenic) == 0]
        rngRest$geneID[idx[!is.na(idx)]] = olsInfo$geneID
        rngRest$geneType[idx[!is.na(idx)]] = olsInfo$geneType
    }
    if (! is.null(selectID) & ! is.null(selectName)) {
        # everything is in place
        olsInfo = data.frame(geneIndex = queryHits(ols),
                             bsIndex = subjectHits(ols),
                             geneID = selectID[queryHits(ols)],
                             geneName = selectName[queryHits(ols)],
                             geneType = selectType[queryHits(ols)]) %>%
            group_by(bsIndex) %>%
            arrange(bsIndex) %>%
            mutate(choice = rule$idx[match(geneType, rule$gene_type)]) %>%
            arrange(choice, .by_group = TRUE) %>%
            slice_head(n = 1)
        # matching index
        idx = match(rng$currIdx, olsInfo$bsIndex)
        # handle intergenic cases
        rngIntergenic = rng[which(is.na(idx))]
        rngIntergenic$geneID = NA
        rngIntergenic$geneName = NA
        rngIntergenic$geneType = "Intergenic"
        # handle normal cases
        rngRest = rng[countOverlaps(rng, rngIntergenic) == 0]
        rngRest$geneID[idx[!is.na(idx)]] = olsInfo$geneID
        rngRest$geneName[idx[!is.na(idx)]] = olsInfo$geneName
        rngRest$geneType[idx[!is.na(idx)]] = olsInfo$geneType
    }

    rng = .sortRanges(c(rngIntergenic, rngRest))

    return(rng)
}


# approximate by one round of inital binding site merging
# approximate by center given by highest number of stacked crosslinks
.approximateBindingSites_medium <- function(rng, sgn, bsSize, minWidth) {
    # approximate binding sites by first iteration of merging
    res = .mergeCrosslinkSites(
        rng = rng,
        sgn = sgn,
        bsSize = bsSize,
        minWidth = minWidth,
        computeOption = "simple"
    )
    rng = res$rng
    rng$bsSize = bsSize
    return(rng)
}

# approximate center by highest pureclip score
#' @importFrom IRanges extractList
.approximateBindingSites_coarse <- function(rng, bsSize, minWidth) {
    # initialize variables locally
    i <- group <- NULL

    ### Merge peaks for given bs size
    redRegion = reduce(rng, min.gapwidth = 1, with.revmap = TRUE)

    redRegion = redRegion[width(redRegion) >= minWidth]
    mcols(redRegion)$id = seq_along(redRegion)

    scoresPerRedRegion = extractList(mcols(rng), redRegion$revmap)
    names(scoresPerRedRegion) = seq_along(scoresPerRedRegion)

    scoresPerRedRegion = unlist(scoresPerRedRegion)
    scoresPerRedRegion$group = rownames(scoresPerRedRegion)

    idx = scoresPerRedRegion %>% as.data.frame() %>%
        group_by(group) %>% summarise(i = which.max(score)) %>%
        mutate(group = as.numeric(group)) %>% arrange(group)

    # Convert the GRanges object to a data frame
    redRegionDf <- data.frame(seqnames = seqnames(redRegion),
                              start = start(redRegion),
                              end = end(redRegion),
                              strand = strand(redRegion),
                              group = redRegion$id,
                              i = idx$i) %>%
        mutate(center = start + i -1)

    rngNew = GRanges(seqnames = redRegionDf$seqnames, ranges = IRanges::IRanges(start = redRegionDf$center, width = floor(bsSize/2)),
                     strand = redRegionDf$strand, bsSize = bsSize)
    return(rngNew)
}

.findFirstMaximum <- function(x, est.minimumStepGain){

    est.option <- df.global <- bsFirstDrop <- idxFirstDrop <- bsSize <- NULL
    localEstimates <- est.option <- est.bsSize <- est.geneFilter <- NULL
    signalToFlankRatio <- ms <- growth_per <- increaseOverMin <- geneWiseFilter <- NULL


    df.global = x %>%
        group_by(bsSize) %>%
        summarise(ms = mean(signalToFlankRatio), sd = sd(signalToFlankRatio)) %>%
        mutate(geneWiseFilter = "mean") %>%
        mutate(growth_per = (ms / dplyr::lag(ms)) - 1) %>%
        mutate(increase = growth_per > 0) %>%
        mutate(increaseOverMin = (growth_per) >= est.minimumStepGain) %>%
        slice(-1)

    # check for minimum increase
    if (all(!df.global$increaseOverMin)) {
        # no size leads to score increase larger than est.minimumStepGain
        # -> trigger local estimation
        localEstimates = .findFirstMaximum_local(x, est.minimumStepGain)
        est.option = localEstimates$option
        est.bsSize = localEstimates$est.bsSize
        est.geneFilter = localEstimates$est.geneFilter

    } else {
        # at lease one size leads to score increase larger than est.minimumStepGain
        # -> trigger global estimation

        # finde die zeile in der ein TRUE auf ein FALSE folgt
        bsFirstDrop = df.global %>%
            filter(increaseOverMin == FALSE & dplyr::lag(increaseOverMin, default = TRUE) == TRUE) %>%
            pull(bsSize)

        # check for no possible maximum
        if (length(bsFirstDrop) == 0) {
            est.option = "error"
            res = list(est.bsSize = est.bsSize, est.geneFilter = est.geneFilter, est.option = est.option)
            return(res)
        } else {
            bsFirstDrop = max(bsFirstDrop)
        }


        # firstMax = which(df.global$increaseOverMin == TRUE)[1]
        # df.current = df.global %>% slice(1:firstMax)

        # firstDrop = which(df.global$increaseOverMin == FALSE)[1]

        if (is.na(bsFirstDrop)) {
            localEstimates = .findFirstMaximum_local(x, est.minimumStepGain)
            est.option = localEstimates$option
            est.bsSize = localEstimates$est.bsSize
            est.geneFilter = localEstimates$est.geneFilter

        } else {
            # pull bsSize
            idxFirstDrop = which(df.global$bsSize == bsFirstDrop)
            est.bsSize = df.global %>% slice(idxFirstDrop-1) %>% pull(bsSize)

            est.geneFilter = x %>%
                filter(bsSize == est.bsSize & signalToFlankRatio >= df.global$ms[df.global$bsSize == est.bsSize]) %>%
                arrange(geneWiseFilter) %>%
                slice_head(n = 1) %>%
                pull(geneWiseFilter)
            # set option
            est.option = "global"
        }

    }

    res = list(est.bsSize = est.bsSize, est.geneFilter = est.geneFilter, est.option = est.option)

    return(res)
}

.findFirstMaximum_local <- function(x, est.minimumStepGain) {

    df.local <- option <- est.bsSize <- est.geneFilter <- NULL
    geneWiseFilter <- signalToFlankRatio <- growth_per <- increase <- bsSize <- NULL

    df.local = x %>%
        group_by(geneWiseFilter) %>%
        mutate(growth_per = (signalToFlankRatio / dplyr::lag(signalToFlankRatio)) - 1) %>%
        mutate(increase = growth_per > 0) %>%
        mutate(increaseOverMin = (growth_per) >= est.minimumStepGain) %>%
        filter(!is.na(growth_per))

    if (any(df.local$increase)) {
        # check if at any point an increase of the curve was found
        # -> if so take it as anchor
        df.curr = df.local %>%
            filter(increase == TRUE) %>%
            ungroup() %>%
            filter(growth_per == max(growth_per))
        # pull bsSize
        est.bsSize = df.curr %>% pull(bsSize)
        # pull geneWiseFilter
        est.geneFilter = df.curr %>% pull(geneWiseFilter)
        # set option
        option = "local"
    }
    if (all(!df.local$increase)) {
        # check if there is not a single increase in the curve at any point
        # -> if so use the overall maximum as anchor
        # pull bsSize
        est.bsSize = df.local$bsSize[which.max(df.local$signalToFlankRatio)]
        # pull geneWiseFilter
        est.geneFilter = df.local$geneWiseFilter[which.max(df.local$signalToFlankRatio)]
        # set option
        option = "fallback"
    }
    # return results
    res = list(option = option, est.bsSize = est.bsSize, est.geneFilter = est.geneFilter)
    return(res)
}

.calcNormalizeFactors <- function(object, rng, trl, normalize.exclude.lower,
                                  normalize.exclude.upper){

    # ---
    # Checks
    # ---
    # compute normalization factors for only those regions that actually have binding sites
    this.trl = trl[toupper(names(trl)) %in% toupper(unique(rng$transcriptRegion))]

    # ---
    # Hosting
    # ---
    # -> compute width for all ranges that overlap with at least one binding site
    w.hosting = lapply(seq_along(this.trl), function(x){
        # get all relevant width
        curr.width = width(reduce(IRanges::subsetByOverlaps(this.trl[[x]], rng)))
        # find upper and lower boundary cutoffs
        currIdx = findInterval(curr.width, vec = quantile(x = curr.width, probs = seq(from = 0, to = 1, by = 0.01)), all.inside = TRUE)
        seqs = seq(from = 0.01, to = 1, by = 0.01)
        names(seqs) = seq_along(seqs)
        idxRemoveLower = as.numeric(names(seqs[seqs <= normalize.exclude.lower]))
        idxRemoveUpper = as.numeric(names(seqs[seqs >= 1-normalize.exclude.upper+0.01]))
        # apply cutoffs
        curr.w.reduced = curr.width[which(!currIdx %in% c(idxRemoveLower, idxRemoveUpper))]
        # convert to normalization values
        width.return = c(mean(curr.w.reduced), median(curr.w.reduced), sum(curr.w.reduced))
        return(width.return)
    })
    w.hosting = do.call(rbind, w.hosting)
    colnames(w.hosting) = c("norm.hosting.mean", "norm.hosting.median", "norm.hosting.sum")
    w.hosting = as.data.frame(w.hosting)
    # set names
    w.hosting$TranscriptRegion = names(this.trl)

    w.total = cbind.data.frame(w.hosting)
    return(w.total)
}

