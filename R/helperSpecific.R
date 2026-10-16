#' @importFrom stats quantile
.selectQuantilesMultipleConditions <-
    function(covDf, userCond, userMeta, userNreps, userCutoff) {
        # returns one row per sample (in the column order of covDf) with the
        # crosslink threshold computed from the cutoff of its condition
        # -> cutoffs and nReps are matched to conditions by position, so that
        # conditions can share the same cutoff
        # -> quantiles are computed directly for the user cutoff (no grid), so
        # that any cutoff can be used (eg. 0.35)
        this.cond = .sampleToCondition(colnames(covDf), userMeta)
        idx = match(this.cond, levels(userCond))
        this.cut = userCutoff[idx]
        this.q = lapply(seq_along(covDf), function(i) {
            quantile(covDf[[i]], probs = this.cut[i])
        })

        qSel = data.frame(cut = this.cut,
                          per = vapply(this.q, names, character(1)),
                          applyTo = this.cond,
                          name = colnames(covDf),
                          value = unname(unlist(this.q)),
                          sel = this.cond,
                          defaultNreps = userNreps[idx])
        return(qSel)
    }

#' @importFrom stats quantile
.selectQuantilesSingleCondtion <-
    function(covDf, userCond, userMeta, userNreps, userCutoff) {
        # returns one row per sample (in the column order of covDf) with the
        # crosslink threshold computed from the single cutoff
        # -> quantiles are computed directly for the user cutoff (no grid), so
        # that any cutoff can be used (eg. 0.35)
        this.cond = .sampleToCondition(colnames(covDf), userMeta)
        idx = match(this.cond, levels(userCond))
        this.q = lapply(seq_along(covDf), function(i) {
            quantile(covDf[[i]], probs = userCutoff)
        })

        qSel = data.frame(cut = userCutoff,
                          per = vapply(this.q, names, character(1)),
                          name = colnames(covDf),
                          value = round(unname(unlist(this.q)), digits = 0),
                          sel = this.cond,
                          defaultNreps = userNreps[idx])
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


.matchMetaData <- function(d){
    # assign local variables
    geneID <- geneName <- geneType <- transcriptRegion <- NULL
    dataset <- signalToFlankRatio <- NULL


    # assign meta data if present
    if ("geneID" %in% colnames(d)) {
        merge.geneID = d %>%
            summarize(geneID = toString(unique(geneID))) %>%
            mutate(geneID = ifelse(geneID == "NA", NA, geneID))
    } else {
        merge.geneID = NA
    }
    if ("geneName" %in% colnames(d)) {
        merge.geneName = d %>%
            summarize(geneName = toString(unique(geneName))) %>%
            mutate(geneName = ifelse(geneName == "NA", NA, geneName))
    } else {
        merge.geneName = NA
    }
    if ("geneType" %in% colnames(d)) {
        merge.geneType = d %>%
            summarize(geneType = toString(unique(geneType))) %>%
            mutate(geneType = ifelse(geneType == "NA", NA, geneType))
    } else {
        merge.geneType = NA
    }
    if ("transcriptRegion" %in% colnames(d)) {
        merge.transcriptRegion = d %>%
            summarize(transcriptRegion = toString(unique(transcriptRegion))) %>%
            mutate(transcriptRegion = ifelse(transcriptRegion == "NA", NA, transcriptRegion))
    } else {
        merge.transcriptRegion = NA
    }
    if ("dataset" %in% colnames(d)) {
        merge.dataset = d %>%
            summarize(dataset = toString(unique(dataset))) %>%
            mutate(dataset = ifelse(dataset == "NA", NA, dataset))
    } else {
        merge.dataset = NA
    }
    if ("score" %in% colnames(d)) {
        merge.score = d %>%
            summarize(score = max(score)) %>%
            mutate(score = ifelse(score == "NA", NA, score))
    } else {
        merge.score = NA
    }
    if ("signalToFlankRatio" %in% colnames(d)) {
        merge.signalToFlankRatio = d %>%
            summarize(signalToFlankRatio = max(signalToFlankRatio)) %>%
            mutate(signalToFlankRatio = ifelse(signalToFlankRatio == "NA", NA, signalToFlankRatio))
    } else {
        merge.signalToFlankRatio = NA
    }

    total.meta = cbind(merge.geneID, merge.geneName, merge.geneType,
                       merge.transcriptRegion, merge.dataset,
                       merge.score, merge.signalToFlankRatio) %>%
        dplyr::distinct_all() %>%
        dplyr::select_if(~sum(!is.na(.)) > 0)
    return(total.meta)
}


.coverageForDifferential <- function(object, ranges, prefix, quiet, match.rangeID){
    this.object = setRanges(object, ranges, quiet = quiet)
    cov = clipCoverage(this.object,
                       match.rangeID = match.rangeID,
                       ranges.merge = TRUE, positions.merge = FALSE,
                       samples.merge = FALSE, out.format = "granges",
                       out.format.overwrite = TRUE, quiet = quiet)
    colnames(mcols(cov)) = paste0(prefix, colnames(mcols(cov)))
    mcols(cov) = cbind(mcols(ranges), mcols(cov))
    return(cov)
}


.sampleToCondition <- function(sample.names, meta) {
    # takes sample names as build by the package ("<id>_<condition>"),
    # optionally with a count prefix ("counts.bs." or "counts.bg.")
    # and the meta data of the object
    # returns the condition for each sample name
    # -> samples are matched by their id and the condition is taken from the
    # meta data and not parsed from the name, so that conditions can contain
    # any character (e.g. "_") and can be renamed after counting
    this.samples = sub("^counts\\.(bs|bg)\\.", "", sample.names)
    this.ids = as.character(meta$id)
    idx = vapply(this.samples, function(x) {
        hit = which(startsWith(x, paste0(this.ids, "_")))
        if (length(hit) == 0) return(NA_integer_)
        # use the longest matching id (e.g. '1_a' over '1')
        hit[which.max(nchar(this.ids[hit]))]
    }, integer(1), USE.NAMES = FALSE)
    if (anyNA(idx)) {
        msg0 = paste0("Sample names do not match the pattern '<id>_<condition>' with an id from the meta data: ",
                      paste(unique(sample.names[is.na(idx)]), collapse = ", "), "\n")
        msg1 = paste0("Make sure that sample names are set as paste0(meta$id, '_', meta$condition). \n")
        stop(c(msg0, msg1))
    }
    return(as.character(meta$condition)[idx])
}

.checkConditionNames <- function(condition) {
    # takes the condition column of the meta data
    # returns an error message if condition names are not valid or NULL
    # -> allowed are letters, numbers and one of the separators '-', '_' or '.'
    # -> the same separator has to be used in all condition names
    this.cond = unique(as.character(condition))
    msg = NULL

    invalidChar = this.cond[!grepl("^[A-Za-z0-9._-]+$", this.cond)]
    if (length(invalidChar) > 0) {
        msg = c(msg, paste0("Condition names can only contain letters, numbers and one of the separators '-', '_' or '.'. Invalid names: ",
                            paste(invalidChar, collapse = ", "), ". "))
    }
    seps = lapply(this.cond, function(x) {
        intersect(c("-", "_", "."), strsplit(x, "")[[1]])
    })
    mixedSep = this.cond[lengths(seps) > 1]
    if (length(mixedSep) > 0) {
        msg = c(msg, paste0("Condition names can only contain one type of separator ('-', '_' or '.'). Invalid names: ",
                            paste(mixedSep, collapse = ", "), ". "))
    }
    usedSep = unique(unlist(seps))
    if (length(mixedSep) == 0 & length(usedSep) > 1) {
        msg = c(msg, paste0("All condition names have to use the same separator. Found: ",
                            paste0("'", usedSep, "'", collapse = ", "), " in ",
                            paste(this.cond, collapse = ", "), ". "))
    }
    return(msg)
}
