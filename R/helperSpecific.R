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
        silent = TRUE))))

    objMod = object
    objMod = suppressMessages(setRanges(objMod, getRanges(object) + flankSize))
    c1 = rowSums(as.data.frame(mcols(coverageOverRanges(
        objMod, returnOptions = "merge_positions_keep_replicates",
        silent = TRUE))))

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
    # returns the inital binding site ranges with matched gene info
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
    idx = match(rng$currIdx, olsInfo$bsIndex)
    # handle intergenic cases
    rngIntergenic = rng[which(is.na(idx))]
    rngIntergenic$geneID = NA
    rngIntergenic$geneName = NA
    rngIntergenic$geneType = "Intergenic"

    # handle normal cases
    # rngRest = rng[!rng %in% rngIntergenic]
    rngRest = rng[countOverlaps(rng, rngIntergenic) == 0]

    rngRest$geneID[idx[!is.na(idx)]] = olsInfo$geneID
    rngRest$geneName[idx[!is.na(idx)]] = olsInfo$geneName
    rngRest$geneType[idx[!is.na(idx)]] = olsInfo$geneType
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

