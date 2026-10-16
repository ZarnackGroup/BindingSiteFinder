test_that("calculateBsFoldChange works", {
    # load clip data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))
    load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])

    # make example data
    bds = makeBindingSites(bds, bsSize = 7)
    bds = assignToGenes(bds, anno.genes = gns)
    bds = imputeBsDifferencesForTestdata(bds)
    bds = expect_warning(calculateBsBackground(bds, anno.genes = gns, use.offset = FALSE))

    c0 = expect_silent(calculateBsFoldChange(bds, veryQuiet = TRUE))

    expect_identical(length(getRanges(c0)), length(getRanges(bds)))

    rngMod = getRanges(bds)
    bdsNoCounts = setRanges(bds, granges(rngMod))
    expect_error(calculateBsFoldChange(bdsNoCounts))

    meta = getMeta(bds)
    meta$condition = factor(c(rep("WT",4)), levels = c("WT"))
    bdsChange = setMeta(bds, meta)
    expect_error(calculateBsFoldChange(bdsChange, veryQuiet = TRUE))

    meta = getMeta(bds)
    meta$condition = factor(c(rep("Ctrl",2), rep("Mut",2)), levels = c("Ctrl", "Mut"))
    bdsChange = setMeta(bds, meta)
    c1 = expect_silent(calculateBsFoldChange(bdsChange, veryQuiet = TRUE))
    expect_identical(getRanges(c0), getRanges(c1))

    meta$condition = factor(c(rep("Ctrl",2), rep("Mut",2)), levels = c("Mut", "Ctrl"))
    bdsChange = setMeta(bds, meta)
    c2 = expect_silent(calculateBsFoldChange(bdsChange, veryQuiet = TRUE))
    resNormal = getRanges(c1)
    resChanged = getRanges(c2)
    expect_equal(resNormal$bs.pvalue, resChanged$bs.pvalue)
    # expect_equal(round(resNormal$bs.log2FoldChange,1), round(resChanged$bs.log2FoldChange,1) * -1)
    # expect_equal(round(resNormal$bg.pvalue, 2), round(resChanged$bg.pvalue, 2))
    # expect_equal(round(resNormal$bg.log2FoldChange,1), round(resChanged$bg.log2FoldChange,1) * -1)

    meta$condition = factor(c(rep("Mut",2), rep("Ctrl",2)), levels = c("Mut", "Ctrl"))
    bdsChange = setMeta(bds, meta)
    c3 = expect_silent(calculateBsFoldChange(bdsChange, veryQuiet = TRUE))
    resNormal = getRanges(c1)
    resChanged = getRanges(c3)
    expect_equal(resNormal$bs.pvalue, resChanged$bs.pvalue)
    expect_equal(resNormal$bs.log2FoldChange, resChanged$bs.log2FoldChange)
    expect_equal(resNormal$bg.pvalue, resChanged$bg.pvalue)
    expect_equal(resNormal$bg.log2FoldChange, resChanged$bg.log2FoldChange)

    # negative counts
    rng = getRanges(bds)
    rng$counts.bs.1_WT[1] = -10
    bdsNeg = setRanges(bds, rng)
    expect_error(calculateBsFoldChange(bdsNeg))
    expect_error(calculateBsFoldChange(bdsNeg, removeNA = TRUE))
    expect_warning(calculateBsFoldChange(bdsNeg, replaceNegative = TRUE))

    # NA counts
    rng = getRanges(bds)
    rng$counts.bs.1_WT[1] = NA
    bdsNA = setRanges(bds, rng)
    expect_error(calculateBsFoldChange(bdsNA))
    expect_error(calculateBsFoldChange(bdsNA, replaceNegative = TRUE))
    expect_warning(calculateBsFoldChange(bdsNA, removeNA = TRUE))
})

test_that("calculateBsFoldChange works with special characters in condition names", {
    # load clip data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))
    load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])

    # set new condition names in meta data and signal
    setConditions = function(object, cond) {
        m = getMeta(object)
        m$condition = factor(rep(cond, each = 2), levels = cond)
        s = getSignal(object)
        names(s$signalPlus) = paste0(m$id, "_", m$condition)
        names(s$signalMinus) = paste0(m$id, "_", m$condition)
        object = setMeta(object, m)
        object = setSignal(object, s)
        return(object)
    }

    bds = makeBindingSites(bds, bsSize = 7, quiet = TRUE)
    bds = assignToGenes(bds, anno.genes = gns, quiet = TRUE)
    bds = imputeBsDifferencesForTestdata(bds)

    c0 = calculateBsBackground(bds, anno.genes = gns, use.offset = FALSE, quiet = TRUE)
    r0 = getRanges(calculateBsFoldChange(c0, quiet = TRUE, veryQuiet = TRUE))

    # hyphens, dots and leading numbers
    # -> DESeq2 converts these to syntactically valid names internally
    for (cond in list(c("WT-xx", "KO-7"), c("1.WT", "KO.15"))) {
        c1 = calculateBsBackground(setConditions(bds, cond), anno.genes = gns,
                                   use.offset = FALSE, quiet = TRUE)
        r1 = getRanges(calculateBsFoldChange(c1, quiet = TRUE, veryQuiet = TRUE))
        expect_identical(r0$bs.log2FoldChange, r1$bs.log2FoldChange)
        expect_identical(r0$bg.log2FoldChange, r1$bg.log2FoldChange)
    }
})
