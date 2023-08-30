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
    expect_equal(round(resNormal$bs.log2FoldChange,1), round(resChanged$bs.log2FoldChange,1) * -1)
    expect_equal(round(resNormal$bg.pvalue, 2), round(resChanged$bg.pvalue, 2))
    expect_equal(round(resNormal$bg.log2FoldChange,1), round(resChanged$bg.log2FoldChange,1) * -1)

    meta$condition = factor(c(rep("Mut",2), rep("Ctrl",2)), levels = c("Mut", "Ctrl"))
    bdsChange = setMeta(bds, meta)
    c3 = expect_silent(calculateBsFoldChange(bdsChange, veryQuiet = TRUE))
    resNormal = getRanges(c1)
    resChanged = getRanges(c3)
    expect_equal(resNormal$bs.pvalue, resChanged$bs.pvalue)
    expect_equal(resNormal$bs.log2FoldChange, resChanged$bs.log2FoldChange)
    expect_equal(resNormal$bg.pvalue, resChanged$bg.pvalue)
    expect_equal(resNormal$bg.log2FoldChange, resChanged$bg.log2FoldChange)
})
