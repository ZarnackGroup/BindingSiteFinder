test_that("Reproducibility filter function works", {
    # load data from which the BSFDataSet is constructed
    csFile <-
        system.file("extdata",
                    "PureCLIP_crosslink_sites_example.bed",
                    package = "BindingSiteFinder")
    cs = rtracklayer::import(con = csFile, format = "BED")
    clipFiles <- system.file("extdata", package = "BindingSiteFinder")
    meta = data.frame(
        condition = factor(c("WT", "WT", "KD", "KD"), levels = c("KD", "WT")),
        clPlus = list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE),
        clMinus = list.files(clipFiles, pattern = "minus.bw$", full.names = TRUE)
    )
    bds = expect_warning(BSFDataSet(ranges = cs, meta = meta))

    bds1 <- makeBindingSites(object = bds, bsSize = 9, minWidth = 2,
                             minCrosslinks = 2, minClSites = 1)

    testRange = getRanges(bds1)
    testMeta = getMeta(bds1)
    # check input vector length
    expect_error(reproducibilityFilter(bds, cutoff = rep(0.1, length(levels(testMeta$condition))+1)))
    expect_error(reproducibilityFilter(bds, n.reps = rep(0.1, length(levels(testMeta$condition))+1)))
    # check output data types
    expect_is(reproducibilityFilter(bds, returnType = "BSFDataSet"), "BSFDataSet")
    expect_is(reproducibilityFilter(bds, returnType = "data.frame"), "data.frame")
})
