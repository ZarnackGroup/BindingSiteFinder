test_that("Reproducibility filter function works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

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
