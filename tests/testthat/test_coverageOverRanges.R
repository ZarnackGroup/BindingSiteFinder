test_that("Coverage function works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    testRange = getRanges(bds)
    testMeta = getMeta(bds)

    # test correct behavoir of input options
    expect_is(coverageOverRanges(bds, merge = TRUE, returnType = "GRanges"), "GRanges")
    expect_is(coverageOverRanges(bds, merge = TRUE, returnType = "matrix"), "matrix")
    expect_is(coverageOverRanges(bds, merge = TRUE, returnType = "data.frame"), "data.frame")
    expect_is(coverageOverRanges(bds, merge = FALSE, returnType = "data.frame"), "data.frame")
    expect_warning(coverageOverRanges(bds, merge = FALSE, returnType = "GRanges"))
    expect_warning(coverageOverRanges(bds, merge = FALSE, returnType = "matrix"))

    # test correct format of output
    expect_equal(nrow(coverageOverRanges(bds, merge = TRUE, returnType = "data.frame")), length(testRange))
    expect_equal(ncol(coverageOverRanges(bds, merge = TRUE, returnType = "data.frame")), nrow(testMeta))
    expect_equal(ncol(coverageOverRanges(bds, merge = FALSE, returnType = "data.frame")), unique(width(testRange)))
})
