test_that("Reproducibility filter function works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))
    # check input vector length
    expect_error(reproducibilityFilter(bds, cutoff = rep(0.1, length(levels(testMeta$condition))+1)))
    expect_error(reproducibilityFilter(bds, n.reps = rep(0.1, length(levels(testMeta$condition))+1)))
    # check output data types
    expect_is(reproducibilityFilter(bds, returnType = "BSFDataSet"), "BSFDataSet")
    expect_is(reproducibilityFilter(bds, returnType = "data.frame"), "data.frame")

    expect_error(reproducibilityFilter(bds, cutoff = c(0.1), n.reps = c(1,1)))
    expect_error(reproducibilityFilter(bds, cutoff = c(0.1, 0.2), n.reps = 1))

    expect_silent(reproducibilityFilter(bds, cutoff = c(0.1), n.reps = c(1)))
    expect_silent(reproducibilityFilter(bds, cutoff = c(0.1, 0.1), n.reps = c(1,1)))

})

