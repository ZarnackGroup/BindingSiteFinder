test_that("Reproducibility filter function works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    expect_error(reproducibilityFilter(bds, cutoff = 0.1))
    expect_error(reproducibilityFilter(bds, n.reps = 2))

    expect_is(reproducibilityFilter(bds, cutoff = c(0.1, 0.1), n.reps = c(2,2), returnType = "BSFDataSet"), "BSFDataSet")
    expect_is(reproducibilityFilter(bds, cutoff = c(0.1, 0.1), n.reps = c(2,2), returnType = "data.frame"), "data.frame")

    expect_error(reproducibilityFilter(bds, cutoff = c(0.1), n.reps = c(1,1)))
    expect_error(reproducibilityFilter(bds, cutoff = c(0.1, 0.2), n.reps = 1))

    expect_silent(reproducibilityFilter(bds, cutoff = c(0.1, 0.1), n.reps = c(1,1)))

})

