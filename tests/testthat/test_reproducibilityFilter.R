test_that("Reproducibility filter function works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    expect_error(reproducibilityFilter(bds, cutoff = c(0.1,0.1)))
    expect_error(reproducibilityFilter(bds, nReps = c(2,2)))

    expect_is(reproducibilityFilter(bds, cutoff = c(0.1), nReps = c(2), returnType = "BSFDataSet"), "BSFDataSet")
    expect_is(reproducibilityFilter(bds, cutoff = c(0.1), nReps = c(2), returnType = "data.frame"), "data.frame")

    expect_error(reproducibilityFilter(bds, cutoff = c(0.1), nReps = c(1,1)))
    expect_error(reproducibilityFilter(bds, cutoff = c(0.1, 0.2), nReps = 1))

    expect_silent(reproducibilityFilter(bds, cutoff = c(0.1), nReps = 1))
    expect_error(expect_warning(reproducibilityFilter(bds, cutoff = c(0.1))))
    expect_error(expect_warning(reproducibilityFilter(bds, nReps = c(1))))

    expect_message(reproducibilityFilter(bds))

})
