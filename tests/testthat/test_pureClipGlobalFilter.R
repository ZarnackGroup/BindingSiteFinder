test_that("Test for pureClipGlobalFilter()", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    # tests with genes as annotion resources
    expect_silent(pureClipGlobalFilter(object = bds, match.score = "score"))
    expect_silent(pureClipGlobalFilter(object = bds, cutoff = 0.1, match.score = "score"))
    expect_error(pureClipGlobalFilter(object = bds, cutoff = 0.1, match.score = "score2"))

    rng = getRanges(bds)
    mcols(rng)$newScore = rng$score
    mcols(rng)$score = NULL
    objNew = setRanges(bds, rng)

    expect_error(pureClipGlobalFilter(object = objNew, cutoff = 0.1, match.score = "score"))
    expect_silent(pureClipGlobalFilter(object = objNew, cutoff = 0.1, match.score = "newScore"))
})
