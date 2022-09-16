test_that("Binding site merging works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    # test range of subset 1
    rng1 = getRanges(bds)[1:100]
    rng2 = getRanges(bds[1:100])
    expect_identical(rng1, rng2)

    # test range of subset 2
    set.seed(1234)
    samp = sample(c(1:length(getRanges(bds))), 10)
    rng1 = getRanges(bds)[samp]
    rng2 = getRanges(bds[samp])
    expect_identical(rng1, rng2)
})
