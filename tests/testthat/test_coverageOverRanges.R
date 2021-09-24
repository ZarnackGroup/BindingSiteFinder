test_that("Coverage function works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    bds1 <- makeBindingSites(object = bds, bsSize = 5, minWidth = 2,
                              minCrosslinks = 2, minClSites = 1)

    expect_warning(coverageOverRanges(
        bds1, returnOptions = "merge_ranges_keep_positions"))

    rngNew = getRanges(bds1)
    names(rngNew) = seq_along(rngNew)
    bds1 = setRanges(bds1, rngNew)
    # merge_ranges_keep_positions
    test = coverageOverRanges(
        bds1, returnOptions = "merge_ranges_keep_positions")
    expect_is(test,"matrix")
    expect_identical(nrow(test), nrow(getMeta(bds1)))
    expect_identical(ncol(test), unique(width(getRanges(bds1))))

    # merge_replicates_per_condition
    test = coverageOverRanges(
        bds1, returnOptions = "merge_replicates_per_condition")
    expect_is(test, "list")
    expect_identical(length(test), length(levels(getMeta(bds1)$condition)))
    expect_identical(dim(test[[1]]), dim(test[[2]]))

    # merge_all_replicates
    test = coverageOverRanges(
        bds1, returnOptions = "merge_all_replicates")
    expect_is(test, "matrix")
    expect_identical(nrow(test), length(getRanges(bds1)))
    expect_identical(ncol(test), unique(width(getRanges(bds1))))

    # merge_positions_keep_replicates
    test = coverageOverRanges(
        bds1, returnOptions = "merge_positions_keep_replicates")
    expect_is(test, "GRanges")
    expect_identical(length(test), length(getRanges(bds1)))
    expect_identical(ncol(mcols(test)), nrow(getMeta(bds1)))

    # test different width ranges
    rngMod = getRanges(bds1)
    rngMod = resize(x = rngMod, fix = "center",
                    width = sample(1:10, size = length(rngMod),
                                   replace = TRUE))
    bdsMod = setRanges(bds1, rngMod)
    expect_warning(
        coverageOverRanges(bdsMod,
                           returnOptions = "merge_positions_keep_replicates"))


})


