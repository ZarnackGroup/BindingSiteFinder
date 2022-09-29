test_that("Binding site merging works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    # testing input type errors
    expect_error(makeBindingSites(object = "bds", bsSize = 3, minWidth = 3,
                                  minCrosslinks = 2, minClSites = 1))
    expect_error(makeBindingSites(object = bds, bsSize = "3", minWidth = 3,
                                  minCrosslinks = 2, minClSites = 1))
    expect_error(makeBindingSites(object = bds, bsSize = 3, minWidth = "3",
                                  minCrosslinks = 2, minClSites = 1))
    expect_error(makeBindingSites(object = bds, bsSize = 3, minWidth = 3,
                                  minCrosslinks = "2", minClSites = 1))
    expect_error(makeBindingSites(object = bds, bsSize = 3, minWidth = 3,
                                  minCrosslinks = 2, minClSites = "1"))
    expect_error(makeBindingSites(object = bds, bsSize = 3, minWidth = 3,
                                  minCrosslinks = 2, minClSites = 4))
    expect_warning(expect_error(makeBindingSites(object = bds, bsSize = 3,
                                                 minWidth = 3,
                                  minCrosslinks = 20, minClSites = 1)))
    expect_warning(makeBindingSites(object = bds, bsSize = 3, minWidth = 4,
                                    minCrosslinks = 1, minClSites = 1))

    # testing output
    bds1 <- makeBindingSites(object = bds, bsSize = 3, minWidth = 2,
                             minCrosslinks = 2, minClSites = 1)
    expect_identical(class(bds), class(bds1))
    expect_identical(getMeta(bds), getMeta(bds1))
    expect_equal(length(unique(width(getRanges(bds1)))), 1)
    expect_equal(unique(width(getRanges(bds1))), 3)

    # test defaults
    bds1 <- makeBindingSites(object = bds, bsSize = 9, minWidth = 3,
                            minCrosslinks = 2, minClSites = 1)
    bds2 <- makeBindingSites(object = bds, bsSize = 9)
    expect_identical(bds1,bds2)

    bds1 <- makeBindingSites(object = bds, bsSize = 9, minWidth = 2,
                             minCrosslinks = 2, minClSites = 1)

})
