test_that("Binding site merging works", {
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

    # testing input type errors
    expect_error(makeBindingSites(object = "bds", bsSize = 3, minWidth = 2,
                                  minCrosslinks = 2, minClSites = 1))
    expect_error(makeBindingSites(object = bds, bsSize = "3", minWidth = 2,
                                  minCrosslinks = 2, minClSites = 1))
    expect_error(makeBindingSites(object = bds, bsSize = 3, minWidth = "2",
                                  minCrosslinks = 2, minClSites = 1))
    expect_error(makeBindingSites(object = bds, bsSize = 3, minWidth = 2,
                                  minCrosslinks = "2", minClSites = 1))
    expect_error(makeBindingSites(object = bds, bsSize = 3, minWidth = 2,
                                  minCrosslinks = 2, minClSites = "1"))
    expect_error(makeBindingSites(object = bds, bsSize = 3, minWidth = 2,
                                  minCrosslinks = 2, minClSites = 4))
    expect_error(makeBindingSites(object = bds, bsSize = 3, minWidth = 2,
                                  minCrosslinks = 20, minClSites = 1))

    # testing output
    bds1 <- makeBindingSites(object = bds, bsSize = 3, minWidth = 2,
                             minCrosslinks = 2, minClSites = 1)
    expect_identical(class(bds), class(bds1))
    expect_identical(getMeta(bds), getMeta(bds1))
    expect_equal(length(unique(width(getRanges(bds1)))), 1)
    expect_equal(unique(width(getRanges(bds1))), 3)
})
