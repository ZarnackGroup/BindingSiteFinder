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
    expect_error(supportRatio(bds, bsWidths = "1"))
    expect_error(supportRatio(bds, bsWidths = c(1.1)))
    expect_error(supportRatio(bds, bsWidths = c(3,5), bsFlank = c(5,5,5)))
    expect_error(supportRatio(bds, bsWidths = c(3,5), bsFlank = c(5)))
    expect_error(supportRatio(bds, bsWidths = c(3,5), bsFlank = c(5,5.1)))
    # testing output
    expect_identical(supportRatio(bds, bsWidths = c(3,5)), supportRatio(bds, bsWidths = c(3,5), bsFlank = c(3,5)))

})
