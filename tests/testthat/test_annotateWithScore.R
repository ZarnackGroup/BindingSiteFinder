test_that("Test that score annotation function works", {
    csFile <- system.file("extdata", "PureCLIP_crosslink_sites_example.bed",
                          package="BindingSiteFinder")
    cs = rtracklayer::import(con = csFile, format = "BED")
    clipFiles <- system.file("extdata", package="BindingSiteFinder")

    # two experimental conditions
    meta = data.frame(condition = factor(c("WT", "WT", "KD", "KD"), levels = c("KD", "WT")),
                      clPlus = list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE),
                      clMinus = list.files(clipFiles, pattern = "minus.bw$", full.names = TRUE))
    bdsOriginal = BSFDataSet(ranges = cs, meta = meta)

    # merge binding sites
    bdsMerge <- makeBindingSites(object = bdsOriginal, bsSize = 9, minWidth = 2,
                            minCrosslinks = 2, minClSites = 1)

    expect_error(annotateWithScore(bdsOriginal, cs))
    expect_is(annotateWithScore(bdsMerge, cs), "BSFDataSet")
})
