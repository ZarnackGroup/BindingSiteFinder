if (.Platform$OS.type != "windows") {
    test_that("Test that score annotation function works", {

        csFile <- system.file("extdata", "PureCLIP_crosslink_sites_examples.bed",
                              package="BindingSiteFinder")
        cs = rtracklayer::import(con = csFile, format = "BED",
                                 extraCols=c("additionalScores" = "character"))
        cs$additionalScores = NULL
        clipFiles <- system.file("extdata", package="BindingSiteFinder")

        # use only if not windows
        meta = data.frame(
            id = c(1,2,3,4),
            condition = factor(c("WT", "WT", "KD", "KD"),
                               levels = c("KD", "WT")),
            clPlus = list.files(clipFiles,
                                pattern = "plus.bw$", full.names = TRUE),
            clMinus = list.files(clipFiles,
                                 pattern = "minus.bw$", full.names = TRUE))

        bdsOriginal = expect_warning(BSFDataSetFromBigWig(ranges = cs, meta = meta))

        # merge binding sites
        bdsMerge <- makeBindingSites(object = bdsOriginal, bsSize = 9,
                                     minWidth = 2,
                                     minCrosslinks = 2, minClSites = 1)

        expect_error(annotateWithScore(bdsOriginal, cs))
        expect_is(annotateWithScore(bdsMerge, cs), "BSFDataSet")
    })
}

