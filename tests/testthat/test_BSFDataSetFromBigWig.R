if (.Platform$OS.type != "windows") {
    test_that("BSFDataSetFromBigWig() can be build correctly", {
        # load data from which the BSFDataSet is constructed
        csFile <-
            system.file("extdata",
                        "PureCLIP_crosslink_sites_example.bed",
                        package = "BindingSiteFinder")
        cs = rtracklayer::import(con = csFile, format = "BED")
        clipFiles <- system.file("extdata", package = "BindingSiteFinder")
        meta = data.frame(
            id = c(1,2,3,4),
            condition = factor(c("WT", "WT", "KD", "KD"), levels = c("KD", "WT")),
            clPlus = list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE),
            clMinus = list.files(clipFiles, pattern = "minus.bw$", full.names = TRUE)
        )
        expect_warning(BSFDataSetFromBigWig(ranges = cs, meta = meta))

        testCs = cs
        strand(testCs) = c(rep("-", 500), rep("+", 500))
        expect_silent(BSFDataSetFromBigWig(
            ranges = testCs, meta = meta, silent = TRUE))

        testMeta = meta
        testMeta$id = c(1,2,3,3)
        expect_error(BSFDataSetFromBigWig(
            ranges = testCs, meta = testMeta, silent = TRUE))

        testMeta = meta[,c(1,2)]
        expect_error(BSFDataSetFromBigWig(
            ranges = testCs, meta = testMeta, silent = TRUE))

        testMeta = meta[c(1,2),]
        expect_silent(BSFDataSetFromBigWig(
            ranges = testCs, meta = testMeta, silent = TRUE))
    })
}
