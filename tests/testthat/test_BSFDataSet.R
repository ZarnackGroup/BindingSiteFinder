test_that("BSFDataSet() can be build correctly", {
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

    # testing input ranges
    testCs = cs
    expect_warning(BSFDataSet(ranges = testCs, meta = meta))
    strand(testCs) = "*"
    expect_error(expect_warning(BSFDataSet(ranges = testCs, meta = meta)))
    strand(testCs) = "+"
    expect_warning(BSFDataSet(ranges = testCs, meta = meta))
    strand(testCs) = "-"
    expect_warning(BSFDataSet(ranges = testCs, meta = meta))
    strand(testCs) = c(rep("-", 500), rep("+", 500))
    BSFDataSet(ranges = testCs, meta = meta)

    # testing input meta data frame
    testMeta = meta
    expect_silent(BSFDataSet(ranges = testCs, meta = testMeta))
    colnames(testMeta)[1] = "a"
    expect_error(BSFDataSet(ranges = testCs, meta = testMeta))
    colnames(testMeta) = c("a", "b", "c")
    expect_error(BSFDataSet(ranges = testCs, meta = testMeta))
})

