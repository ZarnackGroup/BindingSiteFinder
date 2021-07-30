
test_that("BSFDataSet() can be build correctly", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    rng = getRanges(bds)
    sgn = getSignal(bds)
    mta = getMeta(bds)

    bdsNew = BSFDataSet(ranges = rng, signal = sgn, meta = mta)


    # BSFDataSet(ranges = cs, meta = meta, signal = exampleSignal)
    #
    #
    #
    # exampleSignal = getSignal(bds)
    #
    # meta = data.frame(id = c(1,2,3,4),
    #                   condition = factor(c("WT", "WT", "KD", "KD"),
    #                                      levels = c("KD", "WT")))
    #
    # x = BSFDataSet(ranges = cs, meta = meta, signal = exampleSignal)
    #
    # testSignal = exampleSignal
    # testSignal$signalPlus = NULL
    # x = BSFDataSet(ranges = cs, meta = meta, signal = testSignal)
    #
    # testSignal = exampleSignal
    # # names(testSignal) = c("signalPlus", "haha")
    # x = BSFDataSet(ranges = cs, meta = meta, signal = testSignal, forceEqualNames = F)
    #
    #
    # # testing input ranges
    # testCs = cs
    # expect_warning(BSFDataSet(ranges = testCs, meta = meta))
    # strand(testCs) = "*"
    # expect_error(expect_warning(BSFDataSet(ranges = testCs, meta = meta)))
    # strand(testCs) = "+"
    # expect_warning(BSFDataSet(ranges = testCs, meta = meta))
    # strand(testCs) = "-"
    # expect_warning(BSFDataSet(ranges = testCs, meta = meta))
    # strand(testCs) = c(rep("-", 500), rep("+", 500))
    # BSFDataSet(ranges = testCs, meta = meta)
    #
    # # testing input meta data frame
    # testMeta = meta
    # expect_silent(BSFDataSet(ranges = testCs, meta = testMeta, silent = TRUE))
    # colnames(testMeta)[1] = "a"
    # expect_error(BSFDataSet(ranges = testCs, meta = testMeta))
    # colnames(testMeta) = c("a", "b", "c")
    # expect_error(BSFDataSet(ranges = testCs, meta = testMeta))
})





