test_that("BSFDataSet() can be build correctly", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    rng = getRanges(bds)
    sgn = getSignal(bds)
    mta = getMeta(bds)
    # correct input
    bdsNew = expect_message(BSFDataSet(ranges = rng, signal = sgn, meta = mta))

    # meta data table corrupt
    mtaWrong = mta
    mtaWrong$clPlus = rep(mta$clPlus[1],4)
    mtaWrong$clMinus = rep(mta$clMinus[1],4)
    expect_message(BSFDataSet(ranges = rng, signal = sgn, meta = mtaWrong))

    # check correct range sorting
    rng1 = .sortRanges(rng)
    expect_false(identical(rng, rng1))

    rngNew = getRanges(bdsNew)
    expect_silent(BSFDataSet(ranges = rngNew, signal = sgn, meta = mta))

    expect_identical(getRanges(bdsNew), .sortRanges(rngNew))

    # check for inconsistent chromosomes
    rng2 = rng %>% as.data.frame()
    rng2$seqnames[1:10] = "chr1"
    rng2 = makeGRangesFromDataFrame(rng2)

    expect_message(BSFDataSet(ranges = rng2, signal = sgn, meta = mta))
    expect_message(
        expect_warning(
            BSFDataSet(ranges = rng2, signal = sgn, meta = mta,
                       forceEqualNames = FALSE)))
})
