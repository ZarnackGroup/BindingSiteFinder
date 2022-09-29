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

    # test force equal seqnames
    rng = getRanges(bds)
    sgn = getSignal(bds)

    r = as.data.frame(granges(rng))
    d = data.frame(seqnames = "chr6", start = 1, end = 1, width = 1, strand = "+")
    rng2 = rbind(d, r) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

    expect_error(bds2 = setRanges(bds, rng2))

    fixed = forceEqualNames(rng2, sgn)
    expect_silent(setRanges(bds, fixed$ranges))
    expect_silent(setSignal(bds, fixed$signal))

})




