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

})
