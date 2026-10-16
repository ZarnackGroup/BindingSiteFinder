test_that("BSFDataSet() can be build correctly", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    rng = getRanges(bds)
    sgn = getSignal(bds)
    mta = getMeta(bds)
    # correct input
    bdsNew = expect_silent(BSFDataSet(ranges = rng, signal = sgn, meta = mta, silent = TRUE))
    bdsNew = expect_silent(BSFDataSet(ranges = rng, signal = sgn, meta = mta, silent = FALSE))
    bdsNew = expect_silent(BSFDataSet(ranges = rng, signal = sgn, meta = mta))

    # meta data table corrupt
    mtaWrong = mta
    mtaWrong$clPlus = rep(mta$clPlus[1],4)
    mtaWrong$clMinus = rep(mta$clMinus[1],4)
    expect_message(BSFDataSet(ranges = rng, signal = sgn, meta = mtaWrong))
    expect_silent(BSFDataSet(ranges = rng, signal = sgn, meta = mtaWrong, silent = TRUE))

    # check correct range sorting
    rng1 = .sortRanges(rng)
    expect_true(identical(rng, rng1))

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

    expect_message(setRanges(bds, rng2))
    expect_message(setRanges(bds, rng2, dropSeqlevels = TRUE))
    expect_message(expect_warning(setRanges(bds, rng2, dropSeqlevels = FALSE)))

})



test_that("condition names are checked for valid separators", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    setCond = function(cond) {
        m = getMeta(bds)
        m$condition = factor(rep(cond, each = 2), levels = cond)
        return(m)
    }

    # allowed: letters, numbers and one type of separator in all names
    for (cond in list(c("WT", "KO"), c("WT_1", "KO_15"), c("WT", "KO_15"),
                      c("WT-xx", "KO-7"), c("1.WT", "KO.15"))) {
        expect_no_error(setMeta(bds, setCond(cond)))
    }

    # two types of separators in one name
    expect_error(setMeta(bds, setCond(c("WT_1-2", "KO"))), "one type of separator")
    # different separators across names
    expect_error(setMeta(bds, setCond(c("WT_1", "KO-7"))), "same separator")
    # other characters
    expect_error(setMeta(bds, setCond(c("WT 1", "KO"))), "only contain letters")
    expect_error(setMeta(bds, setCond(c("WT/1", "KO"))), "only contain letters")

    # same check when building the object
    expect_error(BSFDataSet(ranges = getRanges(bds), meta = setCond(c("WT_1", "KO-7")),
                            signal = getSignal(bds)), "same separator")
})
