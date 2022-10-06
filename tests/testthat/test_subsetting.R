test_that("Binding site merging works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    # test range of subset 1
    rng1 = getRanges(bds)[1:100]
    rng2 = getRanges(bds[1:100])
    expect_identical(rng1, rng2)

    # test range of subset 2
    set.seed(1234)
    samp = sample(c(1:length(getRanges(bds))), 10)
    rng1 = getRanges(bds)[samp]
    rng2 = getRanges(bds[samp])
    expect_identical(rng1, rng2)

    # check signal for correct subsetting - no drop
    set.seed(1234)
    samp = sample(c(1:length(getRanges(bds))), 10)
    rng1 = getRanges(bds)[samp]
    rng1 = keepSeqlevels(rng1, value = unique(seqnames(rng1)))
    rng1Seqlevels = seqlevels(rng1)
    sgn2 = getSignal(bds[samp])
    sgn2Seqlevels = unique(unlist(lapply(sgn2, function(x){lapply(x, function(y){names(y)})})))
    expect_identical(sort(rng1Seqlevels), sort(sgn2Seqlevels))
    # check signal for correct subsetting - with drop 1
    set.seed(1234)
    samp = sample(c(1:length(getRanges(bds))), 10)
    rng1 = getRanges(bds)[samp, drop = TRUE]
    rng1 = keepSeqlevels(rng1, value = unique(seqnames(rng1)))
    rng1Seqlevels = seqlevels(rng1)
    sgn2 = getSignal(bds[samp, drop = TRUE])
    sgn2Seqlevels = unique(unlist(lapply(sgn2, function(x){lapply(x, function(y){names(y)})})))
    expect_identical(sort(rng1Seqlevels), sort(sgn2Seqlevels))
    # check signal for correct subsetting - with drop 2
    sgn1 = getSignal(bds[samp])
    sgn2 = getSignal(bds[samp, drop = TRUE])
    expect_false(identical(sgn1, sgn2))
    sgn1 = getSignal(bds[samp])
    sgn2 = getSignal(bds[samp, drop = FALSE])
    expect_true(identical(sgn1, sgn2))
    sgn1 = getSignal(bds[samp])
    sgn2 = getSignal(bds[samp])
    expect_true(identical(sgn1, sgn2))
})
