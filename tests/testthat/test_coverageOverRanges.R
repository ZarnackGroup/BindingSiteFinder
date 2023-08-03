test_that("Coverage function works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    # Set artificial KD condition
    metaCond = getMeta(bds)
    metaCond$condition = factor(c(rep("WT", 2), rep("KD", 2)), levels = c("WT", "KD"))
    bdsCond = setMeta(bds, metaCond)
    # Fix replicate names in signal
    namesCond = c("1_WT", "2_WT", "3_KD", "4_KD")
    sgn = getSignal(bdsCond)
    names(sgn$signalPlus) = namesCond
    names(sgn$signalMinus) = namesCond
    bdsCond = setSignal(bdsCond, sgn)
    bds = bdsCond


    bds1 <- makeBindingSites(object = bds, bsSize = 5, minWidth = 2,
                              minCrosslinks = 2, minClSites = 1, quiet = TRUE)

    # 1) Return option results output with same sized input
    # --------------------------------------------------------------------------
    expect_warning(coverageOverRanges(
        bds1, returnOptions = "merge_ranges_keep_positions"))

    rngNew = getRanges(bds1)
    names(rngNew) = seq_along(rngNew)
    bds1 = setRanges(bds1, rngNew)
    # merge_ranges_keep_positions
    test = coverageOverRanges(
        bds1, returnOptions = "merge_ranges_keep_positions")
    expect_is(test,"matrix")
    expect_identical(nrow(test), nrow(getMeta(bds1)))
    expect_identical(ncol(test), unique(width(getRanges(bds1))))

    # merge_replicates_per_condition
    test = coverageOverRanges(
        bds1, returnOptions = "merge_replicates_per_condition")
    expect_is(test, "list")
    expect_identical(length(test), length(levels(getMeta(bds1)$condition)))


    expect_identical(dim(test[[1]]), dim(test[[2]]))

    # merge_all_replicates
    test = coverageOverRanges(
        bds1, returnOptions = "merge_all_replicates")
    expect_is(test, "matrix")
    expect_identical(nrow(test), length(getRanges(bds1)))
    expect_identical(ncol(test), unique(width(getRanges(bds1))))

    # merge_positions_keep_replicates
    test = coverageOverRanges(
        bds1, returnOptions = "merge_positions_keep_replicates")
    expect_is(test, "GRanges")
    expect_identical(length(test), length(getRanges(bds1)))
    expect_identical(ncol(mcols(test)), nrow(getMeta(bds1)))

    # 2) Return option results output with different sized input
    # --------------------------------------------------------------------------
    rngMod = getRanges(bds1)
    rngMod = resize(x = rngMod, fix = "center",
                    width = sample(1:10, size = length(rngMod),
                                   replace = TRUE))
    bdsMod = setRanges(bds1, rngMod)

    expect_warning(
        coverageOverRanges(
            bdsMod, returnOptions = "merge_ranges_keep_positions"))
    expect_warning(
        coverageOverRanges(
            bdsMod, returnOptions = "merge_replicates_per_condition"))
    expect_warning(
        coverageOverRanges(
            bdsMod, returnOptions = "merge_all_replicates"))


    test = coverageOverRanges(bdsMod,
                              returnOptions = "merge_replicates_per_condition",
                              allowNA = TRUE)
    expect_false(isTRUE(is.na(test$KD)))
    test = coverageOverRanges(bdsMod,
                              returnOptions = "merge_all_replicates",
                              allowNA = TRUE)
    expect_false(isTRUE(is.na(test)))
    test = coverageOverRanges(bdsMod,
                              returnOptions = "merge_ranges_keep_positions",
                              allowNA = TRUE)
    expect_false(isTRUE(is.na(test)))

    # 3) test sum/ mean method
    # --------------------------------------------------------------------------
    test = coverageOverRanges(bds1,
        returnOptions = "merge_positions_keep_replicates")
    test1 = coverageOverRanges(bds1, method = "sum",
                              returnOptions = "merge_positions_keep_replicates")
    expect_identical(test,test1)

    test = coverageOverRanges(bds1, method = "mean",
                              returnOptions = "merge_replicates_per_condition")
    test = coverageOverRanges(bds1, method = "mean",
                             returnOptions = "merge_all_replicates")
    # test = coverageOverRanges(bds1, method = "mean",
    #                          returnOptions = "merge_positions_keep_replicates")
    test = coverageOverRanges(bds1, method = "mean",
                             returnOptions = "merge_ranges_keep_positions")

})
