test_that("filterBsBackground works", {
    # load clip data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))
    load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])

    # make testset
    bds = makeBindingSites(bds, bsSize = 7)
    bds = assignToGenes(bds, anno.genes = gns)
    bds = imputeBsDifferencesForTestdata(bds)
    bds = expect_warning(calculateBsBackground(bds, anno.genes = gns, use.offset = FALSE))

    # tests
    c0 = expect_warning(filterBsBackground(bds, flag = FALSE))
    expect_lt(length(getRanges(c0)), length(getRanges(bds)))

    c1 = expect_warning(filterBsBackground(bds, flag = TRUE))
    expect_equal(length(getRanges(c1)), length(getRanges(bds)))
    expect_equal(ncol(mcols(getRanges(c1))), ncol(mcols(getRanges(bds)))+3)

    c2 = expect_warning(filterBsBackground(bds, flag = TRUE, minCounts = FALSE))
    c3 = expect_warning(filterBsBackground(bds, flag = TRUE, balanceBackground = FALSE))
    c4 = expect_warning(filterBsBackground(bds, flag = TRUE, balanceCondition = FALSE))

    expect_equal(ncol(mcols(getRanges(bds)))+2, ncol(mcols(getRanges(c2))))
    expect_equal(ncol(mcols(getRanges(bds)))+2, ncol(mcols(getRanges(c3))))
    expect_equal(ncol(mcols(getRanges(bds)))+2, ncol(mcols(getRanges(c4))))

    c5 = filterBsBackground(bds, flag = TRUE, balanceCondition = FALSE,
                            balanceBackground = FALSE, minCounts = FALSE)
    expect_equal(ncol(mcols(getRanges(bds))), ncol(mcols(getRanges(c5))))
})


test_that("filterBsBackground works with special characters in condition names", {
    # load clip data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))
    load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])

    # set new condition names in meta data and signal
    setConditions = function(object, cond) {
        m = getMeta(object)
        m$condition = factor(rep(cond, each = 2), levels = cond)
        s = getSignal(object)
        names(s$signalPlus) = paste0(m$id, "_", m$condition)
        names(s$signalMinus) = paste0(m$id, "_", m$condition)
        object = setMeta(object, m)
        object = setSignal(object, s)
        return(object)
    }

    bds = makeBindingSites(bds, bsSize = 7, quiet = TRUE)
    bds = assignToGenes(bds, anno.genes = gns, quiet = TRUE)
    bds = imputeBsDifferencesForTestdata(bds)

    c0 = calculateBsBackground(bds, anno.genes = gns, use.offset = FALSE, quiet = TRUE)
    f0 = filterBsBackground(c0, quiet = TRUE, veryQuiet = TRUE)

    # underscores (see issue #22), hyphens and leading numbers
    for (cond in list(c("MYCN_WT", "MYCN_KO"), c("WT-1", "KO-7"), c("1_WT", "2_KO"))) {
        c1 = calculateBsBackground(setConditions(bds, cond), anno.genes = gns,
                                   use.offset = FALSE, quiet = TRUE)
        f1 = expect_no_error(filterBsBackground(c1, quiet = TRUE, veryQuiet = TRUE))
        expect_identical(getRanges(f0)$bsID, getRanges(f1)$bsID)
        expect_identical(colnames(f1@plotData$filterBsBackground$data.balanceCondition)[3:4], cond)
    }
})
