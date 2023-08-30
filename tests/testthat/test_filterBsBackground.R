test_that("filterBsBackground works", {
    # load clip data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))
    load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])

    # make testset
    bds = makeBindingSites(bds, bsSize = 7)
    bds = assignToGenes(bds, anno.genes = gns)
    bds = imputeBsDifferencesForTestdata(bds)
    bds = calculateBsBackground(bds, anno.genes = gns, use.offset = FALSE)

    # tests
    c0 = expect_warning(filterBsBackground(bds, flag = FALSE))
    expect_lt(length(getRanges(c0)), length(getRanges(bds)))

    c1 = expect_warning(filterBsBackground(bds, flag = TRUE))
    expect_equal(length(getRanges(c1)), length(getRanges(bds)))
    expect_equal(ncol(mcols(getRanges(c1))), ncol(mcols(getRanges(bds)))+3)

    c2 = expect_warning(filterBsBackground(bds, flag = TRUE, filter.minCounts = FALSE))
    c3 = expect_warning(filterBsBackground(bds, flag = TRUE, filter.ratio = FALSE))
    c4 = expect_warning(filterBsBackground(bds, flag = TRUE, filter.balance = FALSE))

    expect_equal(ncol(mcols(getRanges(bds)))+2, ncol(mcols(getRanges(c2))))
    expect_equal(ncol(mcols(getRanges(bds)))+2, ncol(mcols(getRanges(c3))))
    expect_equal(ncol(mcols(getRanges(bds)))+2, ncol(mcols(getRanges(c4))))

    c5 = filterBsBackground(bds, flag = TRUE, filter.balance = FALSE,
                            filter.ratio = FALSE, filter.minCounts = FALSE)
    expect_equal(ncol(mcols(getRanges(bds))), ncol(mcols(getRanges(c5))))
})
