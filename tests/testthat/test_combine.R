test_that("Combine function works", {
    # load clip data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    # make binding sites
    bds = makeBindingSites(bds, bsSize = 7)

    # split ranges in two groups
    allRanges = getRanges(bds)
    set.seed(1234)
    idx = sample(1:length(allRanges), size = length(allRanges)/2, replace = FALSE)
    r1 = allRanges[idx]
    r2 = allRanges[-idx]

    # splite meta data
    allMeta = getMeta(bds)
    m1 = allMeta[1:2,]
    m2 = allMeta[3:4,]

    # create new objects
    bds1 = setRanges(bds, r1)
    bds2 = setRanges(bds, r2)
    bds1 = setMeta(bds1, m1)
    bds2 = setMeta(bds2, m2)
    bds1 = setName(bds1, "test1")
    bds2 = setName(bds2, "test2")

    # merge two objects from list
    list = list(bds1, bds2)

    c1 = combineBSF(list = list, overlaps.fix = TRUE, combine.bsSize = NULL, combine.name = NULL, quiet = TRUE)
    r = unique(unique(width(getRanges(list[[1]])), width(getRanges(list[[2]]))))
    expect_identical(unique(width(getRanges(c1))), r)

    inSize = 9
    c2 = combineBSF(list = list, overlaps.fix = TRUE, combine.bsSize = inSize, combine.name = NULL, quiet = TRUE)
    expect_equal(unique(width(getRanges(c2))), inSize)

    c3 = combineBSF(list = list, overlaps.fix = FALSE, combine.bsSize = NULL, combine.name = NULL, quiet = TRUE)
    nr = length(c(getRanges(list[[1]]), getRanges(list[[2]])))
    expect_identical(length(getRanges(c3)), nr)

    inName = "test"
    c4 = combineBSF(list = list, overlaps.fix = FALSE, combine.bsSize = NULL, combine.name = inName, quiet = TRUE)
    expect_identical(getName(c4), inName)


    # tests with different width
    changeRanges = getRanges(bds1)
    bds1C = setRanges(bds1, changeRanges + 2)
    list2 = list(bds1C, bds2)

    expect_error(c1 = combineBSF(list = list2, overlaps.fix = TRUE, combine.bsSize = NULL, combine.name = NULL, quiet = TRUE))

    inSize = 9
    c2 = combineBSF(list = list2, overlaps.fix = TRUE, combine.bsSize = inSize, combine.name = NULL, quiet = TRUE)
    expect_equal(unique(width(getRanges(c2))), inSize)

    c3 = combineBSF(list = list2, overlaps.fix = FALSE, combine.bsSize = NULL, combine.name = NULL, quiet = TRUE)
    nr = length(c(getRanges(list2[[1]]), getRanges(list2[[2]])))
    expect_identical(length(getRanges(c3)), nr)
})


