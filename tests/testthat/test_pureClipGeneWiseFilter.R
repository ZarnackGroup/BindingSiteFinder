test_that("Test for pureClipGeneWiseFilter()", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))
    # load genes
    load(list.files(files, pattern = "gns.rds$", full.names = TRUE))

    # tests with anno.genes as annotion resources
    expect_message(pureClipGeneWiseFilter(object = bds, anno.genes = gns, cutoff = 0.2, overlaps = "keepSingle"))
    expect_warning(pureClipGeneWiseFilter(object = bds, anno.genes = gns, cutoff = 0.2, overlaps = "removeAll"))
    expect_warning(pureClipGeneWiseFilter(object = bds, anno.genes = gns, cutoff = 0.2, overlaps = "keepAll"))

    expect_message(pureClipGeneWiseFilter(object = bds, anno.genes = gns, overlaps = "keepSingle"))
    expect_warning(pureClipGeneWiseFilter(object = bds, anno.genes = gns, overlaps = "removeAll"))
    expect_warning(pureClipGeneWiseFilter(object = bds, anno.genes = gns, overlaps = "keepAll"))

    expect_message(pureClipGeneWiseFilter(object = bds, anno.genes = gns))

    expect_error(pureClipGeneWiseFilter(object = bds))

    # tests with anno.annoDB as annotaiton resource
    # expect_message(pureClipGeneWiseFilter(object = bds, anno.annoDB = Homo.sapiens38, cutoff = 0.2, overlaps = "keepSingle"))
    # expect_message(expect_warning(pureClipGeneWiseFilter(object = bds, anno.annoDB = Homo.sapiens38, cutoff = 0.2, overlaps = "removeAll")))
    # expect_message(expect_warning(pureClipGeneWiseFilter(object = bds, anno.annoDB = Homo.sapiens38, cutoff = 0.2, overlaps = "keepAll")))
    #
    # expect_message(pureClipGeneWiseFilter(object = bds, anno.annoDB = Homo.sapiens38, overlaps = "keepSingle"))
    # expect_message(expect_warning(pureClipGeneWiseFilter(object = bds, anno.annoDB = Homo.sapiens38, overlaps = "removeAll")))
    # expect_message(expect_warning(pureClipGeneWiseFilter(object = bds, anno.annoDB = Homo.sapiens38, overlaps = "keepAll")))
    #
    # expect_message(pureClipGeneWiseFilter(object = bds, anno.annoDB = Homo.sapiens38))
    #
    # expect_error(pureClipGeneWiseFilter(object = bds, anno.genes = gns, anno.annoDB = Homo.sapiens38))

    # check cutoff 0
    testObj1 = expect_message(pureClipGeneWiseFilter(object = bds, anno.genes = gns, cutoff = 0, overlaps = "keepSingle"))
    expect_identical(length(getRanges(testObj1)), length(getRanges(bds)))

    # check quiet mode
    expect_silent(pureClipGeneWiseFilter(object = bds, anno.genes = gns, cutoff = 0.1, quiet = TRUE))

})

