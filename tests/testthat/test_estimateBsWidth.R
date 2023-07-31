test_that("Test for estimateBsWidth()", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))
    # load genes
    load(list.files(files, pattern = "gns.rds$", full.names = TRUE))
    # load transcript region list
    load(list.files(files, pattern = "regions.rds$", full.names = TRUE))

    # Test for correct input
    expect_error(estimateBsWidth(object = bds, anno.genes = gns, est.maxBsWidth = 7, geneResolution = "medium", bsResolution = "medium"))
    expect_error(estimateBsWidth(object = bds, anno.genes = gns, est.maxBsWidth = 7, geneResolution = "medium", bsResolution = "medium",
                                   est.subsetChromosome = "chr21", quiet = TRUE))

    # Test for correct execution
    expect_error(estimateBsWidth(object = bds, anno.genes = gns, est.maxBsWidth = 7, geneResolution = "medium", bsResolution = "coarse",
                                   est.subsetChromosome = "chr22", quiet = TRUE, veryQuiet = TRUE))

    expect_output(estimateBsWidth(object = bds, anno.genes = gns, est.maxBsWidth = 13, geneResolution = "medium", bsResolution = "medium",
                                  est.subsetChromosome = "chr22", quiet = TRUE, veryQuiet = TRUE))

})
