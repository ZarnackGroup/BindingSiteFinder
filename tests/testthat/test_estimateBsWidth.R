test_that("Test for estimateBsWidth()", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))
    # load genes
    load(list.files(files, pattern = "gns.rds$", full.names = TRUE))
    # load transcript region list
    load(list.files(files, pattern = "regions.rds$", full.names = TRUE))

    expect_error(estimateBsWidth(object = bds, anno.genes = gns, est.maxBsWidth = 7, geneResolution = "medium", bsResolution = "medium"))

    expect_warning(estimateBsWidth(object = bds, anno.genes = gns, est.maxBsWidth = 7, geneResolution = "medium", bsResolution = "medium",
                                   est.subsetChromosome = "chr22", quiet = TRUE))

    expect_output(estimateBsWidth(object = bds, anno.genes = gns, est.maxBsWidth = 13, geneResolution = "medium", bsResolution = "medium",
                                   est.subsetChromosome = "chr22", quiet = TRUE))

    # expect_silent(estimateBsWidth(object = bds, anno.genes = gns, est.maxBsWidth = 7, geneResolution = "coarse", bsResolution = "fine",
    #                               est.subsetChromosome = "chr22", veryQuiet = TRUE))
    # expect_silent(estimateBsWidth(object = bds, anno.genes = gns, est.maxBsWidth = 7, geneResolution = "coarse", bsResolution = "medium",
    #                               est.subsetChromosome = "chr22", veryQuiet = TRUE))
    # expect_silent(estimateBsWidth(object = bds, anno.genes = gns, est.maxBsWidth = 7, geneResolution = "coarse", bsResolution = "coarse",
    #                               est.subsetChromosome = "chr22", veryQuiet = TRUE))

})
