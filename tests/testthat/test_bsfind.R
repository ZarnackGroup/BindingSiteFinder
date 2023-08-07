test_that("BSFind wrapper works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))
    # Load genes
    load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
    # load transcript regions
    load(list.files(files, pattern = ".rds$", full.names = TRUE)[2])


    expect_message(BSFind(object = bds, bsSize = 9, anno.genes = gns,
                          anno.transcriptRegionList = regions,
                          est.subsetChromosome = "chr22",
                          quiet = FALSE, veryQuiet = FALSE))

    expect_silent(BSFind(object = bds, bsSize = 9, anno.genes = gns,
                          anno.transcriptRegionList = regions,
                          est.subsetChromosome = "chr22",
                          quiet = TRUE, veryQuiet = TRUE))

    # ------
    my.size = 13
    my.filter = 0.1

    this.bds = BSFind(object = bds, bsSize = my.size, anno.genes = gns,
                     anno.transcriptRegionList = regions,
                     est.subsetChromosome = "chr22",
                     quiet = TRUE, veryQuiet = TRUE)
    expect_identical(my.size, this.bds@params$bsSize)
    expect_false(my.filter == this.bds@params$geneFilter)

    this.bds = BSFind(object = bds, cutoff.geneWiseFilter = my.filter, anno.genes = gns,
                      anno.transcriptRegionList = regions,
                      est.subsetChromosome = "chr22",
                      quiet = TRUE, veryQuiet = TRUE)
    expect_false(my.size == this.bds@params$bsSize)
    expect_identical(my.filter, this.bds@params$geneFilter)

    this.bds = BSFind(object = bds,
                      bsSize = my.size,
                      cutoff.geneWiseFilter = my.filter,
                      anno.genes = gns,
                      anno.transcriptRegionList = regions,
                      est.subsetChromosome = "chr22",
                      quiet = TRUE, veryQuiet = TRUE)
    expect_identical(my.size, this.bds@params$bsSize)

    expect_silent(BSFind(object = bds,
                      bsSize = my.size,
                      cutoff.geneWiseFilter = my.filter,
                      anno.genes = gns,
                      anno.transcriptRegionList = regions,
                      est.subsetChromosome = "chr22",
                      quiet = TRUE, veryQuiet = TRUE))

})
