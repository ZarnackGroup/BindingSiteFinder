test_that("targetGeneSpectrumPlot works with many gene types", {
    # load clip data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))
    load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])

    bds = makeBindingSites(bds, bsSize = 9, quiet = TRUE)
    bds = assignToGenes(bds, anno.genes = gns, quiet = TRUE)

    # add artificial gene types to exceed the size of the color palette
    tmp = bds@plotData$assignToGenes$dataSpectrum
    tmp = rbind(tmp, data.frame(geneType = paste0("Test", 1:15), nBs = 1:15, nGenes = 1:15))
    bds@plotData$assignToGenes$dataSpectrum = tmp

    for (n in c(5, 9, 20)) {
        expect_no_warning(p <- targetGeneSpectrumPlot(bds, showNGroups = n))
        fills = ggplot2::ggplot_build(p)$data[[1]]$fill
        expect_length(fills, n)
        expect_false(anyNA(fills))
        expect_equal(length(unique(fills)), n)
    }
})
