test_that("Binding site merging works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    # testing input type errors
    expect_error(supportRatio(bds, bsWidths = "1"))
    expect_error(supportRatio(bds, bsWidths = c(1.1)))
    expect_error(supportRatio(bds, bsWidths = c(3,5), bsFlank = c(5,5,5)))
    expect_error(supportRatio(bds, bsWidths = c(3,5), bsFlank = c(5)))
    expect_error(supportRatio(bds, bsWidths = c(3,5), bsFlank = c(5,5.1)))

    # testing output
    expect_warning(supportRatio(bds, bsWidths = c(3,5)))
    expect_silent(supportRatio(bds, bsWidths = c(5,7)))
    expect_identical(supportRatio(bds, bsWidths = c(5,7)), supportRatio(bds, bsWidths = c(5,7), bsFlank = c(5,7)))
})
