test_that("Test that score annotation function works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    # prepare objects
    obj0 = bds
    obj1 = makeBindingSites(object = obj0, bsSize = 9, quiet = TRUE)

    # test with OrganismDB as annotation source
    expect_error(annotateWithScore(object = obj1))
    expect_error(annotateWithScore(object = obj1, match.ranges = obj0))
    expect_silent(annotateWithScore(object = obj1, match.ranges = getRanges(obj0)))
    expect_equal(annotateWithScore(object = obj1, match.ranges = getRanges(obj0)),
                 annotateWithScore(object = obj1, match.ranges = getRanges(obj0), match.option = "max"))

    expect_silent(annotateWithScore(object = obj1, match.ranges = getRanges(obj0), match.option = "sum"))
    expect_silent(annotateWithScore(object = obj1, match.ranges = getRanges(obj0), match.option = "mean"))

    expect_error(annotateWithScore(object = obj1, match.ranges = getRanges(obj0), match.option = "test"))

})
