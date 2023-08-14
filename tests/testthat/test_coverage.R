test_that("Coverage function works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    bds = makeBindingSites(bds, bsSize = 7)

    n.samples = nrow(getMeta(bds))
    n.conditions = length(levels(getMeta(bds)$condition))
    n.bs = length(getRanges(bds))
    bs.width = bds@params$bsSize

    # all identical
    t0 = clipCoverage(bds, out.format = "data.frame", positions.merge = TRUE, ranges.merge = TRUE, samples.merge = TRUE, samples.group = "all")
    expect_length(t0, 1)
    t05 = clipCoverage(bds, out.format = "data.frame", positions.merge = TRUE, ranges.merge = TRUE, samples.merge = TRUE, samples.group = "condition")
    expect_equal(nrow(t05), n.conditions)
    expect_equal(ncol(t05), 1)

    t1 = clipCoverage(bds, out.format = "data.frame", positions.merge = FALSE, ranges.merge = FALSE, samples.merge = FALSE)
    expect_equal(length(t1),n.samples)
    expect_equal(nrow(t1[[1]]),n.bs)
    expect_equal(ncol(t1[[1]]),bs.width)

    # single merge
    t2 = clipCoverage(bds, out.format = "data.frame", positions.merge = TRUE, ranges.merge = FALSE, samples.merge = FALSE)
    expect_equal(ncol(t2),bs.width)
    expect_equal(nrow(t2),n.samples)

    t3 = clipCoverage(bds, out.format = "data.frame", positions.merge = FALSE, ranges.merge = FALSE, samples.merge = TRUE, samples.group = "all")
    expect_equal(nrow(t3),n.bs)
    expect_equal(ncol(t3),bs.width)
    t35 = clipCoverage(bds, out.format = "data.frame", positions.merge = FALSE, ranges.merge = FALSE, samples.merge = TRUE, samples.group = "condition")
    expect_equal(length(t35), n.conditions)
    expect_equal(nrow(t35[[1]]), n.bs)
    expect_equal(ncol(t35[[1]]), bs.width)

    t4 = clipCoverage(bds, out.format = "data.frame", positions.merge = FALSE, ranges.merge = TRUE, samples.merge = FALSE)
    expect_equal(length(t4), n.samples)
    expect_equal(length(t4[[1]]), n.bs)

    # double merge
    t5 = clipCoverage(bds, out.format = "data.frame", positions.merge = TRUE, ranges.merge = FALSE, samples.merge = TRUE, samples.group = "all")
    expect_equal(length(t5), bs.width)
    t55 = clipCoverage(bds, out.format = "data.frame", positions.merge = TRUE, ranges.merge = FALSE, samples.merge = TRUE, samples.group = "condition")
    expect_equal(nrow(t55), n.conditions)
    expect_equal(ncol(t55), bs.width)

    t6 = clipCoverage(bds, out.format = "data.frame", positions.merge = FALSE, ranges.merge = TRUE, samples.merge = TRUE, samples.group = "all")
    expect_equal(length(t6), n.bs)
    t65 = clipCoverage(bds, out.format = "data.frame", positions.merge = FALSE, ranges.merge = TRUE, samples.merge = TRUE, samples.group = "condition")
    expect_equal(length(t65), n.conditions)
    expect_equal(length(t65[[1]]), n.bs)

    t7 = clipCoverage(bds, out.format = "data.frame", positions.merge = TRUE, ranges.merge = TRUE, samples.merge = FALSE)
    expect_equal(nrow(t7), n.samples)


})
