test_that("calculateBsBackground works", {
    # load clip data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))
    load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])

    # make binding sites
    bds = makeBindingSites(bds, bsSize = 7)
    bds = assignToGenes(bds, anno.genes = gns)

    m = getMeta(bds)
    m$condition = factor(c("WT", "WT", "KO", "KO"), levels = c("WT", "KO"))

    s = getSignal(bds)
    names(s$signalPlus) = paste0(m$id, "_", m$condition)
    names(s$signalMinus) = paste0(m$id, "_", m$condition)

    bds = setMeta(bds, m)
    bds = setSignal(bds, s)

    # test normal mode
    expect_warning(calculateBsBackground(bds, anno.genes = gns))
    expect_silent(calculateBsBackground(bds, anno.genes = gns, quiet = TRUE))

    c0 = calculateBsBackground(bds, anno.genes = gns, quiet = TRUE)
    r0 = as.data.frame(mcols(getRanges(c0))) %>% select(contains("counts.bs"))
    r2 = as.data.frame(mcols(getRanges(c0))) %>% select(contains("counts.bg"))
    expect_identical(nrow(m), ncol(r0))
    expect_identical(nrow(m), ncol(r2))

    # test with blacklist
    myBlacklist = getRanges(bds)
    set.seed(1234)
    myBlacklist = sample(myBlacklist, size = 500) + 4

    c1 = calculateBsBackground(bds, anno.genes = gns, quiet = TRUE, blacklist = myBlacklist)
    r0 = as.data.frame(mcols(getRanges(c1))) %>% select(contains("counts.bs"))
    r2 = as.data.frame(mcols(getRanges(c1))) %>% select(contains("counts.bg"))
    expect_identical(nrow(m), ncol(r0))
    expect_identical(nrow(m), ncol(r2))

    # test with assignToGenes
    r = getRanges(bds)
    r$geneID = NULL
    bds = setRanges(bds, r)

    expect_error(calculateBsBackground(bds, anno.genes = gns, quiet = TRUE))

    c2 = calculateBsBackground(bds, anno.genes = gns, quiet = TRUE, generate.geneID.bs = TRUE)
    expect_identical(getRanges(c0)$geneID, getRanges(c2)$geneID)

    expect_warning(calculateBsBackground(bds, anno.genes = gns, quiet = TRUE, generate.geneID.bs = TRUE, overlaps = "remove"))
})
