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

    # test if meta data is preserved in the ranges and gene columns are not duplicated
    # add meta data to the ranges
    r = getRanges(bds)
    r$test = "test"
    bds = setRanges(bds, r)

    for (ol in c("remove", "keep")) {
        res = suppressWarnings(calculateBsBackground(bds, anno.genes = gns, quiet = TRUE,
                                                     generate.geneID.bs = TRUE, overlaps = ol))
        cn = colnames(mcols(getRanges(res)))
        expect_false(anyDuplicated(cn) > 0)
        expect_true(all(c("test", "geneID", "geneName", "geneType") %in% cn))
        expect_true(all(getRanges(res)$test == "test"))
    }
})

test_that("calculateBsBackground works with special characters in condition names", {
    # load clip data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))
    load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])

    # set new condition names in meta data and signal
    setConditions = function(object, cond) {
        m = getMeta(object)
        m$condition = factor(rep(cond, each = 2), levels = cond)
        s = getSignal(object)
        names(s$signalPlus) = paste0(m$id, "_", m$condition)
        names(s$signalMinus) = paste0(m$id, "_", m$condition)
        object = setMeta(object, m)
        object = setSignal(object, s)
        return(object)
    }
    # get all counts as plain matrix without names
    getCounts = function(object) {
        df = as.data.frame(mcols(getRanges(object)))
        unname(as.matrix(df[, startsWith(colnames(df), "counts.")]))
    }

    bds = makeBindingSites(bds, bsSize = 7, quiet = TRUE)
    bds = assignToGenes(bds, anno.genes = gns, quiet = TRUE)
    set.seed(1234)
    myBlacklist = sample(getRanges(bds), size = 500) + 4

    c0 = calculateBsBackground(setConditions(bds, c("WT", "KO")), anno.genes = gns,
                               quiet = TRUE, blacklist = myBlacklist)

    # condition names that contain the names of internal count types
    for (cond in list(c("bs_KO", "gene_WT"), c("offset", "black"))) {
        c1 = calculateBsBackground(setConditions(bds, cond), anno.genes = gns,
                                   quiet = TRUE, blacklist = myBlacklist)
        expect_identical(getCounts(c0), getCounts(c1))
    }
})
