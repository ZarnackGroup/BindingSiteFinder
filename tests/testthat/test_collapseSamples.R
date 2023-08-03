test_that("Collapse samples", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    # Set artificial KD condition
    metaCond = getMeta(bds)
    metaCond$condition = factor(c(rep("WT", 2), rep("KD", 2)), levels = c("WT", "KD"))
    bdsCond = setMeta(bds, metaCond)
    # Fix replicate names in signal
    namesCond = c("1_WT", "2_WT", "3_KD", "4_KD")
    sgn = getSignal(bdsCond)
    names(sgn$signalPlus) = namesCond
    names(sgn$signalMinus) = namesCond
    bdsCond = setSignal(bdsCond, sgn)
    bds = bdsCond

    # test for collapseAll = default
    bdsNew = collapseReplicates(bds)
    expect_identical(getRanges(bds), getRanges(bdsNew))
    expect_identical(getMeta(bds), getMeta(bdsNew))

    cov1 = coverageOverRanges(bds, returnOptions = "merge_positions_keep_replicates")
    cov2 = coverageOverRanges(bdsNew, returnOptions = "merge_positions_keep_replicates")

    cov1Counts = mcols(cov1) %>% as.data.frame() %>%
        dplyr::mutate(WT = rowSums(dplyr::across(contains("WT")))) %>%
        dplyr::mutate(KD = rowSums(dplyr::across(contains("KD")))) %>%
        dplyr::select(WT, KD) %>%
        dplyr::relocate(WT, KD)
    cov2Counts = mcols(cov2) %>% as.data.frame()
    expect_identical(cov1Counts, cov2Counts)

    # test for collapseAll = FALSE
    bdsNew = collapseReplicates(bds, collapseAll = FALSE)
    expect_identical(getRanges(bds), getRanges(bdsNew))
    expect_identical(getMeta(bds), getMeta(bdsNew))

    cov1 = coverageOverRanges(bds, returnOptions = "merge_positions_keep_replicates")
    cov2 = coverageOverRanges(bdsNew, returnOptions = "merge_positions_keep_replicates")

    cov1Counts = mcols(cov1) %>% as.data.frame() %>%
        dplyr::mutate(WT = rowSums(dplyr::across(contains("WT")))) %>%
        dplyr::mutate(KD = rowSums(dplyr::across(contains("KD")))) %>%
        dplyr::select(WT, KD) %>%
        dplyr::relocate(WT, KD)
    cov2Counts = mcols(cov2) %>% as.data.frame()
    expect_identical(cov1Counts, cov2Counts)

    # test for collapseAll = TRUE
    bdsNew = collapseReplicates(bds, collapseAll = TRUE)
    expect_identical(getRanges(bds), getRanges(bdsNew))
    expect_identical(getMeta(bds), getMeta(bdsNew))

    cov1 = coverageOverRanges(bds, returnOptions = "merge_positions_keep_replicates")
    cov2 = coverageOverRanges(bdsNew, returnOptions = "merge_positions_keep_replicates")

    cov1Counts = mcols(cov1) %>% as.data.frame() %>%
        dplyr::mutate(All = rowSums(dplyr::across(everything()))) %>%
        dplyr::select(All)
    cov2Counts = mcols(cov2) %>% as.data.frame()
    expect_identical(cov1Counts, cov2Counts)

})

