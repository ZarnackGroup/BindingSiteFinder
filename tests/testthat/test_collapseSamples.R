test_that("Binding site merging works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

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
        dplyr::relocate(KD, WT)
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
        dplyr::relocate(KD, WT)
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

