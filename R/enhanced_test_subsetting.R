library(testthat)
test_that("Binding site merging works", {
    # load data
    TEST_DATASETS_SUBSET = readRDS(file = "/Users/mirko/Projects/PackageDevelopment/BindingSiteDefinition/TestEnvironment/data/TEST_DATASETS_SUBSET.rds")

    # test range of subset 1
    sapply(TEST_DATASETS_SUBSET, function(x){
        rng1 = getRanges(x)[1:100]
        rng2 = getRanges(x[1:100])
        expect_identical(rng1, rng2)
    })
    # test range of subset 2
    sapply(TEST_DATASETS_SUBSET, function(x){
        set.seed(1234)
        samp = sample(c(1:length(getRanges(x))), 10)
        rng1 = getRanges(x)[samp]
        rng2 = getRanges(x[samp])
        expect_identical(rng1, rng2)
    })
    # test subset on coverage 1
    sapply(TEST_DATASETS_SUBSET, function(x){
        cov1 = coverageOverRanges(x, returnOptions = "merge_positions_keep_replicates")
        cov1 = cov1[c(1:100),]
        cov2 = coverageOverRanges(x[1:100], returnOptions = "merge_positions_keep_replicates")
        expect_identical(cov1, cov2)
    })
    # test subset on coverage 2
    sapply(TEST_DATASETS_SUBSET, function(x){
        set.seed(1234)
        samp = sample(c(1:length(getRanges(x))), 10)
        x  = TEST_DATASETS_SUBSET$Package_Test
        cov1 = coverageOverRanges(x, returnOptions = "merge_positions_keep_replicates")
        cov1 = cov1[samp,]
        cov2 = coverageOverRanges(x[samp], returnOptions = "merge_positions_keep_replicates")
        expect_identical(cov1, cov2)
    })

})
