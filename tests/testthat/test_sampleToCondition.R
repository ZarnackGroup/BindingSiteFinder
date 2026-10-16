test_that(".sampleToCondition matches samples by id", {
    meta = data.frame(id = c("1", "2", "11", "12"),
                      condition = factor(c("WT_1", "WT_1", "KO-7", "KO-7")))

    # underscores in conditions and count prefixes
    expect_identical(.sampleToCondition(c("1_WT_1", "counts.bs.11_KO-7", "counts.bg.2_WT_1"), meta),
                     c("WT_1", "KO-7", "WT_1"))
    # id '1' is not confused with id '11' or '12'
    expect_identical(.sampleToCondition(c("11_KO-7", "1_WT_1", "12_KO-7"), meta),
                     c("KO-7", "WT_1", "KO-7"))
    # renamed conditions are taken from the meta data
    expect_identical(.sampleToCondition("counts.bs.1_WT", meta), "WT_1")
    # unknown ids give an error
    expect_error(.sampleToCondition("5_WT", meta))
})
