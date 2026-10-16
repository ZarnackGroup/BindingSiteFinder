test_that("Reproducibility filter function works", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    expect_error(reproducibilityFilter(bds, cutoff = c(0.1,0.1)))
    expect_error(reproducibilityFilter(bds, nReps = c(2,2)))

    expect_is(reproducibilityFilter(bds, cutoff = c(0.1), nReps = c(2), returnType = "BSFDataSet"), "BSFDataSet")
    expect_is(reproducibilityFilter(bds, cutoff = c(0.1), nReps = c(2), returnType = "data.frame"), "data.frame")

    expect_error(reproducibilityFilter(bds, cutoff = c(0.1), nReps = c(1,1)))
    expect_error(reproducibilityFilter(bds, cutoff = c(0.1, 0.2), nReps = 1))

    expect_silent(reproducibilityFilter(bds, cutoff = c(0.1), nReps = 1))
    expect_error(expect_warning(reproducibilityFilter(bds, cutoff = c(0.1))))
    expect_error(expect_warning(reproducibilityFilter(bds, nReps = c(1))))

    expect_message(reproducibilityFilter(bds))

})

test_that("Reproducibility filter works with special characters in condition names", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))

    # one condition per sample, each with its own cutoff
    setConditions = function(object, cond) {
        m = getMeta(object)
        m$condition = factor(cond, levels = cond)
        s = getSignal(object)
        names(s$signalPlus) = paste0(m$id, "_", m$condition)
        names(s$signalMinus) = paste0(m$id, "_", m$condition)
        object = setMeta(object, m)
        object = setSignal(object, s)
        return(object)
    }
    bds = makeBindingSites(bds, bsSize = 7, quiet = TRUE)
    myCutoff = c(0.05, 0.1, 0.2, 0.3)

    b0 = setConditions(bds, c("A", "B", "C", "D"))
    r0 = reproducibilityFilter(b0, cutoff = myCutoff, nReps = c(1, 1, 1, 1), quiet = TRUE)

    b1 = setConditions(bds, c("WT_1", "WT_xx", "KO_15", "KO_7"))
    r1 = expect_no_error(reproducibilityFilter(b1, cutoff = myCutoff, nReps = c(1, 1, 1, 1), quiet = TRUE))
    # support columns are named after the conditions, so compare the values only
    expect_identical(granges(getRanges(r0)), granges(getRanges(r1)))
    expect_identical(unname(as.matrix(mcols(getRanges(r0))[, c("A", "B", "C", "D")])),
                     unname(as.matrix(mcols(getRanges(r1))[, c("WT_1", "WT_xx", "KO_15", "KO_7")])))
    # condition in the plot data matches the condition of each sample
    m1 = getMeta(b1)
    dfPlot = r1@plotData$reproducibilityFilterPlot$data
    expect_identical(dfPlot$condition,
                     as.character(m1$condition)[match(dfPlot$name, paste0(m1$id, "_", m1$condition))])
    expect_no_error(suppressWarnings(reproducibilityFilterPlot(r1)))
    expect_no_error(suppressWarnings(reproducibilityCutoffPlot(b1, cutoff = myCutoff)))
})
