test_that("Test for assignToTranscriptRegions()", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))
    # load genes
    load(list.files(files, pattern = "gns.rds$", full.names = TRUE))
    # load transcript region list
    load(list.files(files, pattern = "regions.rds$", full.names = TRUE))


    # Prepare objects for test
    obj1 = pureClipGlobalFilter(object = bds, cutoff = 0.05)
    obj2_manual = pureClipGeneWiseFilter(object = obj1, anno.genes = gns, quiet = TRUE)
    # obj2_db = pureClipGeneWiseFilter(object = obj1, anno.annoDB = Homo.sapiens38)
    obj3_manual = makeBindingSites(object = obj2_manual, bsSize = 9, quiet = TRUE)
    # obj3_db = makeBindingSites(object = obj2_db, bsSize = 9)
    obj4_manual = assignToGenes(object = obj3_manual, anno.genes = gns, quiet = TRUE)
    # obj4_db = assignToGenes(object = obj3_db, anno.annoDB = Homo.sapiens38)

    # test with OrganismDB as annotation source
    # expect_message(assignToTranscriptRegions(object = obj4_db, anno.annoDB = Homo.sapiens38))
    # expect_message(assignToTranscriptRegions(object = obj4_db, anno.annoDB = Homo.sapiens38, overlappingLoci = "frequency"))
    # expect_message(assignToTranscriptRegions(object = obj4_db, anno.annoDB = Homo.sapiens38, overlappingLoci = "flag"))
    # expect_message(assignToTranscriptRegions(object = obj4_db, anno.annoDB = Homo.sapiens38, overlappingLoci = "remove"))
    #
    # expect_message(assignToTranscriptRegions(object = obj4_db, anno.annoDB = Homo.sapiens38,
    #                                          overlappingLoci = "hierarchy", rule = c("utr3", "cds", "utr5", "intron")))
    # expect_error(assignToTranscriptRegions(object = obj4_db, anno.annoDB = Homo.sapiens38, overlappingLoci = "hierarchy"))
    # expect_error(assignToTranscriptRegions(object = obj4_db, anno.annoDB = Homo.sapiens38,
    #                                        overlappingLoci = "hierarchy", rule = c("hamster", "maus")))

    # test with manual as annotation source
    expect_message(assignToTranscriptRegions(object = obj4_manual, anno.transcriptRegionList = regions))
    expect_message(assignToTranscriptRegions(object = obj4_manual, anno.transcriptRegionList = regions, overlaps = "frequency"))
    expect_message(assignToTranscriptRegions(object = obj4_manual, anno.transcriptRegionList = regions, overlaps = "flag"))
    expect_message(assignToTranscriptRegions(object = obj4_manual, anno.transcriptRegionList = regions, overlaps = "remove"))

    expect_message(assignToTranscriptRegions(object = obj4_manual, anno.transcriptRegionList = regions,
                                             overlaps = "hierarchy", overlaps.rule = c("utr3", "cds", "utr5", "intron")))
    expect_error(assignToTranscriptRegions(object = obj4_manual, anno.transcriptRegionList = regions, overlaps = "hierarchy"))
    expect_error(assignToTranscriptRegions(object = obj4_manual, anno.transcriptRegionList = regions,
                                           overlaps = "hierarchy", overlaps.rule = c("hamster", "maus")))

    x = assignToTranscriptRegions(object = obj4_manual, anno.transcriptRegionList = regions, overlaps = "frequency")
    expect_equal("frequency", x@params$assignToTranscriptRegions$overlaps)
    x = assignToTranscriptRegions(object = obj4_manual, anno.transcriptRegionList = regions, overlaps = "flag")
    expect_equal("flag", x@params$assignToTranscriptRegions$overlaps)

    # test quiet
    expect_silent(assignToTranscriptRegions(object = obj4_manual, anno.transcriptRegionList = regions, overlaps = "frequency", quiet = TRUE))
})
