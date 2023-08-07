test_that("Test for assignToGenes()", {
    # load data
    files <- system.file("extdata", package="BindingSiteFinder")
    load(list.files(files, pattern = ".rda$", full.names = TRUE))
    # load genes
    load(list.files(files, pattern = "gns.rds$", full.names = TRUE))

    # Prepare objects for test
    obj1 = pureClipGlobalFilter(object = bds, cutoff = 0.05)
    obj2_manual = pureClipGeneWiseFilter(object = obj1, anno.genes = gns, quiet = TRUE)
    # obj2_db = pureClipGeneWiseFilter(object = obj1, anno.annoDB = Homo.sapiens38)
    obj3_manual = makeBindingSites(object = obj2_manual, bsSize = 9)
    # obj3_db = makeBindingSites(object = obj2_db, bsSize = 9)

    # Prepare rule
    selectTerms = c("protein_coding", "miRNA", "tRNA", "lincRNA", "snRNA")
    GENE_RULE = unique(gns$gene_type)
    GENE_RULE = GENE_RULE[!GENE_RULE %in% selectTerms]
    GENE_RULE = c(selectTerms, GENE_RULE)
    rule = GENE_RULE

    # test with OrganismDB as annotation source
    # expect_message(assignToGenes(object = obj3_db, anno.annoDB = Homo.sapiens38))
    # expect_message(assignToGenes(object = obj3_db, anno.annoDB = Homo.sapiens38, overlappingLoci = "frequency"))
    # expect_message(expect_error(assignToGenes(object = obj3_db, anno.annoDB = Homo.sapiens38, overlappingLoci = "hierarchy")))
    # expect_message(assignToGenes(object = obj3_db, anno.annoDB = Homo.sapiens38, overlappingLoci = "hierarchy", rule = rule))
    # expect_message(expect_warning(assignToGenes(object = obj3_db, anno.annoDB = Homo.sapiens38, overlappingLoci = "remove")))
    # expect_message(expect_warning(assignToGenes(object = obj3_db, anno.annoDB = Homo.sapiens38, overlappingLoci = "keep")))
    # expect_message(expect_error(assignToGenes(object = obj3_db, anno.annoDB = Homo.sapiens38, overlappingLoci = "test")))

    # test with manual as annotation source
    expect_message(assignToGenes(object = obj3_manual, anno.genes = gns))
    expect_message(assignToGenes(object = obj3_manual, anno.genes = gns, overlaps = "frequency"))
    expect_error(assignToGenes(object = obj3_manual, anno.genes = gns, overlaps = "hierarchy"))
    expect_message(assignToGenes(object = obj3_manual, anno.genes = gns, overlaps = "hierarchy", overlaps.rule = rule))
    expect_warning(assignToGenes(object = obj3_manual, anno.genes = gns, overlaps = "remove"))
    expect_warning(assignToGenes(object = obj3_manual, anno.genes = gns, overlaps = "keep"))
    expect_error(assignToGenes(object = obj3_manual, anno.genes = gns, overlaps = "test"))

    expect_error(assignToGenes(object = obj3_manual, anno.genes = gns, geneMatchID = "GeneID"))
    expect_error(assignToGenes(object = obj3_manual, anno.genes = gns, geneMatchName = "GeneName"))
    expect_error(assignToGenes(object = obj3_manual, anno.genes = gns, geneMatchType = "GeneType"))

    # test quiet
    expect_silent(assignToGenes(object = obj3_manual, anno.genes = gns, quiet = TRUE))

    # test with missing gene annotation meta columns
    gns1 = gns
    gns1$gene_type = NULL
    expect_message(expect_error(assignToGenes(object = obj3_manual, anno.genes = gns1)))
    expect_message(expect_warning(assignToGenes(object = obj3_manual, anno.genes = gns1, overlaps = "remove")))
    expect_message(expect_warning(assignToGenes(object = obj3_manual, anno.genes = gns1, overlaps = "keep")))

    gns1 = gns
    gns1$gene_name = NULL
    expect_message(assignToGenes(object = obj3_manual, anno.genes = gns1))
    expect_message(expect_warning(assignToGenes(object = obj3_manual, anno.genes = gns1, overlaps = "remove")))
    expect_message(expect_warning(assignToGenes(object = obj3_manual, anno.genes = gns1, overlaps = "keep")))

    gns1 = gns
    gns1$gene_id = NULL
    expect_message(assignToGenes(object = obj3_manual, anno.genes = gns1))
    expect_message(expect_warning(assignToGenes(object = obj3_manual, anno.genes = gns1, overlaps = "remove")))
    expect_message(expect_warning(assignToGenes(object = obj3_manual, anno.genes = gns1, overlaps = "keep")))

})

