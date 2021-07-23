library(rtracklayer)
library(GenomicFeatures)

################################################################################
# prepare crosslink sites
################################################################################
# import crosslinks
# csFile = "/Users/mirko/Projects/clip_methods/04_pureclip/PureCLIP.crosslink_sites.bed"
cs = import(con = csFile, format = "BED")
# pureclip score filter
quants = quantile(cs$score, probs = seq(0,1, by = 0.05))
csFilter = cs[cs$score >= quants[2]]


################################################################################
# subset .bw files for those that overlap our example genes
################################################################################
olGen = subsetByOverlaps(gns, csFilter)

# file path
# f = "/Users/mirko/Projects/clip_methods/01_cov/lujh32-combined___bw/DR/raw/"
f = list.files(f, pattern = ".bw$", full.names = TRUE)

# plus strand files
b = rtracklayer::import(f[2], as = "GRanges")
b = subsetByOverlaps(b, olGen)
rtracklayer::export(b, con = "../inst/extdata/rep1_clip_plus.bw", format = "bw")

b = rtracklayer::import(f[4], as = "GRanges")
b = subsetByOverlaps(b, olGen)
rtracklayer::export(b, con = "../inst/extdata/rep2_clip_plus.bw", format = "bw")

b = rtracklayer::import(f[6], as = "GRanges")
b = subsetByOverlaps(b, olGen)
rtracklayer::export(b, con = "../inst/extdata/rep3_clip_plus.bw", format = "bw")

b = rtracklayer::import(f[8], as = "GRanges")
b = subsetByOverlaps(b, olGen)
rtracklayer::export(b, con = "../inst/extdata/rep4_clip_plus.bw", format = "bw")

# minus strand files
b = rtracklayer::import(f[1], as = "GRanges")
b = subsetByOverlaps(b, olGen)
rtracklayer::export(b, con = "../inst/extdata/rep1_clip_minus.bw", format = "bw")

b = rtracklayer::import(f[3], as = "GRanges")
b = subsetByOverlaps(b, olGen)
rtracklayer::export(b, con = "../inst/extdata/rep2_clip_minus.bw", format = "bw")

b = rtracklayer::import(f[5], as = "GRanges")
b = subsetByOverlaps(b, olGen)
rtracklayer::export(b, con = "../inst/extdata/rep3_clip_minus.bw", format = "bw")

b = rtracklayer::import(f[7], as = "GRanges")
b = subsetByOverlaps(b, olGen)
rtracklayer::export(b, con = "../inst/extdata/rep4_clip_minus.bw", format = "bw")



################################################################################
# prepare BindingSiteFinder example dataset
################################################################################
library(BindingSiteFinder)
# clipFiles = "/Users/mirko/Projects/clip_methods/01_cov/lujh32-combined___bw/DR/raw/"
clipFilesP <- list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE)
clipFilesM <- list.files(clipFiles, pattern = "minus.bw$", full.names = TRUE)
meta = data.frame(id = c(1,2,3,4),
                  condition = c("WT", "WT", "KD", "KD"),
                  clPlus = clipFilesP, clMinus = clipFilesM)

# bdsFull = BSFDataSet(ranges = csFilter, meta = meta)
bdsFull = BSFDataSetFromBigWig(ranges = csFilter, meta = meta)
bds = .subsetByChr(bdsFull, chr = "chr22")

save(bds, file = "./bds.rda")


################################################################################
# perpare gene annotaiton data
################################################################################
# annoFile = "/Users/mirko/Projects/PackageDevelopment/BindingSiteDefinition/data/gencode.v37.annotation.gff3"
anno <- makeTxDbFromGFF(file = annoFile)
gns = genes(anno)

# Make annotation database from gff3 file
# annoFile = "/Users/mirko/Projects/PackageDevelopment/BindingSiteDefinition/data/gencode.v37.annotation.chr22.header.gff3"
annoDb = makeTxDbFromGFF(file = annoFile, format = "gff3")
annoInfo = import(annoFile, format = "gff3")
# Get genes as GRanges
gns = genes(annoDb)
idx = match(gns$gene_id, annoInfo$gene_id)
elementMetadata(gns) = cbind(elementMetadata(gns), elementMetadata(annoInfo)[idx,])
# Clean gene object for storage
mcols(gns)$ccdsid = NULL
mcols(gns)$protein_id = NULL
mcols(gns)$havana_transcript = NULL
mcols(gns)$ont = NULL
mcols(gns)$hgnc_id = NULL
mcols(gns)$exon_id = NULL
mcols(gns)$exon_number = NULL
mcols(gns)$transcript_support_level = NULL
mcols(gns)$transcript_name = NULL
mcols(gns)$transcript_id = NULL
mcols(gns)$transcript_type = NULL
mcols(gns)$score = NULL
mcols(gns)$phase = NULL
mcols(gns)$gene_id.1 = NULL
mcols(gns)$Parent = NULL
mcols(gns)$level = NULL
mcols(gns)$tag = NULL
mcols(gns)$havana_gene = NULL
mcols(gns)$source = NULL
mcols(gns)$ID = NULL
mcols(gns)$type = NULL
# Save genes object for load
save(gns, file="gns.rds")

# Count the overlaps of each binding site fore each part of the transcript.
cdseq = cds(annoDb)
intrns = unlist(intronsByTranscript(annoDb))
utrs3 = unlist(threeUTRsByTranscript(annoDb))
utrs5 = unlist(fiveUTRsByTranscript(annoDb))
regions = list(CDS = cdseq, Intron = intrns, UTR3 = utrs3, UTR5 = utrs5)
save(regions, file = "./regions.rds")
