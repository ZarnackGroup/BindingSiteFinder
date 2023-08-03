# Prepare BSFDataSet for testing

# load packages
library(BindingSiteFinder)

# signal reduction function
# -> this function is part of the package, but not exported to the user
# -> we need it here once ot reduce the file size of the test dataset
.reduceSignalToFrame <- function(object, frame, quiet = FALSE) {
    rng = getRanges(object)
    # extend ranges by desired frame on which signal should be kept
    rng = rng + frame
    newObj = setRanges(object, rng, quiet = quiet)
    # drop signal outside of protective frame
    newObj = newObj[seq_along(rng), drop=TRUE]
    # put original sized ranges back in place
    rng = rng-frame
    newObj = setRanges(newObj, rng, quiet = quiet)

    return(newObj)
}

# load all datasets
TestSets = readRDS(file = "/Users/mirko/Projects/Datasets/01_data/TEST_DATASETS_FULL.rds")
# select the example case
this.set = TestSets$SRSF6_Diabetes

# reduce crosslink sites to chr22 for size
cs = getRanges(this.set)
cs = subset(cs, seqnames(cs) == "chr22")
seqlevels(cs) = "chr22"
this.new = setRanges(this.set, cs)

# reduce signal for size
this.new = .reduceSignalToFrame(this.new, frame = 50)

# export object
bds = setName(bds, "Test set")
save(bds, file = "./bds.rda")



# --- old
# import example crosslinks
# csFile = "/Users/mirko/Projects/clip_methods/04_pureclip/PureCLIP.crosslink_sites.bed"
# cs = import(con = csFile, format = "BED")
#
# # reduce sites to chr22 for size
# cs = subset(cs, seqnames(cs) == "chr22")
# seqlevels(cs) = "chr22"
#
# # pureclip score filter -> reduce file size
# quants = quantile(cs$score, probs = seq(0,1, by = 0.05))
# csFilter = cs[cs$score >= quants[2]]
#
# # make meta data
# clipFiles = "/Users/mirko/Projects/clip_methods/01_cov/lujh32-combined___bw/DR/raw/"
# clipFilesP <- list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE)
# clipFilesM <- list.files(clipFiles, pattern = "minus.bw$", full.names = TRUE)
# meta = data.frame(id = c(1,2,3,4),
#                   condition = c("WT", "WT", "KD", "KD"),
#                   clPlus = clipFilesP, clMinus = clipFilesM)
#
# # make BSFDataSet
# bdsFull = BSFDataSetFromBigWig(ranges = csFilter, meta = meta)
#
# # reduce signal for size
# bds = .reduceSignalToFrame(bdsFull, frame = 50)
# bds = setName(bds, "Test set")
#
# # export object
# save(bds, file = "./bds.rda")
