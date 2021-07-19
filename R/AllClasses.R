#' @rdname BSFDataSet
#' @export
BSFDataSet <- setClass(
    "BSFDataSet",
    representation = list(
        ranges = "GRanges",
        meta = "data.frame",
        signal = "list",
        summary = "data.frame"
    )
)
setValidity("BSFDataSet", function(object) {
    msg <- NULL

    # check input ranges
    if (any(strand(object@ranges) == "*")) {
        msg = c(msg, "Strand can not be of type'*'.")
    }
    # check input meta data
    if (ncol(object@meta) < 3) {
        msg = c(msg, "Metadata must contain at least 3 columns.")
    }
    if (!c("clPlus") %in% colnames(object@meta) ||
        !c("clMinus") %in% colnames(object@meta)) {
        msg = c(msg,
                "Metadata must contain columns named 'clPlus' and 'clMinus' ")
    }
    if (!c("condition") %in% colnames(object@meta)) {
        msg = c(msg, "Metadata must contain a column named 'condtion'.")
    }
    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }
})

#' BSFDataSet object and constructors
#'
#' \code{BSFDataSet} contains the class \code{GenomicRanges}, which is used to
#' store input ranges. It enforces for all ranges to have a "+" or "-" strand
#' annotation,"*" is not allowed. Alongside these ranges meta data is stored
#' as \code{data.frame}.This dataframe needs to contain three columns which
#' must be named "condition", "clPlus" and "clMinus". It is used to provide the
#' location to the iCLIP coverage files to the import function. On object
#' initialization these files are loaded and internally represented
#' as RLE-Lists. See the vignette for different construction examples.
#'
#' @param ranges a \code{GenomicRanges} with the desired ranges to process. The
#' strand slot must be either + or -.
#' @param meta a \code{data.frame} with a minimum of three columns. The first
#' column holds sample type information, such as the condition. The second and
#' third column must be named 'clPlus' and 'clMinus' and must contain the path
#' to the files with the individual crosslink events. This can be in
#' .bw (bigwig) format for example.
#' See the vignette for how such files can be obtained from sequenced reads.
#' @param forceEqualNames to maintain the integrity of chromosome
#' names (TRUE/ FALSE). The option ensures that chromosome names present in
#' the GRanges are also all present in the signal list and vice versa.
#' Chromosomes names present in only the signal list or the ranges are removed.
#' @param silent suppress loading message (TRUE/ FALSE)
#'
#' @return A BSFDataSet object.
#'
#' @docType class
#'
#' @importFrom rtracklayer import
#' @import GenomicRanges
#'
#' @examples
#' if (.Platform$OS.type != "windows") {
#'    csFile <- system.file("extdata", "PureCLIP_crosslink_sites_example.bed",
#'                          package="BindingSiteFinder")
#'    cs = rtracklayer::import(con = csFile, format = "BED")
#'    clipFiles <- system.file("extdata", package="BindingSiteFinder")
#'
#'    # do only if not windows
#'    meta = data.frame(condition = factor(c("WT", "WT", "KD", "KD"),
#'     levels = c("KD", "WT")),
#'                      clPlus = list.files(clipFiles, pattern = "plus.bw$",
#'                       full.names = TRUE),
#'                      clMinus = list.files(clipFiles, pattern = "minus.bw$",
#'                       full.names = TRUE))
#'    bds = BSFDataSet(ranges = cs, meta = meta)
#'
#'    # one experimental condition
#'    meta = data.frame(condition = c("WT", "WT", "WT", "WT"),
#'                      clPlus = list.files(clipFiles, pattern = "plus.bw$",
#'                       full.names = TRUE),
#'                      clMinus = list.files(clipFiles, pattern = "minus.bw$",
#'                       full.names = TRUE))
#'    bds = BSFDataSet(ranges = cs, meta = meta, forceEqualNames = TRUE)
#'
#' }
#'
#' @rdname BSFDataSet
#' @export
BSFDataSet <- function(ranges, meta, forceEqualNames = TRUE, silent = FALSE) {
    # check input ranges
    if (!all(c(any(strand(ranges) == "-"), any(strand(ranges) == "+")))) {
        warning("Input ranges are only on one strand. ")
    }
    if (length(unique(width(ranges))) > 1) {
        warning("Ranges are of differnt width. ")
    }
    if (any(width(ranges) > 1)) {
        message("Input ranges are larger than 1 nt. ")
    }
    # check input metadata
    if(!all(c(any(colnames(meta) == "condition"),
              any(colnames(meta) == "clPlus"),
              any(colnames(meta) == "clMinus")))){
        stop("Meta data columns must contain 'condition', 'clPlus', 'clMinus'")
    }

    if (!is.factor(meta$condition)) {
        message("Condition column is not factor, converting to factor. ")
        meta$condition = factor(meta$condition)
    }

    # Loading message
    if (!isTRUE(silent)){
        message("Importing ranges. ")
    }
    # build colSignal by importing files as RLE
    signalPlus = unlist(lapply(meta$clPlus, function(x) {
        rtracklayer::import(con = x, as = "RleList")
    }))
    signalMinus = unlist(lapply(meta$clMinus, function(x) {
        abs(rtracklayer::import(con = x, as = "RleList"))
    }))
    names(signalPlus) = paste0(seq_len(nrow(meta)), "_", meta$condition)
    names(signalMinus) = paste0(seq_len(nrow(meta)), "_", meta$condition)
    signal = list(signalPlus = signalPlus, signalMinus = signalMinus)

    # check input signal
    if (isTRUE(forceEqualNames)) {
        rngChrs = as.character(unique(seqnames(ranges)))
        # fix input signal
        signal = lapply(signal, function(selStrand) {
            lapply(selStrand, function(chrList) {
                chrList[names(chrList) %in% rngChrs]
            })
        })
        # fix input ranges
        sgnChrs = names(signal[[1]][[1]])
        rngChrs = rngChrs[rngChrs %in% sgnChrs]
        ranges = ranges[as.character(seqnames(ranges)) %in% rngChrs]
    }
    if (!isTRUE(forceEqualNames)) {
        rngChrs = unique(seqnames(ranges))
        sgnChrs = names(signal[[1]][[1]])
        # check signal
        if (!all(rngChrs %in% sgnChrs) & all(sgnChrs %in% rngChrs)) {
            warning("forceEqualNames is FALSE and chromosome names in the
                    ranges and signal objects do not match. ")
        }
    }

    # set placeholder for summary slot
    summary = data.frame()

    # construct final object
    obj = new(
        "BSFDataSet",
        ranges = ranges,
        meta = meta,
        signal = signal,
        summary = summary
    )
    return(obj)
}
