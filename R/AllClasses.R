#' @rdname BSFDataSet
#' @export
BSFDataSet <- setClass(
    "BSFDataSet",
    representation = list(
        ranges = "GRanges",
        meta = "data.frame",
        signal = "list",
        summary = "data.frame",
        params = "list",
        plotData = "list",
        results = "data.frame"
    )
)
setValidity("BSFDataSet", function(object) {
    msg <- NULL

    # check input ranges
    if (any(strand(object@ranges) == "*")) {
        msg = c(msg, "Strand can not be of type'*'.")
    }
    # check input meta data
    if (!c("id") %in% colnames(object@meta) |
        !c("condition") %in% colnames(object@meta)) {
        msg = c(msg, "Metadata must contain a column named 'id' and
                a column named 'condtion'.")
    }
    if (length(object@meta$id) != length(unique(object@meta$id))) {
        msg = c(msg, "IDs in meta data are not unique. ")
    }
    # check signal
    if (!c("signalPlus") %in% names(object@signal) &
        !c("signalMinus") %in% names(object@signal)) {
        msg = c(msg, "Signal slot must at least contain one element named
                'signalPlus', or 'signalMinus'. ")
    }
    if (is.null(object@signal$signalPlus) &
        is.null(object@signal$signalMinus)) {
        msg = c(msg, "Signal slot can not be empty. ")
    }
    if (isTRUE(!all(names(object@signal) %in% c("signalPlus") |
                    names(object@signal) %in% c("signalMinus")))) {
        msg = c(msg, "Incorrect name in signal list. ")
    }
    # check signal list structure
    if(length((match(c("signalPlus", "signalMinus"), names(object@signal), nomatch = 0) > 0)) != 2){
        msg = c(msg, "Signal must contain exactly two sublists. ")
    }
    if(!all(match(c("signalPlus", "signalMinus"), names(object@signal), nomatch = 0) > 0)) {
        # check if signal list contains plus/ minus sublists
        msg = c(msg, "Signal must contain lists named 'signalPlus' and 'signalMinus'. ")
    }
    if (!is.null(object@signal$signalPlus)) {
        if (!all(unlist(lapply(object@signal$signalPlus,
                               function(x){is(x, "SimpleRleList")})))) {
            msg = c(msg, "Signal plus elements are not of type RleList")
        }
    }
    if (!is.null(object@signal$signalMinus)) {
        if (!all(unlist(lapply(object@signal$signalMinus,
                               function(x){is(x, "SimpleRleList")})))) {
            msg = c(msg, "Signal minus elements are not of type RleList")
        }
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
#' store input ranges. Alongside with the iCLIP signal in \code{list} structure
#' and additional meta data as \code{data.frame}.
#'
#' The ranges are enforced to have to have a "+" or "-" strand
#' annotation,"*" is not allowed. They are expected to be of the same width
#' and a warning is thrown otherwise.
#'
#' The meta information is stored as \code{data.frame} with at least two
#' required columns, 'id' and 'condition'. They are used to build the unique
#' identifier for each replicate split by '_' (eg. id = 1 and condition = WT
#' will result in 1_WT).
#'
#' The meta data needs to have the additional columns 'clPlus' and 'clMinus'
#' to be present if \code{BSFDataSetFromBigWig} is called.
#' It is used to provide the location to the iCLIP coverage
#' files to the import function. On object initialization these files are
#' loaded and internally represented in the signal slot of the object
#' (see \code{\link{BSFDataSet}}).
#'
#' The iCLIP signal is stored in a special list structure.
#' At the lowest level crosslink counts per nucleotide are stored as
#' \code{Rle} per chromosome summarized as a \code{SimpleRleList}. Such a list
#' exits for each replicate and must be named by the replicate identifier
#' (eg. 1_WT). Therefore this list contains always exactly the same number of
#' entries as the number of replicates in the dataset. Since we handle strands
#' initially seperated from each other this list must be given twice, once
#' for each strand. The strand specific entries must be named 'signalPlus' and
#' 'signalMinus'.
#'
#' The option \code{dropSeqlevels} forces the seqnames of the ranges and
#' the signal to be the same. If for a specific chromosome in the ranges
#' no respective entry in the signal list can be found, then entries with
#' that chromosome are dropped This behavior is needed to keep the
#' \link{BSFDataSet} object in sync, which is required for downstream functions
#' such as \code{coverageOverRanges}
#'
#' @param ranges a \code{GenomicRanges} with the desired ranges to process. The
#' strand slot must be either + or -.
#' @param meta a \code{data.frame} with at least two columns. The first column
#' should be a unique numeric id. The second column holds sample type
#' information, such as the condition.
#' @param signal a \code{list} with the two entries 'signalPlus' and
#' 'signalMinus', following a special representation of \code{SimpleRleList}
#' for counts per replicates (see details for more information).
#' @param silent suppress messages but not warnings (TRUE/ FALSE)
#' @param dropSeqlevels enforce seqnames to be the same in ranges and signal,
#' by dropping unused seqlevels which is required for most downstream functions
#' such as \code{coverageOverRanges}
#'
#' @return A BSFDataSet object.
#'
#' @aliases BSFDataSet, BSFDataSet-class, BSFDataSetFromBigWig
#'
#' @docType class
#'
#' @import GenomicRanges
#' @importFrom rtracklayer import
#'
#' @examples
#'
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' rng = getRanges(bds)
#' sgn = getSignal(bds)
#' mta = getMeta(bds)
#' bdsNew = BSFDataSet(ranges = rng, signal = sgn, meta = mta)
#'
#' @rdname BSFDataSet
#' @export
BSFDataSet <- function(ranges, meta, signal, dropSeqlevels = TRUE, silent = FALSE) {
    msg = NULL
    # check input ranges
    if (!all(c(any(strand(ranges) == "-"), any(strand(ranges) == "+")))) {
        warning("Input ranges are only on one strand. \n")
    }
    if (length(unique(width(ranges))) > 1) {
        msg = c(msg, "Ranges are of different width. \n")
        # warning("Ranges are of differnt width. ")
    }
    if (any(width(ranges) > 1)) {
        msg = c(msg, "Input ranges are larger than 1 nt. \n")
        # message("Input ranges are larger than 1 nt. ")
    }
    # check input metadata
    if(!all(c(any(colnames(meta) == "id"),
              any(colnames(meta) == "condition")))){
        msg = c(msg, "Meta data columns must contain 'id' and 'condition'. \n")
        # stop("Meta data columns must contain 'id' and 'condition'. ")
    }
    if (!is.factor(meta$condition)) {
        msg = c(msg, "Condition column is not a factor, converting to factor. \n")
        # message("Condition column is not factor, converting to factor. ")
        meta$condition = factor(meta$condition)
    }
    if (any(duplicated(meta$clPlus)) |
        any(duplicated(meta$clMinus))) {
        msg = c(msg, "Given paths are duplicated. Please check your input. \n")
        # message("Given path are duplicated. Please check your input.")
    }
    if (! identical(.sortRanges(ranges), ranges)){
        msg = c(msg, "Input ranges are not sorted, sorting for you. \n")
        # message('Input ranges are not sorted, sorting for you.')
        ranges = .sortRanges(ranges)
    }

    # check for signal and ranges integrity
    fixed = .checkForDropSeqlevels(ranges = ranges, signal = signal,
                                dropSeqlevels = dropSeqlevels)
    ranges = fixed$ranges
    signal = fixed$signal
    if (!is.null(fixed$msg)) {
        msg = c(msg, fixed$msg)
    }

    if (!isTRUE(silent)) {
        if(!is.null(msg)) {
            message(msg)
        }
    }

    # set placeholder for summary slot
    summary = data.frame()
    params = list()
    plotData = list()
    results = data.frame()

    # construct final object
    obj = new(
        "BSFDataSet",
        ranges = ranges,
        meta = meta,
        signal = signal,
        summary = summary,
        params = params,
        plotData = plotData,
        results = results
    )
    return(obj)
}

#' @rdname BSFDataSet
#' @export
BSFDataSetFromBigWig <- function(ranges, meta, silent = FALSE,
                                 dropSeqlevels = TRUE) {
    # check the meta data dataframe for additional info where to find
    # the big wig files
    if(!all(c(
        any(colnames(meta) == "id"),
        any(colnames(meta) == "condition"),
        any(colnames(meta) == "clPlus"),
        any(colnames(meta) == "clMinus")))){
        stop("Meta data columns must contain 'id', 'condition' and
             'clPlus', 'clMinus'. ")
    }
    # check if given path point to actual files
    if (all(!file.exists(meta$clPlus) &
            !file.exists(meta$clMinus))) {
        stop("Given path do not point to existing files. ")
    }
    if (any(duplicated(meta$clPlus)) |
        any(duplicated(meta$clMinus))) {
        warning("Given path are duplicated. Please check your input.")
    }
    # build colSignal by importing files as RLE
    signalPlus = unlist(lapply(meta$clPlus, function(x) {
        rtracklayer::import(con = x, as = "RleList")
    }))
    signalMinus = unlist(lapply(meta$clMinus, function(x) {
        abs(rtracklayer::import(con = x, as = "RleList"))
    }))
    names(signalPlus) = paste0(meta$id, "_", meta$condition)
    names(signalMinus) = paste0(meta$id, "_", meta$condition)
    signal = list(signalPlus = signalPlus, signalMinus = signalMinus)

    # construct BindingSiteFinder data set
    object = BSFDataSet(ranges = ranges, meta = meta, signal = signal,
                        dropSeqlevels = dropSeqlevels, silent = silent)
    return(object)
}

