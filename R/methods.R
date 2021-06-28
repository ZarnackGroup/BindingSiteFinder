#' Show method to for the BSFDataSet
#'
#' Prints the information for each of the slots in the \code{BSFDataSet} object.
#' Ranges of the \code{\link{getRanges}} slot are shown, as well as the number of
#' crosslinks per strand \code{\link{getSignal}} and the levels of the experimental
#' conditions (\code{\link{getMeta}}).
#'
#' @docType methods
#' @name show
#' @rdname show
#' @aliases show show,BSFDataSet-method
#'
#' @param object a \code{BSFDataSet} object
#'
#' @return shows the current object state
#'
#' @seealso \code{\link{BSFDataSet}}
#'
#' @export
setMethod("show",
          "BSFDataSet",
          function(object) {
              cat("Object of class BSFDataSet \n")
              cat("Contained ranges: ",
                  format(
                      length(object@ranges),
                      big.mark = ".",
                      decimal.mark = ","
                  ),
                  "\n")
              cat("----> Number of chromosomes: ",
                  length(levels(seqnames(object@ranges))), "\n")
              if(length(unique(width(object@ranges))) > 6) {
                  cat("----> Ranges width: ", unique(width(object@ranges))[seq(1,6)], "...\n")
              } else {
                  cat("----> Ranges width: ", unique(width(object@ranges)), "\n")
              }
              cat("Contained conditions: ",
                  levels(object@meta$condition),
                  "\n")
          })

#' Accessor method for the ranges of the BSFDataSet object
#'
#' The ranges slot holds the genomic ranges information of the sites currently
#' in the object. They are encoded as a GRanges object with each binding site
#' having a single ranges entry.
#'
#' @docType methods
#' @name getRanges
#' @rdname getRanges
#' @aliases getRanges getRanges,BSFDataSet-method
#'
#' @param object a \code{BSFDataSet} object
#'
#' @return returns the genomic ranges (GRanges) of the associated ranges
#'
#' @seealso \code{\link{BSFDataSet}}
#'
#' @examples
#'
#' # load data
#' csFile <- system.file("extdata", "PureCLIP_crosslink_sites_example.bed",
#'  package="BindingSiteFinder")
#' cs = rtracklayer::import(con = csFile, format = "BED")
#' clipFiles <- system.file("extdata", package="BindingSiteFinder")
#'
#' # two experimental conditions
#' meta = data.frame(condition = factor(c("WT", "WT", "KD", "KD"), levels = c("KD", "WT")),
#' clPlus = list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE),
#' clMinus = list.files(clipFiles, pattern = "minus.bw$", full.names = TRUE))
#' bds = BSFDataSet(ranges = cs, meta = meta)
#'
#' getRanges(bds)
#'
#' @export
setMethod(
    f = "getRanges",
    signature(object = "BSFDataSet"),
    definition = function(object) {
        validObject(object)
        return(object@ranges)
    }
)

#' Setter method for the ranges of the BSFDataSet object
#' The GRanges object that holds the genomic ranges information can be replaced.
#'
#' @docType methods
#' @name setRanges
#' @rdname setRanges
#' @aliases setRanges setRanges,BSFDataSet-method
#'
#' @param object a \code{BSFDataSet} object
#' @param newRanges an object of type \code{GRanges}
#' @param ... additional arguments
#'
#' @return object of type \code{\link{BSFDataSet}} with updated ranges
#'
#' @seealso \code{\link{BSFDataSet}}
#'
#' @examples
#'
#' # load data
#' csFile <- system.file("extdata", "PureCLIP_crosslink_sites_example.bed",
#'  package="BindingSiteFinder")
#' cs = rtracklayer::import(con = csFile, format = "BED")
#' clipFiles <- system.file("extdata", package="BindingSiteFinder")
#'
#' # two experimental conditions
#' meta = data.frame(condition = factor(c("WT", "WT", "KD", "KD"), levels = c("KD", "WT")),
#' clPlus = list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE),
#' clMinus = list.files(clipFiles, pattern = "minus.bw$", full.names = TRUE))
#' bds = BSFDataSet(ranges = cs, meta = meta)
#'
#' rng = getRanges(bds)
#' rng = rng + 10
#' bdsNew = setRanges(bds, rng)
#'
#' @export
setMethod(
    f = "setRanges",
    signature(object = "BSFDataSet"),
    definition = function(object, newRanges) {
        object@ranges <- newRanges
        validObject(object)
        return(object)
    }

)

#' Accessor method for the meta data of the BSFDataSet object
#'
#' Meta data is stored as a \code{data.frame} and must contain the columns
#' "condition", "clPlus" and "clMinus".
#'
#' @docType methods
#' @name getMeta
#' @rdname getMeta
#' @aliases getMeta getMeta,BSFDataSet-method
#'
#' @param object a BSFDataSet object
#'
#' @return returns the meta data \code{data.frame} with the columns "condition",
#' "clPlus" and "clMinus".
#'
#' @seealso \code{\link{BSFDataSet}}
#'
#' @examples
#'
#' # load data
#' csFile <- system.file("extdata", "PureCLIP_crosslink_sites_example.bed",
#'  package="BindingSiteFinder")
#' cs = rtracklayer::import(con = csFile, format = "BED")
#' clipFiles <- system.file("extdata", package="BindingSiteFinder")
#'
#' # two experimental conditions
#' meta = data.frame(condition = factor(c("WT", "WT", "KD", "KD"), levels = c("KD", "WT")),
#' clPlus = list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE),
#' clMinus = list.files(clipFiles, pattern = "minus.bw$", full.names = TRUE))
#' bds = BSFDataSet(ranges = cs, meta = meta)
#'
#' getMeta(bds)
#'
#' @export
setMethod(
    f = "getMeta",
    signature(object = "BSFDataSet"),
    definition = function(object) {
        validObject(object)
        return(object@meta)
    }
)


#' Accessor method for the signal data of the BSFDataSet object
#'
#' Signal data is loaded from the path specified in \code{\link{getMeta}} columns
#' "clPlus" and "clMinus" and stored as a list of RLE lists.
#'
#' @docType methods
#' @name getSignal
#' @rdname getSignal
#' @aliases getSignal getSignal,BSFDataSet-method
#'
#' @param object a BSFDataSet object
#'
#' @return returns the signal data, as list of RLE list for each strand, named
#' after the meta data columns "clPlus" and "clMinus"
#'
#' @seealso \code{\link{getMeta}} \code{\link{BSFDataSet}}
#'
#' @examples
#'
#' # load data
#' csFile <- system.file("extdata", "PureCLIP_crosslink_sites_example.bed",
#'  package="BindingSiteFinder")
#' cs = rtracklayer::import(con = csFile, format = "BED")
#' clipFiles <- system.file("extdata", package="BindingSiteFinder")
#'
#' # two experimental conditions
#' meta = data.frame(condition = factor(c("WT", "WT", "KD", "KD"), levels = c("KD", "WT")),
#' clPlus = list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE),
#' clMinus = list.files(clipFiles, pattern = "minus.bw$", full.names = TRUE))
#' bds = BSFDataSet(ranges = cs, meta = meta)
#'
#' getSignal(bds)
#'
#' @export
setMethod(
    f = "getSignal",
    signature(object = "BSFDataSet"),
    definition = function(object) {
        validObject(object)
        return(object@signal)
    }
)

#' Setter method for the signal data of the BSFDataSet object
#'
#' Signal data is loaded from the path specified in \code{\link{getMeta}} columns
#' "clPlus" and "clMinus" and stored as a list of RLE lists.
#'
#' @docType methods
#' @name setSignal
#' @rdname setSignal
#' @aliases setSignal setSignal,BSFDataSet-method
#'
#' @param object a BSFDataSet object
#' @param newSignal list of RLE lists
#' @param ... additional arguments
#'
#' @return an object of type \code{\link{BSFDataSet}} with updated signal
#'
#' @seealso \code{\link{BSFDataSet}}
#'
#' @examples
#'
#' # load data
#' csFile <- system.file("extdata", "PureCLIP_crosslink_sites_example.bed",
#'  package="BindingSiteFinder")
#' cs = rtracklayer::import(con = csFile, format = "BED")
#' clipFiles <- system.file("extdata", package="BindingSiteFinder")
#'
#' # two experimental conditions
#' meta = data.frame(condition = factor(c("WT", "WT", "KD", "KD"), levels = c("KD", "WT")),
#' clPlus = list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE),
#' clMinus = list.files(clipFiles, pattern = "minus.bw$", full.names = TRUE))
#' bds = BSFDataSet(ranges = cs, meta = meta)
#'
#' sgn = getSignal(bds)
#' sgn = lapply(sgn, function(selStrand){
#'    lapply(selStrand, function(chrList){
#'        chrList[names(chrList) == "chr1"]
#'    })
#' })
#' bdsNew = setSignal(bds, sgn)
#'
#' @export
setMethod(
    f = "setSignal",
    signature(object = "BSFDataSet"),
    definition = function(object, newSignal) {
        object@signal <- newSignal
        validObject(object)
        return(object)
    }
)

#' Accessor method for the summary slot of the BSFDataSet object
#'
#' The summary slot is used to track information of the filtering steps applied in
#' the \code{\link{makeBindingSites}} function
#'
#' @docType methods
#' @name getSummary
#' @rdname getSummary
#' @aliases getSummary getSummary,BSFDataSet-method
#'
#' @param object a \code{BSFDataSet} object
#' @param ... additional arguments
#'
#' @return returns the summary information storted in the summary slot after
#' \code{\link{makeBindingSites}} was run
#'
#' @seealso \code{\link{BSFDataSet}} \code{\link{makeBindingSites}}
#'
#' @examples
#'
#' # load data
#' csFile <- system.file("extdata", "PureCLIP_crosslink_sites_example.bed",
#'  package="BindingSiteFinder")
#' cs = rtracklayer::import(con = csFile, format = "BED")
#' clipFiles <- system.file("extdata", package="BindingSiteFinder")
#'
#' # two experimental conditions
#' meta = data.frame(condition = factor(c("WT", "WT", "KD", "KD"), levels = c("KD", "WT")),
#' clPlus = list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE),
#' clMinus = list.files(clipFiles, pattern = "minus.bw$", full.names = TRUE))
#' bds = BSFDataSet(ranges = cs, meta = meta)
#'
#' bds <- makeBindingSites(object = bds, bsSize = 9, minWidth = 2,
#' minCrosslinks = 2, minClSites = 1)
#'
#' getSummary(bds)
#'
#' @export
setMethod(
    f = "getSummary",
    signature(object = "BSFDataSet"),
    definition = function(object) {
        validObject(object)
        if (nrow(object@summary) == 0) {
            warning(
                "Function makeBindingSites() was not run.
                 Run makeBindingSites() to get a summary on the merging steps."
            )
        }
        if (nrow(object@summary) > 0) {
            return(object@summary)
        }
    }
)

#' Setter method for the summary slot of the BSFDataSet object
#'
#' The summary slot is used to track information of the filtering steps applied in
#' the \code{\link{makeBindingSites}} function
#'
#' @docType methods
#' @name setSummary
#' @rdname setSummary
#' @aliases setSummary setSummary,BSFDataSet-method
#'
#' @param object a BSFDataSet object
#' @param summary a data.frame with the summary information to be stored in \code{BSFDataSet}
#' @param ... additional arguments
#'
#' @return an object of type \code{\link{BSFDataSet}} with updated summary info
#'
#' @seealso \code{\link{BSFDataSet}}
#'
#' @examples
#'
#' # load data
#' csFile <- system.file("extdata", "PureCLIP_crosslink_sites_example.bed",
#'  package="BindingSiteFinder")
#' cs = rtracklayer::import(con = csFile, format = "BED")
#' clipFiles <- system.file("extdata", package="BindingSiteFinder")
#'
#' # two experimental conditions
#' meta = data.frame(condition = factor(c("WT", "WT", "KD", "KD"), levels = c("KD", "WT")),
#' clPlus = list.files(clipFiles, pattern = "plus.bw$", full.names = TRUE),
#' clMinus = list.files(clipFiles, pattern = "minus.bw$", full.names = TRUE))
#' bds = BSFDataSet(ranges = cs, meta = meta)
#'
#' df = data.frame(processingStep = c(1,2),
#' parameter = c(3,4))
#' bds = setSummary(bds, df)
#'
#' @export
setMethod(
    f = "setSummary",
    signature(object = "BSFDataSet"),
    definition = function(object, summary) {
        object@summary <- summary
        validObject(object)
        return(object)
    }

)
