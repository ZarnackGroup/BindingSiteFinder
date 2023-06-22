#' Show method to for the BSFDataSet
#'
#' Prints the information for each of the slots in the \code{BSFDataSet} object.
#' Ranges of the \code{\link{getRanges}} slot are shown, as well as the number
#' of crosslinks per strand \code{\link{getSignal}} and the levels of the
#' experimental conditions (\code{\link{getMeta}}).
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
#' @examples
#'
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' show(bds)
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
                  cat("----> Ranges width: ",
                      unique(width(object@ranges))[seq(1,6)], "...\n")
              } else {
                  cat("----> Ranges width: ",
                      unique(width(object@ranges)), "\n")
              }
              cat("Contained Signal:",
                  format(sum(sapply(object@signal$signalPlus, sum) +
                                 sapply(object@signal$signalMinus, sum)),
                         big.mark = ","), "\n")
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
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
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
#' @param dropSeqlevels enforce seqnames to be the same in ranges and signal,
#' by dropping unused seqlevels which is required for most downstream functions
#' such as \code{coverageOverRanges}
#' @param quiet logical; whether to print messages
#' @param ... additional arguments
#'
#' @return object of type \code{\link{BSFDataSet}} with updated ranges
#'
#' @seealso \code{\link{BSFDataSet}}
#'
#' @examples
#'
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' rng = getRanges(bds)
#' rng = rng + 10
#' bdsNew = setRanges(bds, rng)
#'
#' @export
setMethod(
    f = "setRanges",
    signature(object = "BSFDataSet"),
    definition = function(object, newRanges, dropSeqlevels = TRUE, quiet = FALSE) {
        msg = NULL
        # check for ranges and signal integrity
        fixed = .checkForDropSeqlevels(ranges = newRanges, signal = object@signal,
                                    dropSeqlevels = dropSeqlevels)

        object@ranges <- fixed$ranges
        object@signal <- fixed$signal
        if (!is.null(fixed$msg)) {
            if(!quiet) message(fixed$msg)
        }
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
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
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


#' Setter method for the meta data of the BSFDataSet object
#'
#' Meta data is stored as a \code{data.frame} and must contain the columns
#' "condition", "clPlus" and "clMinus".
#'
#' @docType methods
#' @name setMeta
#' @rdname setMeta
#' @aliases setMeta setMeta,BSFDataSet-method
#'
#' @param object a BSFDataSet object
#' @param newMeta the replacement meta data table
#' @param ... additional arguments
#'
#' @return an object of type \code{\link{BSFDataSet}} with updated meta data
#'
#' @seealso \code{\link{BSFDataSet}}
#'
#' @examples
#'
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' nMeta = getMeta(bds)
#' setMeta(bds, nMeta)
#'
#' @export
setMethod(
    f = "setMeta",
    signature(object = "BSFDataSet"),
    definition = function(object, newMeta) {
        object@meta = newMeta

        validObject(object)
        return(object)
    }
)



#' Accessor method for the signal data of the BSFDataSet object
#'
#' Signal data is loaded from the path specified in \code{\link{getMeta}}
#' columns "clPlus" and "clMinus" and stored as a list of RLE lists.
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
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
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
#' Signal data is loaded from the path specified in \code{\link{getMeta}}
#' columns "clPlus" and "clMinus" and stored as a list of RLE lists.
#'
#' @docType methods
#' @name setSignal
#' @rdname setSignal
#' @aliases setSignal setSignal,BSFDataSet-method
#'
#' @param object a BSFDataSet object
#' @param newSignal list of RLE lists
#' @param dropSeqlevels enforce seqnames to be the same in ranges and signal,
#' by dropping unused seqlevels which is required for most downstream functions
#' such as \code{coverageOverRanges}
#' @param quiet logical; whether to print messages
#' @param ... additional arguments
#'
#' @return an object of type \code{\link{BSFDataSet}} with updated signal
#'
#' @seealso \code{\link{BSFDataSet}}
#'
#' @examples
#'
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' sgn = getSignal(bds)
#' sgn = lapply(sgn, function(selStrand){
#'    lapply(selStrand, function(chrList){
#'        chrList[names(chrList) == "chr22"]
#'    })
#' })
#' bdsNew = setSignal(bds, sgn)
#'
#' @export
setMethod(
    f = "setSignal",
    signature(object = "BSFDataSet"),
    definition = function(object, newSignal, dropSeqlevels = TRUE, quiet = FALSE) {
        msg = NULL
        # check for ranges and signal integrity
        fixed = .checkForDropSeqlevels(ranges = object@ranges, signal = newSignal,
                                       dropSeqlevels = dropSeqlevels)

        object@ranges <- fixed$ranges
        object@signal <- fixed$signal
        if (!is.null(fixed$msg)) {
            if(!quiet) message(fixed$msg)
        }
        validObject(object)
        return(object)
    }
)

#' Accessor method for the summary slot of the BSFDataSet object
#'
#' The summary slot is used to track information of the filtering steps applied
#' in the \code{\link{makeBindingSites}} function
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
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
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
#' The summary slot is used to track information of the filtering steps applied
#' in the \code{\link{makeBindingSites}} function
#'
#' @docType methods
#' @name setSummary
#' @rdname setSummary
#' @aliases setSummary setSummary,BSFDataSet-method
#'
#' @param object a BSFDataSet object
#' @param summary a data.frame with the summary information to be stored in
#' \code{BSFDataSet}
#' @param ... additional arguments
#'
#' @return an object of type \code{\link{BSFDataSet}} with updated summary info
#'
#' @seealso \code{\link{BSFDataSet}}
#'
#' @examples
#'
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
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


#' Subset a BSFDataSet object
#'
#' You can subset \link{BSFDataSet} by identifier or by position using the
#' \code{`[`} operator. Empty seqlevels are being droppend after the subset.
#'
#' @param x A \link{BSFDataSet} object.
#' @param i Position of the identifier or the name of the identifier itself.
#' @param j Not used.
#' @param ... Additional arguments not used here.
#' @param drop if the signal not covered by the subsetted ranges should be
#' dropped or not
#'
#' @return A \link{BSFDataSet} object.
#'
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomeInfoDb keepSeqlevels
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' bdsNew = bds[1:10]
#'
#' @name subset-BSFDataSet
NULL

#' @rdname subset-BSFDataSet
#' @export
setMethod("[", signature(x = "BSFDataSet", i = "ANY", j = "ANY", drop = "ANY"),
          function(x, i, j, ..., drop=FALSE)
          {
              rng = x@ranges
              sgn = x@signal
              rngSub = rng[i]

              # drop empty seqlevels from ranges
              rngSub = GenomeInfoDb::keepSeqlevels(rngSub,
                                                   value = unique(seqnames(rngSub)))

              # reduce signal to only those seqlevels that are still present in
              # the ranges after subsetting
              sgnSub = lapply(sgn, function(currStrand){
                  lapply(currStrand, function(currSample){
                      currSample = currSample[match(names(currSample), seqnames(rngSub), nomatch = 0) > 0]
                  })
              })

              # drop unused seqlevels
              if (isTRUE(drop)) {
                  # drop signal on all regions not covered by binding sites
                  sgnSub = lapply(sgn, function(currStrand){
                      lapply(currStrand, function(currSample){
                          currSample = currSample[match(names(currSample), seqnames(rngSub), nomatch = 0) > 0]
                          as(lapply(currSample, function(currChr){
                              rngCov = coverage(ranges(rngSub), width = length(currChr))
                              currChr[rngCov == 0] = 0
                              currChr
                          }), "SimpleRleList")
                      })
                  })
              }
              initialize(x, ranges = rngSub, signal = sgnSub)
          }
)
