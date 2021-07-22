#' Coverage function for BSFDataSet objects
#'
#' The crosslink coverage is computed for all ranges in the the given
#' \code{BSFDataSet} object (see \code{\link{BSFDataSet}} for details).
#' Depending on the \code{returnOptions} the resulting coverage information is
#' summarized, suitable for diverse computation and plotting tasks. The
#' coverage can only be compute for objects with identical ranges.
#'
#' If \code{returnOptions} is set to merge_ranges_keep_positions:
#' Returns a matrix with ncol being the nucleotides of the ranges (equal to the
#' width of the input ranges) and nrow being the number of replicates in the
#' meta information.
#'
#' If \code{returnOptions} is set to merge_replicates_per_condition:
#' Returns a list of matrices. Each list corresponds to one condition set in
#' the meta information. The matrix in each entry has ncols equal to the
#' ranges width and nrow equal to the number of ranges. Counts per ranges
#' and position are summed.
#'
#' If \code{returnOptions} is set to merge_all_replicates:
#' Returns a matrix with ncols equal to the range width and nrow equal to the
#' number of ranges. Counts per range and position are summed.
#'
#' If \code{returnOptions} is set to merge_positions_keep_replicates:
#' Returns a GRanges object where the counts are summed for each replicate
#' and added to the original granges object.
#'
#' @param object a BSFDataSet object
#' @param returnOptions one of merge_ranges_keep_positions,
#' merge_replicates_per_condition, merge_all_replicates,
#' merge_positions_keep_replicates
#' @param silent TRUE/ FALSE, suppress warning messages
#'
#' @return an object of class specified in \code{returnOptions}
#' @import GenomicRanges
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' rng = coverageOverRanges(
#' bds, returnOptions = "merge_ranges_keep_positions", silent = TRUE)
#' rng = coverageOverRanges(
#' bds, returnOptions = "merge_replicates_per_condition", silent = TRUE)
#' rng = coverageOverRanges(
#' bds, returnOptions = "merge_all_replicates", silent = TRUE)
#' rng = coverageOverRanges(
#' bds, returnOptions = "merge_positions_keep_replicates", silent = TRUE)
#'
#' @export
coverageOverRanges <- function(
    object,
    returnOptions = c("merge_ranges_keep_positions",
                      "merge_replicates_per_condition",
                      "merge_all_replicates",
                      "merge_positions_keep_replicates"),
    silent = FALSE
    ) {
    # Check input
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    # split by strand
    rng = getRanges(object)
    # check for length of ranges
    if (length(unique(width(getRanges(object)))) > 1) {
        stop("width of ranges are not identical")
    }
    # check for unique range names
    if (!isTRUE(silent)) {
        if(any(duplicated(names(rng)))) {
            warning("duplicate names in ranges, adding preceding ID")
            names(rng) = paste0(c(seq_along(rng)), "_", names(rng))
        }
    }
    # split by strand
    rngPlus = rng[strand(rng) == "+"]
    rngMinus = rng[strand(rng) == "-"]
    # prepare signal
    sgn = getSignal(object)
    # extract present conditions
    condition = levels(getMeta(object)$condition)
    # manage return option
    returnOptions = match.arg(
        returnOptions,
        choices = c("merge_ranges_keep_positions",
                    "merge_replicates_per_condition",
                    "merge_all_replicates",
                    "merge_positions_keep_replicates"))

    # Forward to helper functions
    # --------------------------------------------------------------------------
    if (returnOptions == "merge_replicates_per_condition") {
        covRet = .coverageOverRanges.merge_replicates_per_condition(
            sgn = sgn, rng = rng, condition = condition
        )
    }
    if (returnOptions == "merge_all_replicates") {
        covRet = .coverageOverRanges.merge_all_replicates(
            sgn = sgn, rng = rng
        )
    }
    if (returnOptions == "merge_positions_keep_replicates") {
        covRet = .coverageOverRanges.merge_positions_keep_replicates(
            sgn = sgn, rng = rng
        )
    }
    if (returnOptions == "merge_ranges_keep_positions") {
        covRet = .coverageOverRanges.merge_ranges_keep_positions(
            sgn = sgn, rng = rng
        )
    }
    return(covRet)
}



.coverageOverRanges.merge_positions_keep_replicates <- function(sgn, rng) {

    # split by strand
    rngPlus = rng[strand(rng) == "+"]
    rngMinus = rng[strand(rng) == "-"]

    # Plus strand
    # --------------------------------------------------------------------------
    if (length(rngPlus) > 0) {
        # compute initial coverage
        mcols(rngPlus) = as.matrix(
            do.call(cbind, lapply(sgn$signalPlus, function(x) {
                sum(x[rngPlus])
            })))
    }

    # Minus strand
    # --------------------------------------------------------------------------
    if (length(rngMinus) > 0) {
        # compute initial coverage
        mcols(rngMinus) = as.matrix(
            do.call(cbind, lapply(sgn$signalMinus, function(x) {
                sum(x[rngMinus])
            })))
    }

    # Combine strands for return
    # --------------------------------------------------------------------------
    if (length(rngMinus) > 0 & length(rngPlus) > 0) {
        retCov = c(rngPlus, rngMinus)
    }
    if (length(rngMinus) == 0) {
        retCov = rngPlus
    }
    if (length(rngPlus) == 0) {
        retCov = rngMinus
    }
    # sort ranges
    rngCov = c(retCov)
    rngCov = GenomeInfoDb::sortSeqlevels(rngCov)
    rngCov = sort(rngCov)

    return(rngCov)
}

.coverageOverRanges.merge_ranges_keep_positions <- function(sgn, rng) {
    # split by strand
    rngPlus = rng[strand(rng) == "+"]
    rngMinus = rng[strand(rng) == "-"]

    # Plus strand
    # --------------------------------------------------------------------------
    if (length(rngPlus) > 0) {
        # compute initial coverage
        matPlus = lapply(sgn$signalPlus, function(x) {
            y = as.matrix(x[rngPlus])
            rownames(y) = names(rngPlus)
            return(y)
        })
        # summarize
        covPlus = do.call(rbind, lapply(matPlus, colSums))
    }
    if (length(rngPlus) == 0) {
        # set cov on plus strand to zero if no range is present
        covPlus = 0
    }

    # Minus strand
    # --------------------------------------------------------------------------
    if (length(rngMinus) > 0) {
        # compute initial coverage
        matMinus = lapply(sgn$signalMinus, function(x) {
            y = as.matrix(x[rngMinus])
            rownames(y) = names(rngMinus)
            return(y)
        })
        # summarize
        covMinus = do.call(rbind, lapply(matMinus, colSums))
        # reverse orientation
        covMinus = as.matrix(rev(as.data.frame(covMinus)))
    }
    if (length(rngMinus) == 0) {
        # set cov on plus strand to zero if no range is present
        covMinus = 0
    }

    # Combine strands for return
    # --------------------------------------------------------------------------
    if (length(rngMinus) > 0 & length(rngPlus) > 0) {
        retCov = covPlus + covMinus
    }
    if (length(rngMinus) == 0) {
        retCov = covPlus
    }
    if (length(rngPlus) == 0) {
        retCov = covMinus
    }
    return(retCov)
}

.coverageOverRanges.merge_all_replicates <- function(sgn, rng) {
    # split by strand
    rngPlus = rng[strand(rng) == "+"]
    rngMinus = rng[strand(rng) == "-"]

    # Plus strand
    # --------------------------------------------------------------------------
    if (length(rngPlus) > 0) {
        # compute initial coverage
        matPlus = lapply(sgn$signalPlus, function(x) {
            y = as.matrix(x[rngPlus])
            rownames(y) = names(rngPlus)
            return(y)
        })
        # summarize
        covPlus = Reduce('+', matPlus)
    }
    if (length(rngPlus) == 0) {
        # set cov on plus strand to zero if no range is present
        covPlus = 0
    }

    # Minus strand
    # --------------------------------------------------------------------------
    if (length(rngMinus) > 0) {
        # compute initial coverage
        matMinus = lapply(sgn$signalMinus, function(x) {
            y = as.matrix(x[rngMinus])
            rownames(y) = names(rngMinus)
            return(y)
        })
        # summarize
        covMinus = Reduce('+', matMinus)
        # reverse orientation
        covMinus = as.matrix(rev(as.data.frame(covMinus)))
    }
    if (length(rngMinus) == 0) {
        # set cov on plus strand to zero if no range is present
        covMinus = 0
    }

    # Combine strands for return
    # --------------------------------------------------------------------------
    if (length(rngMinus) > 0 & length(rngPlus) > 0) {
        retCov = rbind(covPlus, covMinus)
    }
    if (length(rngMinus) == 0) {
        retCov = covPlus
    }
    if (length(rngPlus) == 0) {
        retCov = covMinus
    }
    return(retCov)
}

.coverageOverRanges.merge_replicates_per_condition <- function(sgn, rng, condition) {
    # split by strand
    rngPlus = rng[strand(rng) == "+"]
    rngMinus = rng[strand(rng) == "-"]

    # Plus strand
    # --------------------------------------------------------------------------
    if (length(rngPlus) > 0) {
        # compute initial coverage
        matPlus = lapply(sgn$signalPlus, function(x) {
            y = as.matrix(x[rngPlus])
            rownames(y) = names(rngPlus)
            return(y)
        })
        # summarize per condition
        matPlusCond = lapply(condition, function(x){
            matPlus[c(grep(x, names(matPlus)))]})
        names(matPlusCond) = condition
        covPlus = lapply(matPlusCond, function(x){Reduce('+', x)})
    }
    if (length(rngPlus) == 0) {
        # set cov on plus strand to zero if no range is present
        covPlus = 0
    }

    # Minus strand
    # --------------------------------------------------------------------------
    if (length(rngMinus) > 0) {
        # compute initial coverage
        matMinus = lapply(sgn$signalMinus, function(x) {
            y = as.matrix(x[rngMinus])
            rownames(y) = names(rngMinus)
            return(y)
        })
        # summarize per condition
        matMinusCond = lapply(condition, function(x){
            matMinus[c(grep(x, names(matMinus)))]})
        names(matMinusCond) = condition
        covMinus = lapply(matMinusCond, function(x){Reduce('+', x)})
        # reverse orientation
        covMinus = lapply(covMinus, function(x){
            as.matrix(rev(as.data.frame(x)))})
    }
    if (length(rngMinus) == 0) {
        # set cov on plus strand to zero if no range is present
        covMinus = 0
    }

    # Combine strands for return
    # --------------------------------------------------------------------------
    if (length(rngMinus) > 0 & length(rngPlus) > 0) {
        retCov = lapply(seq_along(covPlus), function(x){
            rbind(covPlus[[x]],covMinus[[x]])
        })
        names(retCov) = names(covPlus)
    }
    if (length(rngMinus) == 0) {
        retCov = covPlus
    }
    if (length(rngPlus) == 0) {
        retCov = covMinus
    }
    return(retCov)
}
