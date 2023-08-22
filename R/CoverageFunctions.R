

#' Coverage function for BSFDataSet objects
#'
#' Function that computes a crosslink coverage with all samples over all ranges
#' given in the \code{\link{BSFDataSet}}. The coverage can be summarized over
#' all combinations of the three dimension (samples, ranges, positions).
#'
#' When summarizing the crosslink coverage over samples (\code{samples.merge=TRUE})
#' one can decide whether to summarize all samples or whether to keep conditions
#' separate (\code{samples.group}). This either reduces the samples dimension
#' to a single matrix, or a list. For a binding site set with 100 binding sites
#' of width=7 and 4 replicates with 2 conditions, the following options are
#' possible. With merging enabled and \code{samples.group='all'} the coverage of
#' all samples is combined. With \code{samples.group='condition'} only samples
#' of the same condition are grouped.
#'
#' When summarizing the crosslink coverage over ranges, all ranges are combined
#' which reduces the ranges dimension to a single vector. This turns eg. a
#' binding site set of 100 binding sites with width=7 into a vector of length
#' 100 with exactly one column. Depending on how the samples were summarized, the
#' result can be a single such vector, or a list.
#'
#' When summarizing the crosslink coverage over positions, all positions are
#' combined which reduces the positions dimension to a single vector. This turns
#' eq. a binding site set of 100 binding sites with width=7 into a vector of
#' length 1 with 7 columns. Depending on how the samples were summarized, the
#' result can be a single such vector, or a list.
#'
#' For all summarizing operations options \code{sum} and \code{mean} exists. This
#' allows for normalization by the eg. the number of binding sites, size of the
#' range, number of sample, etc..
#'
#' If the resulting object does have a dimension that fits to the number of
#' input ranges the result can be directly attached to them. Basically extending
#' the \code{GRanges} object (\code{out.format}).
#'
#' @param object a BSFDataSet object
#' @param ranges.merge logical; whether to merge ranges
#' @param ranges.merge.method character; how to combine ranges ('sum' or 'mean')
#' @param positions.merge logical; whether to merge positions
#' @param positions.merge.method character; how to combine positions ('sum' or 'mean')
#' @param samples.merge logical: whether to merge samples
#' @param samples.group character; how samples should be grouped when combining
#' ('all', 'condition')
#' @param samples.merge.method charater; how to combine positions ('sum' or 'mean')
#' @param out.format character; how the coverage should be returned
#' ('data.frame' or 'granges'). Note that option 'granges' only exists if the
#' output coverage is of the same rows as the input ranges.
#' @param out.format.overwrite logical; if \code{out.format='granges'}, then
#' decide wheter the meta columns should be extended by the coverage information
#' or be overwritten
#' @param match.rangeID character; unique internal identifier.
#' Name of the meta column of the input ranges
#' that should be used as identifier to match the coverage back to the input
#' ranges. Is 'bsID' as default, since that ID exists for all binding sites
#' after \code{makeBindingSites} was called.
#' @param quiet logical; whether to print messages
#'
#' @return an object of class specified in \code{out.format}
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' bds = makeBindingSites(bds, bsSize = 7)
#' # sum of each replicate over each binding site position
#' c1 = clipCoverage(bds, out.format = "data.frame", positions.merge = TRUE,
#'  ranges.merge = FALSE, samples.merge = FALSE)
#' # total signal per binding site from all samples
#' c2 = clipCoverage(bds, out.format = "granges", positions.merge = FALSE,
#'  ranges.merge = TRUE, samples.merge = TRUE, samples.group = "all")
#' # total signal per binding site from all samples - split by condition
#' c3 = clipCoverage(bds, out.format = "granges", positions.merge = FALSE,
#'  ranges.merge = TRUE, samples.merge = TRUE, samples.group = "condition")
#'
#' @export
clipCoverage <- function(
        object,
        ranges.merge = FALSE,
        ranges.merge.method = c("sum", "mean"),
        positions.merge = FALSE,
        positions.merge.method = c("sum", "mean"),
        samples.merge = TRUE,
        samples.group = c("all", "condition"),
        samples.merge.method = c("sum", "mean"),
        out.format = c("granges", "data.frame"),
        out.format.overwrite = FALSE,
        match.rangeID = "bsID",
        quiet = FALSE
) {

    # init local variables
    # out.format <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is.logical(quiet))

    # check ranges
    rng = getRanges(object)
    if (length(rng) == 0) {
        msg = paste0("No ranges present. Please provide valid input.")
        stop(msg)
    }
    # initialize identical width check
    rng.width.identical = TRUE
    if (length(unique(width(rng))) > 1) {
        msg0 = paste0("Ranges have ", length(unique(width(rng))), " different length.\n")
        msg1 = paste0("Coverage is computed as a single value per range and sample.\n")
        msg2 = paste0("(merge.range = FALSE, merge.sample = FALSE, merge.positions = TRUE)")
        if (!quiet) warning(c(msg0, msg1, msg2))
        rng.width.identical = FALSE
    }

    # match input for positions
    positions.merge.method = match.arg(positions.merge.method, choices = c("sum", "mean"))
    # match input for ranges
    ranges.merge.method = match.arg(ranges.merge.method, choices = c("sum", "mean"))
    # match input for samples
    samples.merge.method = match.arg(samples.merge.method, choices = c("sum", "mean"))
    samples.group = match.arg(samples.group, choices = c("all", "condition"))
    # match output format
    out.format = match.arg(out.format, choices = c("granges", "data.frame"))

    # manage range names
    metaNames = colnames(mcols(rng))
    if (!match.rangeID %in% metaNames) {
        # match id is not present in meta data
        msg0 = paste0(match.rangeID, " was defined as unique ID for the ranges, but is not present in the meta data of the ranges. \n")
        msg1 = paste0("Specify another column or create an appropriate ID. \n")
        stop(c(msg0, msg1))
    } else {
        # match id is present in meta data
        selectID = mcols(rng)[match(match.rangeID, colnames(mcols(rng)))][[1]]
        names(rng) = selectID
    }

    # get signal
    sgn = getSignal(object)

    # get conditions
    conditions = levels(getMeta(object)$condition)

    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    # ranges have all the same size
    if (isTRUE(rng.width.identical)) {
        # compute coverage for all replicates positions and ranges
        # -> this creates list of matrices with the dimensions equal to the
        # --> number of the ranges
        # --> width of the ranges
        cov.final = .makeCoverageMatrix(rng = rng, sgn = sgn)


        # merge output as defined
        cov.final = .mergeSamples(cov.final = cov.final,
                                  samples.merge = samples.merge,
                                  samples.group = samples.group,
                                  samples.merge.method = samples.merge.method,
                                  conditions = conditions)

        cov.final = .mergeRanges(cov.final = cov.final,
                                 ranges.merge = ranges.merge,
                                 ranges.merge.method = ranges.merge.method)

        cov.final = .mergePositions(cov.final = cov.final,
                                    positions.merge = positions.merge,
                                    positions.merge.method = positions.merge.method)
    }
    # ranges have different sizes
    if (!isTRUE(rng.width.identical)) {
        cov.final = .makeCoverageVector(sgn = sgn, rng = rng, method = ranges.merge.method)
    }

    # MANAGE OUTPUT
    # --------------------------------------------------------------------------
    # check output format granges
    if (out.format == "granges") {
        if (is(cov.final, "data.frame")) {
            if (! nrow(cov.final) == length(rng)) {
                msg0 = paste0("Coverage matrix and input ranges do not have the same dimensions. \n")
                msg1 = paste0("Output defaults to 'data.frame' or `vector`.\n")
                if (!quiet) warning(c(msg0, msg1))
                out.format = "data.frame"
            }
        }
        if (is(cov.final, "numeric")) {
            if (!length(cov.final) == length(rng)) {
                msg0 = paste0("Coverage vector and input ranges do not have the same dimensions. \n")
                msg1 = paste0("Output defaults to 'data.frame' or `vector`.\n")
                if (!quiet) warning(c(msg0, msg1))
                out.format = "data.frame"
            }
        }
    }

    # make output
    if (out.format == "granges"){
        rngOut = rng
        if (isTRUE(out.format.overwrite)) {
            mcols(rngOut) = cov.final
        }
        if (!isTRUE(out.format.overwrite)) {
            mcols(rngOut) = cbind(mcols(rngOut), cov.final)
        }
        out = rngOut
    }

    if (out.format == "data.frame") {
        out = cov.final
    }

    return(out)
}

.makeCoverageMatrix <- function(sgn, rng) {

    # split by strand
    rngPlus = rng[strand(rng) == "+"]
    rngMinus = rng[strand(rng) == "-"]

    # Plus strand
    if (length(rngPlus) > 0) {
        # compute initial coverage
        matPlus = lapply(sgn$signalPlus, function(x) {
            y = as.matrix(x[rngPlus])
            rownames(y) = names(rngPlus)
            return(y)
        })
    } else {
        matPlus = NULL
    }

    # Minus strand
    if (length(rngMinus) > 0) {
        # compute initial coverage
        matMinus = lapply(sgn$signalMinus, function(x) {
            y = as.matrix(x[rngMinus])
            rownames(y) = names(rngMinus)
            # reverse minus strand matrices to fit to plus orientation
            y = y[, ncol(y):1]
            return(y)
        })
    } else {
        matMinus = NULL
    }

    # Combine strands
    if (length(rngMinus) > 0 & length(rngPlus) > 0) {
        # if both strands had ranges
        matCov = lapply(seq_along(matPlus), function(x){
            y = dplyr::bind_rows(as.data.frame(matPlus[[x]]), as.data.frame(matMinus[[x]]))
            y = y[match(names(rng), rownames(y)),]
            y
        })
        names(matCov) = names(matPlus)
    }
    if (length(rngMinus) == 0) {
        # only plus strand with ranges
        matCov = lapply(seq_along(matPlus), function(x){
            y = as.data.frame(matPlus[[x]])
            y
        })
        names(matCov) = names(matPlus)
    }
    if (length(rngPlus) == 0) {
        # only minus strand with ranges
        matCov = lapply(seq_along(matMinus), function(x){
            y = as.data.frame(matMinus[[x]])
            y
        })
        names(matCov) = names(matMinus)
    }

    return(matCov)
}

.makeCoverageVector <- function(sgn, rng, method) {

    # split by strand
    rngPlus = rng[strand(rng) == "+"]
    rngMinus = rng[strand(rng) == "-"]

    # Plus strand
    # --------------------------------------------------------------------------
    if (length(rngPlus) > 0) {
        # compute initial coverage
        if (method == "sum") {
            covPlus = as.matrix(
                do.call(cbind, lapply(sgn$signalPlus, function(x) {
                    sum(x[rngPlus])
                })))
            rownames(covPlus) = names(rngPlus)
        }
        if (method == "mean") {
            covPlus = as.matrix(
                do.call(cbind, lapply(sgn$signalPlus, function(x) {
                    mean(x[rngPlus])
                })))
            rownames(covPlus) = names(rngPlus)
        }
    }

    # Minus strand
    # --------------------------------------------------------------------------
    if (length(rngMinus) > 0) {
        # compute initial coverage
        if (method == "sum") {
            covMinus = as.matrix(
                do.call(cbind, lapply(sgn$signalMinus, function(x) {
                    sum(x[rngMinus])
                })))
            rownames(covMinus) = names(rngMinus)
        }
        if (method == "mean") {
            covMinus = as.matrix(
                do.call(cbind, lapply(sgn$signalMinus, function(x) {
                    mean(x[rngMinus])
                })))
            rownames(covMinus) = names(rngMinus)
        }
    }

    # Combine strands for return
    # --------------------------------------------------------------------------
    if (length(rngMinus) > 0 & length(rngPlus) > 0) {
        cov.final = rbind(covPlus, covMinus)
    }
    if (length(rngMinus) == 0) {
        cov.final = rngPlus
    }
    if (length(rngPlus) == 0) {
        cov.final = rngMinus
    }
    # match results to manage order
    idx = match(rownames(cov.final), names(rng))
    cov.final = cov.final[order(idx),]

    return(cov.final)
}

.mergeRanges <- function(cov.final,
                         ranges.merge,
                         ranges.merge.method
){
    # merge ranges
    if (isTRUE(ranges.merge)) {
        if (is(cov.final, "list")) {
            # merge a list of matrices
            if (ranges.merge.method == "sum") {
                cov.final = lapply(cov.final, rowSums, na.rm = TRUE)
            }
            if (ranges.merge.method == "mean") {
                cov.final = lapply(cov.final, rowMeans, na.rm = TRUE)
            }
        }
        if (is(cov.final, "data.frame")) {
            # merge a single matrix
            if (ranges.merge.method == "sum") {
                cov.final = rowSums(cov.final, na.rm = TRUE)
            }
            if (ranges.merge.method == "mean") {
                cov.final = rowMeans(cov.final, na.rm = TRUE)
            }
        }

    }
    return(cov.final)
}

.mergeSamples <- function(cov.final,
                          samples.merge,
                          samples.group,
                          samples.merge.method,
                          conditions
){
    if (isTRUE(samples.merge)) {
        if (samples.group == "condition") {
            # merge samples by condition
            # -> return list of matrices with length equal to number of conditions
            # -> should work with NAs
            mat.list = lapply(conditions, function(x){cov.final[c(grep(x, names(cov.final)))]})
            names(mat.list) = conditions

            if (samples.merge.method == "sum") {
                # merge by sum
                cov.final = lapply(mat.list, function(x){Reduce('+', x)})
            }
            if (samples.merge.method == "mean") {
                # merge by mean
                cov.final = lapply(mat.list, function(x){Reduce('+', x) / length(x)})
            }
        }
        if (samples.group == "all") {
            # merge all samples
            # -> returns a single matrix
            if (samples.merge.method == "sum") {
                # merge by sum
                cov.final = Reduce("+", cov.final)
            }
            if (samples.merge.method == "mean") {
                # merge by mean
                cov.final = Reduce("+", cov.final) / length(cov.final)
            }
        }

    }
    return(cov.final)
}

.mergePositions <- function(cov.final,
                            positions.merge,
                            positions.merge.method
){
    # merge positions
    if (isTRUE(positions.merge)) {
        if (is(cov.final, "numeric")) {
            # merge a single vector
            # -> from merge ranges
            if (positions.merge.method == "sum") {
                cov.final = sum(cov.final)
            }
            if (positions.merge.method == "mean") {
                cov.final = mean(cov.final)
            }
        }
        if (is(cov.final, "data.frame")) {
            # merge a single matrix
            # -> from samples option 'all'
            if (positions.merge.method == "sum") {
                cov.final = colSums(cov.final)
            }
            if (positions.merge.method == "mean") {
                cov.final = colMeans(cov.final)
            }
        }
        if (is(cov.final, "list")) {
            # merge a list
            if (all(unlist(lapply(cov.final, class)) == "numeric")) {
                # the list is made from vectors
                if (positions.merge.method == "sum") {
                    cov.final = lapply(cov.final, sum)
                    cov.final = as.data.frame(do.call(rbind, cov.final))
                }
                if (positions.merge.method == "mean") {
                    cov.final = lapply(cov.final, mean)
                    cov.final = as.data.frame(do.call(rbind, cov.final))
                }
            }
            if (all(unlist(lapply(cov.final, class)) == "data.frame")) {
                # the list is made from matrices
                if (positions.merge.method == "sum") {
                    cov.final = lapply(cov.final, colSums)
                    cov.final = as.data.frame(do.call(rbind, cov.final))
                }
                if (positions.merge.method == "mean") {
                    cov.final = lapply(cov.final, colMeans)
                    cov.final = as.data.frame(do.call(rbind, cov.final))
                }
            }
        }
    }
    return(cov.final)
}



#-------------------------------------------------------------------------------
# OLD Function
#-------------------------------------------------------------------------------

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
#' @param method sum/ mean, select how replicates/ ranges should be summarized
#' @param allowNA TRUE/ FALSE, allow NA values in the output if input ranges
#' are of different width
#' @param quiet logical, whether to print messages
#'
#' @return an object of class specified in \code{returnOptions}
#' @import GenomicRanges
#' @importFrom plyr rbind.fill
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' rng = coverageOverRanges(bds, returnOptions = "merge_ranges_keep_positions")
#' rng = coverageOverRanges(bds, returnOptions = "merge_replicates_per_condition")
#' rng = coverageOverRanges(bds, returnOptions = "merge_all_replicates")
#' rng = coverageOverRanges(bds, returnOptions = "merge_positions_keep_replicates")
#'
#' @export
coverageOverRanges <- function(
    object,
    returnOptions = c("merge_ranges_keep_positions",
                      "merge_replicates_per_condition",
                      "merge_all_replicates",
                      "merge_positions_keep_replicates"),
    method = "sum",
    allowNA = FALSE,
    quiet = TRUE
    ) {

    # deprecation notice
    # .Deprecated("clipCoverage")
    # first exchange everywhere it is used

    # Check input
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))

    # manage return option
    returnOptions = match.arg(
        returnOptions,
        choices = c("merge_ranges_keep_positions",
                    "merge_replicates_per_condition",
                    "merge_all_replicates",
                    "merge_positions_keep_replicates"))
    method = match.arg(method,choices = c("sum", "mean"))
    # get range
    rng = getRanges(object)

    # check if lenght is all the same or different
    if (length(unique(width(getRanges(object)))) > 1 &
        allowNA == FALSE &
        returnOptions != "merge_positions_keep_replicates") {
        # ranges do not have all the same size
        warning("Width of ranges are not identical. This will cause NAs in the
                output. Option `merge_positions_keep_replicates` is enforced.
                Avoid this behavior by setting allowNA = TRUE.")
        returnOptions = "merge_positions_keep_replicates"
    }

    # check for unique range names
    if(any(duplicated(names(rng)))) {
        msg0 = paste0("duplicate names in ranges, adding preceding ID")
        if(!quiet) message(c(msg0))
        names(rng) = paste0(c(seq_along(rng)), "_", names(rng))
    }

    # prepare signal
    sgn = getSignal(object)
    # extract present conditions
    condition = levels(getMeta(object)$condition)


    # Forward to helper functions
    # --------------------------------------------------------------------------
    if (returnOptions == "merge_replicates_per_condition") {
        covRet = .coverageOverRanges.merge_replicates_per_condition(
            sgn = sgn, rng = rng, condition = condition, method = method, allowNA = allowNA
        )
    }
    if (returnOptions == "merge_all_replicates") {
        covRet = .coverageOverRanges.merge_all_replicates(
            sgn = sgn, rng = rng, method = method, allowNA = allowNA
        )
    }
    if (returnOptions == "merge_positions_keep_replicates") {
        covRet = .coverageOverRanges.merge_positions_keep_replicates(
            sgn = sgn, rng = rng, method = method
        )
    }
    if (returnOptions == "merge_ranges_keep_positions") {
        covRet = .coverageOverRanges.merge_ranges_keep_positions(
            sgn = sgn, rng = rng, method = method
        )
    }
    return(covRet)
}

.coverageOverRanges.merge_positions_keep_replicates <- function(sgn, rng, method) {

    # split by strand
    rngPlus = rng[strand(rng) == "+"]
    rngMinus = rng[strand(rng) == "-"]

    # Plus strand
    # --------------------------------------------------------------------------
    if (length(rngPlus) > 0) {
        # compute initial coverage
        if (method == "sum") {
            mcols(rngPlus) = as.matrix(
                do.call(cbind, lapply(sgn$signalPlus, function(x) {
                    sum(x[rngPlus])
                })))
        }
        if (method == "mean") {
            mcols(rngPlus) = as.matrix(
                do.call(cbind, lapply(sgn$signalPlus, function(x) {
                    mean(x[rngPlus])
                })))
        }

    }

    # Minus strand
    # --------------------------------------------------------------------------
    if (length(rngMinus) > 0) {
        # compute initial coverage
        if (method == "sum") {
            mcols(rngMinus) = as.matrix(
                do.call(cbind, lapply(sgn$signalMinus, function(x) {
                    sum(x[rngMinus])
                })))
        }
        if (method == "mean") {
            mcols(rngMinus) = as.matrix(
                do.call(cbind, lapply(sgn$signalMinus, function(x) {
                    mean(x[rngMinus])
                })))
        }
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
    rngCov = .sortRanges(rngCov)

    return(rngCov)
}

.coverageOverRanges.merge_ranges_keep_positions <- function(sgn, rng, method) {
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
        if (method == "sum") {
            covPlus = do.call(rbind, lapply(matPlus, colSums))
        }
        if (method == "mean") {
            covPlus = do.call(rbind, lapply(matPlus, colMeans))
        }
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
        if (method == "sum") {
            covMinus = do.call(rbind, lapply(matMinus, colSums))
        }
        if (method == "mean") {
            covMinus = do.call(rbind, lapply(matMinus, colMeans))
        }
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

.coverageOverRanges.merge_all_replicates <- function(sgn, rng, method, allowNA) {
    # set names
    names(rng) = seq_along(rng)
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
        if (method == "sum") {
            covPlus = Reduce('+', matPlus)
        }
        if (method == "mean") {
            covPlus = Reduce('+', matPlus) / length(matPlus)
        }
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
        if (method == "sum") {
            covMinus = Reduce('+', matMinus)
        }
        if (method == "mean") {
            covMinus = Reduce('+', matMinus) / length(matMinus)
        }
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
        if (allowNA == TRUE) {
            retCov = plyr::rbind.fill(as.data.frame(covPlus), as.data.frame(covMinus))
        }
        if (allowNA == FALSE) {
            retCov = rbind(covPlus, covMinus)
        }
    }
    if (length(rngMinus) == 0) {
        retCov = covPlus
    }
    if (length(rngPlus) == 0) {
        retCov = covMinus
    }
    # match results to manage order
    idx = match(rownames(retCov), names(rng))
    retCov = retCov[order(idx),]
    return(retCov)
}

.coverageOverRanges.merge_replicates_per_condition <- function(sgn, rng, condition, method, allowNA) {
    # set names
    names(rng) = seq_along(rng)
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
        # manage combine method
        if (method == "sum") {
            covPlus = lapply(matPlusCond, function(x){Reduce('+', x)})
        }
        if (method == "mean") {
            covPlus = lapply(matPlusCond, function(x){Reduce('+', x) / length(x)})
        }
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
        # mange combine method
        if (method == "sum") {
            covMinus = lapply(matMinusCond, function(x){Reduce('+', x)})
        }
        if (method == "mean") {
            covMinus = lapply(matMinusCond, function(x){Reduce('+', x) / length(x)})
        }
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
        if (allowNA == TRUE) {
            retCov = lapply(seq_along(covPlus), function(x){
                retCov = plyr::rbind.fill(as.data.frame(covPlus[[x]]), as.data.frame(covMinus[[x]]))
            })
        }
        if (allowNA == FALSE) {
            retCov = lapply(seq_along(covPlus), function(x){
                rbind(covPlus[[x]],covMinus[[x]])
            })
        }
        names(retCov) = names(covPlus)
    }
    if (length(rngMinus) == 0) {
        retCov = covPlus
    }
    if (length(rngPlus) == 0) {
        retCov = covMinus
    }
    # match results to manage order
    retCov = lapply(retCov, function(x){
        idx = match(rownames(x), names(rng))
        x = x[order(idx),]
    })
}

