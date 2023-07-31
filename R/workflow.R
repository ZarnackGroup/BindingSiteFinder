################################################################################
# Workflow functions
# - these functions are used in the main bindingSiteFinde workflow
# - each function updates the bindingSiteFinder object
# - functions are either designed to work on single nt wide crosslink sites or
# on merged binding sites
# - for each function a related plotting funciton exists to visualize the results
# - the wrapper function BSFind() executes the functions below in the desired order
################################################################################


#' Filter PureCLIP sites by their score distribution
#'
#' Function that applies a filter on the global crosslink site score distribution.
#' The \code{\link{GenomicRanges}} contained in the \code{\link{BSFDataSet}} need to
#' have a meta-column that holds a numeric score value, which is used for filtering.
#' The name of the column can be set with \code{match.score}.
#'
#' @param object a \code{\link{BSFDataSet}} object with stored crosslink ranges
#' of width=1
#' @param cutoff numeric; defines the cutoff for which sites to keep, the
#' smallest step is 1\% (0.01)
#' @param match.score character; meta column name of the crosslink site
#' \code{\link{GenomicRanges}} object that holds the score which is used for
#' sub-setting
#' @param quiet logical; whether to print messages
#'
#' @return an object of class  \code{\link{BSFDataSet}} with its ranges filtered
#' by those that passed the threshold set with \code{cutoff}
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' # apply 5% filter
#' pureClipGlobalFilter(object = bds, cutoff = 0.05)
#'
#' @export
pureClipGlobalFilter <- function(object, # bindingSiteFinder
                                 cutoff = 0.01, # defines the cutoff for which PureCLIP sites to keep; smallest are 1% (0.01); a cutoff of 0.05 means to remove the lowest 5% PureCLIP sites base on the score
                                 match.score = "score",
                                 quiet = FALSE
) {
    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is.logical(quiet))
    # check for correct cutoff value
    possibleCutoffs = seq(from = 0, to = 1, by = 0.01)
    if (! cutoff %in% possibleCutoffs) {
        msg = paste0("Specified cutoff (", cutoff,
                     "), does not refer to a percentage value. Try cutoffs of the form 0.01, 0.05, ... instead.")
        stop(msg)
    }
    # Check for correct score column
    metaNames = colnames(mcols(getRanges(object)))
    if (! match.score %in% metaNames) {
        msg = paste0("Specified match.score (", match.score, "), is present in the ranges. Add respective column. ")
        stop(msg)
    }

    # Create matching vectors for score column from input source
    # --------------------------------------------------------------------------
    rngInitial = getRanges(object)
    score = mcols(rngInitial)[match(match.score, colnames(mcols(rngInitial)))][[1]]

    # ---
    # Store function parameters in list
    optstr = list(cutoff = cutoff, match.score = match.score)
    object@params$pureClipGlobalFilter = optstr

    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    # compute quantile based on PureCLIP score
    qCut = quantile(score, probs = seq(0,1, by = 0.01))

    # apply user set cutoff
    names(qCut) = as.numeric(sub("%", "", names(qCut))) / 100
    idx = match(cutoff, names(qCut))
    rngFiltered = rngInitial[score >= qCut[idx]]

    # ---
    # Store for plotting
    dfPlot = score
    cutoffPoint = qCut[idx]
    object@plotData$pureClipGlobalFilter$data = dfPlot
    object@plotData$pureClipGlobalFilter$cutoffPoint = cutoffPoint

    # ---
    # Store for results
    resultLine = data.frame(
        funName = "pureClipGlobalFilter()", class = "crosslink sites",
        nIn = length(rngInitial), nOut = length(rngFiltered),
        per = paste0(round(length(rngFiltered)/ length(rngInitial), digits = 2)*100,"%"),
        options = paste0("Cutoff=", optstr$cutoff, ", scoreMatchCol=", optstr$match.score)
    )
    object@results = rbind(object@results, resultLine)

    # return BSF object with only ranges above cutoff
    objectNew = setRanges(object, rngFiltered, quiet = quiet)
    return(objectNew)
}


#' Filter PureCLIP sites by their score distribution per gene
#'
#' Function that applies a filter on the crosslink site score distribution at
#' gene level. This allows to filter for those sites with the strongest signal
#' on each gene. Since scores are tied to the expression level of the hosting
#' transcript this function allows a fair filter for all genes partially
#' independet of the expression level.
#'
#' The \code{\link{GenomicRanges}} contained in the \code{\link{BSFDataSet}} need to
#' have a meta-column that holds a numeric score value, which is used for filtering.
#' The name of the column can be set with \code{scoreCol}.
#'
#' In the case of overlapping gene annotation, a single crosslink site will be
#' attributed to multiple genes. The \code{\link{overlaps}} parameter allows
#' to control these cases. Option `keepSingle` will only keep a single instance
#' of the site; `removeAll` will remove both sites; `keepAll` will keep both
#' sites.
#'
#' @param object a \code{\link{BSFDataSet}} object with stored crosslink ranges
#' of width=1
#' @param cutoff numeric; defines the cutoff for which sites to remove, the
#' smallest step is 1\% (0.01). A cutoff of 5\% will remove the lowest 5\% sites,
#' given their score, on each gene, thus keeping the strongest 95\%.
#' @param overlaps character; how overlapping gene loci should be handled.
#' @param anno.annoDB an object of class \code{OrganismDbi} that contains
#' the gene annotation.
#' @param anno.genes an object of class \code{\link{GenomicRanges}} that represents
#' the gene ranges directly
#' @param match.score character; meta column name of the crosslink site
#' \code{\link{GenomicRanges}} object that holds the score which is used for
#' sub-setting
#' @param match.geneID character; meta column name of the genes
#' \code{\link{GenomicRanges}} object that holds a unique geneID
#' @param quiet logical; whether to print messages
#'
#' @return an object of class \code{\link{BSFDataSet}} with its ranges filtered
#' by those that passed the gene-wise threshold set with \code{cutoff}
#'
#' @importFrom dplyr filter reframe group_by mutate mutate_all arrange row_number
#' @importFrom GenomicFeatures genes
#' @importFrom stats density
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' # Load GRanges with genes
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#' # apply 5% gene-wise filter
#' pureClipGeneWiseFilter(object = bds, anno.genes = gns, cutoff = 0.5, overlaps = "keepSingle")
#'
#' @export
pureClipGeneWiseFilter <- function(object, # bindingSiteFinder
                                   cutoff = 0.05, # defines the cutoff for which PureCLIP sites to keep; smallest steps are in 1% (0.01); cutoff of 0.7 means that the top 70% of PureCLIP sites per gene are retained.
                                   overlaps = c("keepSingle", "removeAll", "keepAll"), # Options to deal with PureCLIP sites on overlapping loci genes; use keepSingle as default; keepAll will inflate the number of binding sites by the number of overlaps, removeAll will kill all duplicated instances
                                   anno.annoDB = NULL, # annotation data base of class OrganismDbi; can be NULL -> but then anno.genes must be provided
                                   anno.genes = NULL, # ranges of the genes of class GenomicRanges; can be NULL -> but then anno.annoDB must be provided
                                   match.score = "score",
                                   match.geneID = "gene_id",
                                   quiet = FALSE) {
    # Local variables
    datasource <- cutoffs <- quantiles <- name <- NULL
    # INPUT CHECKS
    # --------------------------------------------------------------------------
    # type checks
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is.logical(quiet))

    # Check if none is specified
    if (is.null(anno.annoDB) & is.null(anno.genes)) {
        msg = paste0("None of the required annotation sources anno.annoDB or anno.genes was specified. ")
        stop(msg)
    }
    # Check if both are specified
    if (!is.null(anno.annoDB) & !is.null(anno.genes)) {
        msg = paste0("Both of the required annotation sources anno.annoDB or anno.genes are specified. Please provide only one of the two. ")
        stop(msg)
    }
    # Checks if anno.annoDB should be used
    if (!is.null(anno.annoDB) & is.null(anno.genes)) {
        stopifnot(is(anno.annoDB, "OrganismDb"))
        datasource = "anno.annoDB"
        # extract relevant annotation
        anno.genes = genes(anno.annoDB, columns = c("GENEID"))
        # Create matching vectors for columns from input annotation
        selectID = as.character(anno.genes$GENEID)
    }
    # Checks if anno.genes should be used
    if (is.null(anno.annoDB) & !is.null(anno.genes)) {
        stopifnot(is(anno.genes, "GenomicRanges"))

        datasource = "anno.genes"
        # extract relevant annotation
        anno.genes = anno.genes
        # check correct annotation columns
        annoColNames = colnames(mcols(anno.genes))
        if (!all(match.geneID %in% annoColNames)) {
            msg = paste0("The provided matching column for the geneID (",
                         match.geneID, ") is not present in the provided annotation. \n")
            stop(msg)
        }
        # Create matching vectors for columns from input annotation
        selectID = mcols(anno.genes)[match(match.geneID, colnames(mcols(anno.genes)))][[1]]

    }
    # check if cutoff was estimated before or if default should be used
    if (!is.null(object@params$bsSize) &
        !is.null(object@params$geneFilter)) {
        cutoff = object@params$geneFilter
    }

    # handle options (remove, keep)
    overlaps = match.arg(overlaps, choices = c("keepSingle", "removeAll", "keepAll"))
    # Check for correct score column
    # get all ranges
    rngInitial = getRanges(object)
    metaNames = colnames(mcols(rngInitial))
    if (! match.score %in% metaNames) {
    # if (! any(metaNames == match.score)) {
        msg = paste0("Specified match.score (", match.score, "), is present in the ranges. Add respective column. ")
        stop(msg)
    }

    # ---
    # Store function parameters in list
    optstr = list(cutoff = cutoff, overlaps = overlaps,
                  source = datasource, match.geneID = match.geneID)
    object@params$geneWiseFilter = optstr

    # MAIN COMPUTE
    # --------------------------------------------------------------------------

    # Create matching vectors for score column from input source
    # --------------------------------------------------------------------------
    score = mcols(rngInitial)[match(match.score, colnames(mcols(rngInitial)))][[1]]
    rngInitial$score = score
    # get all intergenic cases
    rngIntergeneic = rngInitial[countOverlaps(rngInitial, anno.genes) == 0]
    mcols(rngIntergeneic)$geneID = "Intergenic"
    # get all cases to solve
    # rngToSolve = subset(rngInitial, !rngInitial %in% rngIntergeneic)
    # rngToSolve = rngInitial[!rngInitial %in% rngIntergeneic]
    rngToSolve = rngInitial[countOverlaps(rngInitial, rngIntergeneic) == 0]

    # Gene specific cutoffs
    # ---------------------
    # set cutpoints
    potentialCutoffs = seq(0, 1, by = 0.01)
    # group scores by gene range
    ols = findOverlaps(anno.genes, rngToSolve) %>%
        as.data.frame() %>%
        mutate(score = rngToSolve$score[subjectHits]) %>%
        group_by(queryHits)
    # calculate gene range specific quantile cutoff
    quants = ols %>%
        dplyr::reframe(quantiles = quantile(score, probs = potentialCutoffs)) %>%
        group_by(queryHits) %>%
        mutate(cutoffs = potentialCutoffs) %>%
        dplyr::filter(as.character(cutoffs) == as.character(cutoff))
    # apply cutoff to every PureCLIP peak per gene
    peakIdToKeep = ols
    peakIdToKeep$quantiles = quants$quantiles[match(peakIdToKeep$queryHits, quants$queryHits)]
    peakIdToKeep = peakIdToKeep %>%
        dplyr::filter(score >= quantiles)
    # add gene information
    filteredPerRegion = rngToSolve[peakIdToKeep$subjectHits]
    mcols(filteredPerRegion)$geneID = selectID[peakIdToKeep$queryHits]
    # ---
    # Store for plotting
    dfPlot = countOverlaps(filteredPerRegion) %>% table() %>% as.data.frame() %>% rename("#N overlaps" = ".")
    object@plotData$geneWiseFilter$data = dfPlot

    # Deal with multiple loci
    # -----------------------
    # check for overlapping loci peaks
    potentialOverlaps = findOverlaps(filteredPerRegion)
    overlappingLociHits = potentialOverlaps[queryHits(potentialOverlaps) != subjectHits(potentialOverlaps)]
    # reduce to a single pair
    singlePair = overlappingLociHits %>%
        as.data.frame() %>%
        dplyr::filter(queryHits > subjectHits)
    overlappingLociHits = overlappingLociHits[queryHits(overlappingLociHits) %in% singlePair$queryHits]
    # get numbers
    overlappingLociGenePairs = paste0(unique(names(filteredPerRegion[queryHits(overlappingLociHits)])),
                                      ",",
                                      unique(names(filteredPerRegion[subjectHits(overlappingLociHits)])), " | ")
    fractionOfRangesRemoved = paste0(round(length(overlappingLociHits) / length(rngInitial) * 100, digits = 2), "%")
    # Handle cases
    if (length(overlappingLociHits) > 0) {
        msg0 = paste0(fractionOfRangesRemoved, " (", length(overlappingLociHits),"/", length(rngInitial), ")",
                      " peaks overlap with multiple anno.genes in the given gene annotation. \n")
        if (overlaps == "keepSingle") {
            msg1 = "A single instance of each peak is kept. This is recommended. \n "
            filteredPerRegion = filteredPerRegion[!duplicated(filteredPerRegion)]
            # report message
            if (!quiet) message(c(msg0, msg1))
        }
        if (overlaps == "removeAll") {
            msg1 = "All peaks on any duplicated range are removed. This is not recommended, since many peaks are lost. \n"
            filteredPerRegion = filteredPerRegion[countOverlaps(filteredPerRegion) == 1]
            warning(c(msg0, msg1))
        }
        if (overlaps == "keepAll") {
            msg1 = "Duplicated peaks due to overlapping loci are not removed. This is not recommended, since peak numbers are inflated. \n"
            warning(msg0, msg1)
        }
    }

    # add intergenic cases again
    filteredPerRegion = c(filteredPerRegion, rngIntergeneic)

    # ---
    # Store for results
    resultLine = data.frame(
        funName = "pureClipGeneWiseFilter()", class = "crosslink sites",
        nIn = length(rngInitial), nOut = length(filteredPerRegion),
        per = paste0(round(length(filteredPerRegion)/ length(rngInitial), digits = 2)*100,"%"),
        options = paste0("Cutoff=", optstr$cutoff, ", overlapp option=", optstr$overlaps)
    )
    object@results = rbind(object@results, resultLine)

    objectOut = setRanges(object, filteredPerRegion, quiet = quiet)
    return(objectOut)
}


#' Assign binding sites to their hosting genes
#'
#' Function that assigns each binding site in the \code{\link{BSFDataSet}} to its
#' hosting gene given a gene annotation (\code{anno.annoDB}, \code{anno.genes}).
#'
#' Regardless of the annotation source that is being used, the respective meta
#' information about the genes have to be present. They can be set by the
#' \code{match.geneID}, \code{match.geneName} and \code{match.geneType} arguments.
#'
#' In the case of overlapping gene annotation, a single binding site will be
#' associated with multiple genes. The \code{\link{overlaps}} parameter allows
#' to decide in these cases. Option `frequency` will take the most frequently
#' observed gene type, option `hierarchy` works in conjunction with a user defined
#' rule (\code{overlaps.rule}). Options `remove` and `keep` will remove or
#' keep all overlapping cases, respectively.
#'
#' Note that if an overlaps exists, but gene types are identical options
#' `frequency` and `hierarchy` will cause the gene that was seen first to be
#' selected as representative.
#'
#' @param object a \code{\link{BSFDataSet}} object with stored binding sites. This
#' means that ranges should be > 1
#' @param overlaps character; how overlapping gene loci should be handled.
#' @param overlaps.rule character vector; a vector of gene type that should
#' be used to handle overlapping cases in a hierarchical manor. The order of the
#' vector is the order of the hierarchy.
#' @param anno.annoDB an object of class \code{OrganismDbi} that contains
#' the gene annotation.
#' @param anno.genes an object of class \code{\link{GenomicRanges}} that represents
#' the gene ranges directly
#' @param match.geneID character; meta column name of the gene ID
#' @param match.geneName character; meta column name of the gene name
#' @param match.geneType character; meta column name of the gene type
#' @param quiet logical; whether to print messages
#'
#' @return an object of class \code{\link{BSFDataSet}} with binding sites having
#' hosting gene information added to their meta columns.
#'
#' @importFrom dplyr select slice_head group_by count ungroup arrange desc left_join
#' @importFrom GenomicFeatures genes
#' @importFrom dplyr desc
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' # Load GRanges with genes
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#' bds = makeBindingSites(object = bds, bsSize = 9)
#' bds = assignToGenes(bds, anno.genes = gns)
#'
#' @export
assignToGenes <- function(object,
                          overlaps = c("frequency", "hierarchy", "remove", "keep"),
                          overlaps.rule = NULL,
                          anno.annoDB = NULL,
                          anno.genes = NULL,
                          match.geneID = "gene_id",
                          match.geneName = "gene_name",
                          match.geneType = "gene_type",
                          quiet = FALSE
) {
    # local variables
    datasource <- geneIndex <- geneID <- geneName <- value <- bsIndex <- Freq <- geneType <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is.logical(quiet))

    # Check if none is specified
    if (is.null(anno.annoDB) & is.null(anno.genes)) {
        msg = paste0("None of the required annotation sources anno.annoDB or anno.genes was specified. ")
        stop(msg)
    }
    # Check if both are specified
    if (!is.null(anno.annoDB) & !is.null(anno.genes)) {
        msg = paste0("Both of the required annotation sources anno.annoDB or anno.genes are specified. Please provide only one of the two. ")
        stop(msg)
    }
    # Checks if anno.annoDB should be used
    if (!is.null(anno.annoDB) & is.null(anno.genes)) {
        stopifnot(is(anno.annoDB, "OrganismDb"))
        if (!is.null(anno.genes)) {
            msg = paste0("Parameter anno.annoDB and anno.genes are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            datasource = "anno.annoDB"
            # extract relevant annotation
            anno.genes = genes(anno.annoDB, columns=c("ENSEMBL", "GENETYPE", "SYMBOL", "GENEID"))
            # Create matching vectors for columns from input annotation
            # --------------------------------------------------------------------------
            selectID = as.character(anno.genes$GENEID)
            selectName = as.character(anno.genes$SYMBOL)
            selectType = as.character(anno.genes$GENETYPE)
        }
    }
    # Checks if anno.genes should be used
    if (is.null(anno.annoDB) & !is.null(anno.genes)) {
        stopifnot(is(anno.genes, "GenomicRanges"))
        if (!is.null(anno.annoDB)) {
            msg = paste0("Parameter anno.annoDB and anno.genes are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            datasource = "anno.genes"
            # extract relevant annotation
            anno.genes = anno.genes
            # check correct annotation columns
            inNames = c(match.geneID, match.geneName, match.geneType)
            annoColNames = colnames(mcols(anno.genes))
            if (!all(inNames %in% annoColNames)) {
                msg = paste0("One or multiple of the provided matching columns (",
                             match.geneID, ", ", match.geneName, ", ", match.geneType,
                             ") is not present in the provided annotation. \n")
                stop(msg)
            }
            # Create matching vectors for columns from input annotation
            # --------------------------------------------------------------------------
            selectID = mcols(anno.genes)[match(match.geneID, colnames(mcols(anno.genes)))][[1]]
            selectName = mcols(anno.genes)[match(match.geneName, colnames(mcols(anno.genes)))][[1]]
            selectType = mcols(anno.genes)[match(match.geneType, colnames(mcols(anno.genes)))][[1]]

        }
    }

    # handle options (hierarchy, frequency, remove, keep)
    overlaps = match.arg(overlaps, choices = c("frequency", "hierarchy", "remove", "keep"))
    # Check multiple loci options
    if (overlaps == "hierarchy" & is.null(overlaps.rule)) {
        msg1 = paste0("Binding sites on anno.genes with overlapping annotaitons is set to be handled by 'hierarchy', but no rule is provided. \n")
        msg2 = paste0("Change the 'overlaps' option or provide a valid rule. \n")
        stop(c(msg1, msg2))
    }
    # --------------------------------------------------------------------------
    # MAIN
    # --------------------------------------------------------------------------

    # ---
    # Store function parameters in list
    optstr = list(source = datasource, match.geneID = match.geneID,
                  match.geneName = match.geneName, match.geneType = match.geneType,
                  overlaps = overlaps, overlaps.rule = overlaps.rule)
    object@params$assignToGenes = optstr

    # Match Binding sites and anno.genes
    # --------------------------------------------------------------------------
    # split ranges in duplicated event and non duplicated events
    rngInitial = getRanges(object)
    rngInitial$currIdx = seq_along(rngInitial)
    ols = findOverlaps(anno.genes, rngInitial)
    # Deal with multiple loci
    # --------------------------------------------------------------------------
    # Count out how many cases there are
    countOlsSize = as.data.frame(table(duplicated(subjectHits(ols))))
    totalBS = countOlsSize$Freq[1]
    duplicatedBS = countOlsSize$Freq[2]
    duplicatedFraction = paste0(round(duplicatedBS / totalBS * 100, digits = 2), "%")
    # ----
    # Store results for plotting
    dfPlot = data.frame(geneIndex = queryHits(ols), bsIndex = subjectHits(ols),
                        geneID = selectID[queryHits(ols)],
                        geneName = selectName[queryHits(ols)],
                        geneType = selectType[queryHits(ols)]) %>%
        mutate(value = 1) %>%
        dplyr::select(-geneIndex, -geneID, -geneName) %>%
        pivot_wider(names_from = geneType, values_from = value, values_fn = length, values_fill = 0) %>%
        dplyr::select(-bsIndex) %>%
        mutate_all(~ ifelse(. > 1, 1, .))
    object@plotData$assignToGenes$dataOverlaps = dfPlot
    # Case where there are no overlaps
    if (duplicatedBS == 0) {
        msg = paste0("No binding sites (",
                     format(totalBS, big.mark = ',', decimal.mark = "."),
                     ") on overlapping gene loci. \n")
        msg2 = paste0("Parameter 'overlaps' set to '",
                      overlaps, "' not in use. \n")
        if(!quiet) message(c(msg, msg2))
    }
    # Cases where there are multiple overlaps
    if (length(duplicatedBS) > 0) {
        msg0 = paste0(duplicatedFraction, " (", duplicatedBS,"/", totalBS, ")",
                      " of binding sites overlap with multiple anno.genes in the given gene annotation. \n")
        if (overlaps == "hierarchy") {
            # apply hierarchy solution
            msg1 = paste0("Apply 'hierarchy' solution")
            if(!quiet) message(c(msg0, msg1))
            ruleMod = data.frame(gene_type = overlaps.rule, idx = seq_along(overlaps.rule))
            rngResolved = .resolveGeneOverlapsWithRule(rng = rngInitial, ols = ols, rule = ruleMod,
                                                       selectID = selectID, selectName = selectName,
                                                       selectType = selectType)
            if (length(rngInitial) != length(rngResolved)) {
                if(!quiet) message(paste0(format(length(rngInitial)-length(rngResolved), big.mark = ",", decimal.mark = "."),
                                          " binding sites could not be assigned to a gene range.\n"))
            }

        }
        if (overlaps == "frequency") {
            # apply frequency solution
            msg1 = paste0("Apply 'frequency' solution")
            if (!quiet) message(c(msg0, msg1))
            ruleFreq = selectType %>% table() %>% as.data.frame() %>%
                rename(gene_type = 1) %>% arrange(desc(Freq)) %>%
                mutate(idx = row_number())  %>% rename(n = 2)
            rngResolved = .resolveGeneOverlapsWithRule(rng = rngInitial, ols = ols, rule = ruleFreq,
                                                       selectID = selectID, selectName = selectName, selectType = selectType)
            if (length(rngInitial) != length(rngResolved)) {
                if (!quiet) message(paste0(format(length(rngInitial)-length(rngResolved), big.mark = ",", decimal.mark = "."),
                                           " binding sites could not be assigned to a gene range.\n"))
            }

        }
        if (overlaps == "remove") {
            # apply remove option
            msg1 = "Binding sites from overlapping loci are removed. This is not recommended. Please see options 'heirarchy' and 'frequency'. \n"
            warning(c(msg0, msg1))
            rngResolved = rngInitial[subjectHits(ols)]
            mcols(rngResolved) = cbind(mcols(rngResolved), geneID = selectID[queryHits(ols)],
                                       geneType = selectType[queryHits(ols)], geneName = selectName[queryHits(ols)] )
            rngResolved = rngResolved[countOverlaps(rngResolved) == 1]
        }
        if (overlaps == "keep") {
            msg1 = "Binding sites from overlapping loci are kept. This is not recommended. Please see options 'heirarchy' and 'frequency'. \n"
            warning(c(msg0, msg1))
            rngResolved = rngInitial[subjectHits(ols)]
            mcols(rngResolved) = cbind(mcols(rngResolved), geneID = selectID[queryHits(ols)],
                                       geneType = selectType[queryHits(ols)], geneName = selectName[queryHits(ols)] )
        }
    }
    # ---
    # Store results for plotting
    nCountGenes = rngResolved %>% mcols() %>% as.data.frame() %>%
        group_by(geneID) %>% count(geneType) %>% ungroup() %>%
        count(geneType) %>% arrange(desc(n))
    nCountBs = rngResolved %>% mcols() %>% as.data.frame() %>%
        count(geneType) %>% arrange(desc(n))
    nCountAll = left_join(nCountBs, nCountGenes, by = "geneType") %>%
        rename("nBs" = "n.x", "nGenes" = "n.y")

    object@plotData$assignToGenes$dataSpectrum = nCountAll

    # ---
    # Store for results
    resultLine = data.frame(
        funName = "assignToGenes()", class = "binding sites",
        nIn = length(rngInitial), nOut = length(rngResolved),
        per = paste0(round(length(rngResolved)/ length(rngInitial), digits = 2)*100,"%"),
        options = paste0("Overlaps=", optstr$overlaps, ", Source=", optstr$source,
                         ifelse(!rlang::is_empty(optstr$overlaps.rule),
                                paste0(", overlaps.rule=",
                                       paste(optstr$overlaps.rule, collapse = ">")), ""))
    )
    object@results = rbind(object@results, resultLine)

    objectOut = setRanges(object, rngResolved, quiet = quiet)
    return(objectOut)
}


#' Assign binding sites to their hosting transcript regions
#'
#' Function that assigns each binding site in the \code{\link{BSFDataSet}} to its
#' hosting transcript region given an annotation database (\code{anno.annoDB}), or
#' a GRanges list / \code{\link{CompressedGRangesList}} (\code{anno.transcriptRegionList})
#' that holds the ranges for the transcript regions of interest.
#'
#' Since the assignment is based on the overlaps of annotated transcript ranges
#' and binding sites, no matching meta data is needed.
#'
#' In the case of transcript regions overlaps are very frequent. To resolve these
#' cases the \code{\link{overlaps}} argument can be used. Option `frequency`
#' will take the most frequently observed transcript region, option `hierarchy`
#' works in conjunction with a user defined rule (\code{overlaps.rule}).
#' Options `flag` and `remove` will label binding sites with an ambiguous tag or
#' remove all overlapping cases, respectively.
#'
#' @param object a \code{\link{BSFDataSet}} object with stored binding sites. This
#' means that ranges should be > 1
#' @param overlaps character; how overlapping transcript regions should be handled.
#' @param overlaps.rule character vector; a vector of transcript region
#' names that should be used to handle overlapping cases in a hierarchical manor.
#' The order of the vector is the order of the hierarchy.
#' @param anno.annoDB an object of class \code{OrganismDbi} that contains
#' the transcript region annotation.
#' @param anno.transcriptRegionList an object of class \code{\link{CompressedGRangesList}}
#' that holds an ranges for each transcript region
#'
#' @param normalize.exclude.upper numeric; percentage value that indicates the
#' upper boundary for transcript region length to be considered when calculating
#' normalization factors for regions.
#' @param normalize.exclude.lower numeric; percentage value that indicates the
#' lower boundary for transcript region length to be considered when calculating
#' normalization factors for regions.
#'
#' @param quiet logical; whether to print messages
#'
#' @return an object of class \code{\link{BSFDataSet}} with binding sites having
#' hosting transcript region information added to their meta columns.
#'
#' @importFrom dplyr rename_with desc
#' @importFrom GenomicFeatures cds intronsByTranscript threeUTRsByTranscript fiveUTRsByTranscript
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[2])
#'
#' bds = makeBindingSites(object = bds, bsSize = 9)
#' bds = assignToGenes(bds, anno.genes = gns)
#' bds = assignToTranscriptRegions(object = bds, anno.transcriptRegionList = regions)
#'
#' @export
assignToTranscriptRegions <- function(object, # bindingSiteFinder
                                      overlaps = c("frequency", "hierarchy", "flag", "remove"), # overlapping loci solution options
                                      overlaps.rule = NULL, # rule to apply when option hierarchy is selected
                                      anno.annoDB = NULL, # annotation data base of class OrganismDbi; can be NULL -> but then anno.transcriptRegionList must be provided
                                      anno.transcriptRegionList = NULL, # GRangesList object with the transcript features to be used; can be NULL -> but then anno.annoDB must be provided
                                      normalize.exclude.upper = 0.02,
                                      normalize.exclude.lower = 0.02,
                                      quiet = FALSE
) {
    # local variables
    datasource <- `.` <- Freq <- transcript_type <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is.logical(quiet))

    # Check if none is specified
    if (is.null(anno.annoDB) & is.null(anno.transcriptRegionList)) {
        msg = paste0("None of the required annotation sources anno.annoDB or anno.transcriptRegionList was specified. ")
        stop(msg)
    }
    # Check if both are specified
    if (!is.null(anno.annoDB) & !is.null(anno.transcriptRegionList)) {
        msg = paste0("Both of the required annotation sources anno.annoDB or anno.transcriptRegionList are specified. Please provide only one of the two. ")
        stop(msg)
    }
    # anno.annoDB should be used
    if (!is.null(anno.annoDB) & is.null(anno.transcriptRegionList)) {
        stopifnot(is(anno.annoDB, "OrganismDb"))
        if (!is.null(anno.transcriptRegionList)) {
            msg = paste0("Parameter anno.annoDB and anno.transcriptRegionList are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            datasource = "anno.annoDB"
            # extract relevant annotation - transcripts
            cdseq = cds(anno.annoDB)
            intrns = unlist(intronsByTranscript(anno.annoDB))
            utrs3 = unlist(threeUTRsByTranscript(anno.annoDB))
            utrs5 = unlist(fiveUTRsByTranscript(anno.annoDB))
            anno.transcriptRegionList = GRangesList(CDS = cdseq, Intron = intrns, UTR3 = utrs3, UTR5 = utrs5)
            names(anno.transcriptRegionList) = toupper(names(anno.transcriptRegionList))
        }
    }
    # anno.transcriptRegionList should be used
    if (is.null(anno.annoDB) & !is.null(anno.transcriptRegionList) ) {
        stopifnot(is(anno.transcriptRegionList, "GenomicRangesList"))
        if (!is.null(anno.annoDB)) {
            msg = paste0("Parameter anno.annoDB and anno.transcriptRegionList are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            datasource = "anno.transcriptRegionList"
            # extract relevant annotation - transcripts
            anno.transcriptRegionList = anno.transcriptRegionList
            names(anno.transcriptRegionList) = toupper(names(anno.transcriptRegionList))
        }
    }
    # handle options (hierarchy, frequency, remove, keep)
    overlaps = match.arg(overlaps, choices = c("frequency", "hierarchy", "flag", "remove"))
    # Check multiple loci options
    if (overlaps == "hierarchy" & is.null(overlaps.rule)) {
        msg1 = paste0("Multiple region problem is set to be handled by 'hierarchy', but no overlaps.rule is provided. \n")
        msg2 = paste0("Change the 'overlaps' option or provide a valid overlaps.rule. \n")
        stop(c(msg1, msg2))
    }
    overlaps.rule = toupper(overlaps.rule)
    if (! all(overlaps.rule %in% names(anno.transcriptRegionList))){
        msg = paste0("Regions defined in overlaps.rule (", paste(overlaps.rule, collapse = " > "),
                     ") does not match the input region names (", paste(names(anno.transcriptRegionList), collapse = ", "), ")")
        stop(msg)
    }
    if (overlaps == "hierarchy" & is.null(overlaps.rule)) {
        msg1 = paste0("Overlap are set to be solved by 'hierarchy', but no overlaps.rule is provided. \n")
        msg2 = paste0("Change the 'overlaps' option or provide a valid overlaps.rule. \n")
        stop(c(msg1, msg2))
    }

    # ---
    # Store function parameters in list
    optstr = list(source = datasource, overlaps = overlaps, overlaps.rule = overlaps.rule)
    object@params$assignToTranscriptRegions = optstr

    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    # Match transcript annotation with binding sites
    # --------------------------------------------------------------------------
    rngInitial = getRanges(object)
    rng = getRanges(object)

    currBsSize = object@params$bsSize

    cRange = rng - unique(floor(currBsSize / 2))
    countDf = lapply(anno.transcriptRegionList, function(x) {
        countOverlaps(cRange, x)
    }) %>% as.data.frame() %>% rename_with(toupper)
    # Add intergenic cases
    countDf = countDf %>% mutate(OTHER = ifelse(rowSums(.) == 0, 1, 0))

    # ---
    # Store results for plotting
    plotDf = countDf
    plotDf[plotDf >= 1] = 1
    object@plotData$assignToTranscriptRegions$dataOverlaps = plotDf

    # Deal with multiple loci
    # --------------------------------------------------------------------------
    # Count out how many cases there are
    totalBs = nrow(countDf)
    numberOfLociHit = rowSums(countDf != 0)
    ambigousBs = nrow(countDf[numberOfLociHit != 1,])
    ambigousBsFraction = paste0(round(ambigousBs / totalBs * 100, digits = 2), "%")

    # Case where there are no overlaps
    if (ambigousBs == 0) {
        msg0 = paste0("No binding sites (", format(totalBs, big.mark = ',', decimal.mark = "."), ") on overlapping transcript loci. \n")
        msg1 = paste0("Parameter 'overlaps' set to '", overlaps, "' not in use. \n")
        if(!quiet) message(c(msg0, msg1))
        # assign regions
        countDfSorted = countDf
        countDfSorted[countDfSorted > 0] = 1
        col_names = colnames(countDfSorted)
        transcriptRegion = ifelse(rowSums(countDfSorted) == 1, col_names[max.col(countDfSorted, ties.method = "first")], "Ambiguous")
        mcols(rng) = cbind.data.frame(mcols(rng), transcriptRegion = transcriptRegion)
    }
    # Cases where we have to resolve overlapping loci
    if (ambigousBs > 0) {
        msg0 = paste0(ambigousBsFraction, " (", ambigousBs,"/", totalBs, ")",
                      " of binding sites overlap with multiple different transcript regions in the given annotation. \n")
        if (overlaps == "hierarchy") {
            # apply hierarchy solution
            msg1 = paste0("Apply 'hierarchy' solution. \n ")
            if(!quiet) message(c(msg0, msg1))
            countDfSorted = countDf[, overlaps.rule]
            countDfSorted[countDfSorted > 0] = 1
            transcriptRegion = names(countDfSorted)[max.col(countDfSorted, ties.method = "first")]
            mcols(rng) = cbind.data.frame(mcols(rng), transcriptRegion = transcriptRegion)
        }
        if (overlaps == "frequency") {
            msg1 = paste0("Apply 'frequency' solution. \n ")
            if(!quiet) message(c(msg0, msg1))
            countDfSorted = countDf
            # countDfSorted[countDfSorted > 0] = 1
            # get overlaps.rule by frequency
            overlaps.rule = names(countDfSorted)[max.col(countDfSorted, ties.method = "first")] %>%
                table() %>% as.data.frame() %>% rename(transcript_type = 1) %>%
                arrange(desc(Freq)) %>% pull(transcript_type) %>% as.character()
            # apply frequency solution
            countDfSorted = countDfSorted[, overlaps.rule]
            transcriptRegion = names(countDfSorted)[max.col(countDfSorted, ties.method = "first")]
            mcols(rng) = cbind.data.frame(mcols(rng), transcriptRegion = transcriptRegion)
        }
        if (overlaps == "flag") {
            msg1 = paste0("Binding sites with multiple different overlapping transcript regions are being flaged with the tag 'Ambiguous'. \n")
            if(!quiet) message(c(msg0, msg1))
            countDfSorted = countDf
            countDfSorted[countDfSorted > 0] = 1
            # flag multiple loci binding sites
            col_names = colnames(countDfSorted)
            transcriptRegion = ifelse(rowSums(countDfSorted) == 1, col_names[max.col(countDfSorted, ties.method = "first")], "Ambiguous")
            mcols(rng) = cbind.data.frame(mcols(rng), transcriptRegion = transcriptRegion)
            # remove mulitple loci binding sitges
            mcols(rng) = cbind.data.frame(mcols(rng), transcriptRegion = transcriptRegion)
        }
        if (overlaps == "remove") {
            msg1 = paste0("Binding sites with multiple different overlapping transcript regions are being removed. This is not recommended! \n")
            if(!quiet) message(c(msg0, msg1))
            countDfSorted = countDf
            countDfSorted[countDfSorted > 0] = 1
            # flag multiple loci binding sites
            col_names = colnames(countDfSorted)
            transcriptRegion = ifelse(rowSums(countDfSorted) == 1, col_names[max.col(countDfSorted, ties.method = "first")], "Ambiguous")
            mcols(rng) = cbind.data.frame(mcols(rng), transcriptRegion = transcriptRegion)
            # remove mulitple loci binding sitges
            mcols(rng) = cbind.data.frame(mcols(rng), transcriptRegion = transcriptRegion)
            rng = rng[rng$transcriptRegion != "Ambiguous"]
        }
    }

    # ---
    # Calculate width for length based normalization
    ## normalize by all expressed regions
    debugList = list(object = object, rng = rng)
    # saveRDS(debugList, file = "./debugList.rds")
    norm.factors = .calcNormalizeFactors(object, rng, anno.transcriptRegionList,
                                         normalize.exclude.lower = normalize.exclude.lower,
                                         normalize.exclude.upper = normalize.exclude.upper)

    # ---
    # Store results for plotting
    dfPlot = rng$transcriptRegion %>%
        table() %>%
        as.data.frame() %>%
        rename("TranscriptRegion" = ".") %>%
        arrange(desc(Freq)) %>%
        left_join(., y = norm.factors, by = "TranscriptRegion")
    object@plotData$assignToTranscriptRegions$dataSpectrum = dfPlot

    # ---
    # Store for results
    resultLine = data.frame(
        funName = "assignToTranscriptRegions()", class = "binding sites",
        nIn = length(rngInitial), nOut = length(rng),
        per = paste0(round(length(rng)/ length(rngInitial), digits = 2)*100,"%"),
        options = paste0("Overlaps=", optstr$overlaps, ", Source=", datasource,
                         ifelse(!rlang::is_empty(optstr$overlaps.rule), paste0(", overlaps.rule=", paste(optstr$overlaps.rule, collapse = ">")), ""))
    )
    object@results = rbind(object@results, resultLine)

    outObject = setRanges(object, rng, quiet = quiet)
    return(outObject)
}


#' Function to estimate the appropriate binding site width together with the
#' optimal gene-wise filter level.
#'
#' This function tests different width of binding sites for different gene-wise
#' filtering steps. For each test the signal-to-score ratio is calculated. The
#' mean over all gene-wise filterings at each binding site width is used to
#' extract the optimal width, which serves as anchor to select the optimal
#' gene-wise filter.
#'
#' Parameter estimation is done on a subset of all crosslink sites
#' (\code{est.subsetChromosome}).
#'
#' Gene-level filter can be tested with varying
#' levels of accuracy ranging from `finest` to `coarse`, representing 1% and
#' 20% steps, respectively.
#'
#' Binding site computation at each step can be done on three different accuracy
#' level (\code{bsResolution}). Option `fine` is equal to a normal run
#' of the \code{\link{makeBindingSites}} function. `medium` will perform
#' a shorter version of the binding site computation, skipping some of the
#' refinement steps. Option `coarse` will approximate binding sites by merged
#' crosslinks regions, aligning the center at the site with the highest score.
#'
#' For each binding site in each set given the defined resolutions a signal-to-
#' flank score ratio is calculated and the mean of this score per set is returned.
#' Next a mean of means is created which results in a single score for each
#' binding site width that was tested. The width that yielded the highest score
#' is selected as optimal. In addtion the \code{minimumStepGain} option
#' allows control over the minimum additional gain in the score that a tested
#' width has to have to be selected as the best option.
#'
#' To enhance the sensitivity of the binding site estimation, the sensitivity
#' mode exists. In this mode crosslink sites undergo a pre-filtering and merging
#' step, to exclude potential artifical peaks (experimental-, mapping-biases).
#' If sensitivity mode is activated the \code{est.minWidth} option should be set
#' to 1.
#'
#' The optimal geneFilter is selected as the first one that passes the merged
#' mean of the selected optimal binding site width.
#'
#'
#' @param object a \code{\link{BSFDataSet}} object with stored crosslink sites.
#' This means that ranges should have a width = 1.
#' @param bsResolution character; level of resolution at which different binding
#' site width should be tested
#' @param geneResolution character; level of resolution at which gene-wise filtering
#' steps should be tested
#' @param est.maxBsWidth numeric; the largest binding site width which should
#' considered in the testing
#' @param est.minimumStepGain numeric; the minimum additional gain in the score
#' in percent the next binding site width has to have, to be selected as best option
#' @param est.maxSites numeric; maximum number of PureCLIP sites that are used;
#' @param est.subsetChromosome character; define on which chromosome the
#' estimation should be done in function \code{\link{estimateBsWidth}}
#' @param est.minWidth the minimum size of regions that are subjected to the
#' iterative merging routine, after the initial region concatenation.
#' @param est.offset constant added to the flanking count in the signal-to-flank
#' ratio calculation to avoid division by Zero
#'
#' @param sensitive logical; whether to enable sensitive pre-filtering before
#' binding site merging or not
#' @param sensitive.size numeric; the size (in nucleotides) of the merged
#' sensitive region
#' @param sensitive.minWidth numeric; the minimum size (in nucleoties) of the
#' merged sensitive region
#'
#' @param anno.annoDB an object of class \code{OrganismDbi} that contains
#' the gene annotation.
#' @param anno.genes an object of class \code{\link{GenomicRanges}} that represents
#' the gene ranges directly
#' @param bsResolution.steps numeric vector; option to use a user defined threshold
#' for binding site width directly. Overwrites \code{bsResolution}
#' @param geneResolution.steps numeric vector; option to use a user defined threshold
#' vector for gene-wise filtering resolution. Overwrites \code{geneResolution}
#' @param quiet logical; whether to print messages
#' @param veryQuiet logical; whether to suppress all messages
#' @param reportScoresPerBindingSite report the ratio score for each binding site
#' separately. Warning! This is for debugging and testing only. Downstream
#' functions can be impaired.
#' @param ... additional arguments passed to \code{\link{pureClipGeneWiseFilter}}
#'
#' @return an object of class \code{\link{BSFDataSet}} with binding sites with
#' the `params` slots `bsSize` and `geneFilter` being filled
#'
#' @importFrom dplyr summarise slice_head slice_tail
#' @importFrom utils head
#' @importFrom GenomicFeatures genes
#' @importFrom stats sd
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[2])
#' estimateBsWidth(bds, anno.genes = gns, est.maxBsWidth = 19,
#'  geneResolution = "coarse", bsResolution = "coarse", est.subsetChromosome = "chr22")
#'
#' @export
estimateBsWidth <- function(object, # BindingSiteFinder object
                            bsResolution = c("medium", "fine", "coarse"),
                            geneResolution = c("medium", "coarse", "fine", "finest"),
                            est.maxBsWidth = 13,
                            est.minimumStepGain = 0.02,
                            est.maxSites = Inf,
                            est.subsetChromosome = "chr1", #
                            est.minWidth = 2,
                            est.offset = 1,
                            sensitive = FALSE,
                            sensitive.size = 5,
                            sensitive.minWidth = 2,
                            anno.annoDB = NULL,
                            anno.genes = NULL,
                            bsResolution.steps = NULL,
                            geneResolution.steps = NULL,
                            quiet = TRUE,
                            veryQuiet = FALSE,
                            reportScoresPerBindingSite = FALSE,
                            ...
) {

    # initialize local variables
    bsSize <- ms <- growth_per <- sel <- geneWiseFilter <- est.option <- signalToFlankRatio <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is.logical(quiet))

    if (c(est.maxBsWidth %% 2) == 0) {
        stop("est.maxBsWidth is even. An odd number is required to have a distinct binding site center.")
    }

    # Check if none is specified
    if (is.null(anno.annoDB) & is.null(anno.genes)) {
        msg = paste0("None of the required annotation sources anno.annoDB or anno.genes was specified. ")
        stop(msg)
    }
    # Check if both are specified
    if (!is.null(anno.annoDB) & !is.null(anno.genes)) {
        msg = paste0("Both of the required annotation sources anno.annoDB or anno.genes are specified. Please provide only one of the two. ")
        stop(msg)
    }
    # Checks if anno.annoDB should be used
    if (!is.null(anno.annoDB) & is.null(anno.genes)) {
        stopifnot(is(anno.annoDB, "OrganismDb"))
        if (!is.null(anno.genes)) {
            msg = paste0("Parameter anno.annoDB and anno.genes are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            # extract relevant annotation
            anno.genes = genes(anno.annoDB)
        }
    }
    # Checks if anno.genes should be used
    if (is.null(anno.annoDB) & !is.null(anno.genes)) {
        stopifnot(is(anno.genes, "GenomicRanges"))
        if (!is.null(anno.annoDB)) {
            msg = paste0("Parameter anno.annoDB and anno.genes are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            # extract relevant annotation
            anno.genes = anno.genes
        }
    }

    # handle gene filter resolution
    geneResolution = match.arg(geneResolution, choices = c("medium", "coarse", "fine", "finest"))
    if (!is.null(geneResolution.steps)) {
        geneResolution.steps = geneResolution.steps
    } else {
        if (geneResolution == "coarse") {
            geneResolution.steps = seq(from = 0, to = 1, by = 0.2)
            geneResolution.steps = geneResolution.steps[1:length(geneResolution.steps)-1]
        }
        if (geneResolution == "medium") {
            geneResolution.steps = seq(from = 0, to = 1, by = 0.1)
            geneResolution.steps = geneResolution.steps[1:length(geneResolution.steps)-1]
        }
        if (geneResolution == "fine") {
            geneResolution.steps = seq(from = 0, to = 1, by = 0.05)
            geneResolution.steps = geneResolution.steps[1:length(geneResolution.steps)-1]
        }
        if (geneResolution == "finest") {
            geneResolution.steps = seq(from = 0, to = 1, by = 0.01)
            geneResolution.steps = geneResolution.steps[1:length(geneResolution.steps)-1]
        }
    }

    # handle binding site compute resolution
    bsResolution = match.arg(bsResolution, choices = c("medium", "fine", "coarse"))

    # handle est.maxBsWidth
    if (!is.null(bsResolution.steps)) {
        bsResolution.steps = bsResolution.steps
    } else {
        bsResolution.steps = seq(from = 3, to = est.maxBsWidth, by = 2)
    }

    # ---
    # Store function parameters in list
    optstr = list(bsResolution = bsResolution,
                  geneResolution = geneResolution,
                  est.maxBsWidth = est.maxBsWidth,
                  est.maxSites = est.maxSites,
                  est.minimumStepGain = est.minimumStepGain,
                  est.subsetChromosome = est.subsetChromosome,
                  est.minWidth = est.minWidth,
                  est.offset = est.offset,
                  sensitive = sensitive,
                  sensitive.size = sensitive.size,
                  sensitive.minWidth = sensitive.minWidth)
    object@params$estimateBsWidth = optstr

    # PREPARE TEST RANGES + SIGNAL
    # --------------------------------------------------------------------------
    checkRng = getRanges(object)
    # check if subset is part of the seqnames from ranges
    if (!all(est.subsetChromosome %in% levels(seqnames(checkRng)))){
        msg = paste0("Chromosome to estimate on (", est.subsetChromosome, "), is not included in the ranges: ", paste(levels(seqnames(checkRng)), collapse = ", "))
        stop(msg)
    }

    # limit estimation to a specific chromosome
    if (!is.null(est.subsetChromosome)) {
        redObj = .subsetByChr(object, chr = est.subsetChromosome, quiet = quiet)
    } else {
        redObj = object
    }

    # limit estimation to a maximum number of sites
    estRng = getRanges(redObj)
    if (length(estRng) > est.maxSites) {
        estRng = head(estRng, est.maxSites)
        redObj = setRanges(redObj, estRng, quiet = quiet)
    }

    # limit the clip signal to the ranges used for estimation (plus some extra offset)
    if (!is.null(est.subsetChromosome)) {
       # reduce signal to frame
        maxFrame = ceiling((max(bsResolution.steps) *3))
        redObj = .reduceSignalToFrame(redObj, frame = maxFrame, quiet = quiet)
    } else {
        # don't reduce signal
        redObj = redObj
    }

    # collapse signal from replicates
    sgnMerge = .collapseSamples(getSignal(redObj))


    # MAIN COMPUTE
    # --------------------------------------------------------------------------

    # a counter that shows to how many percent the computation is done
    counterTotalIterations = length(geneResolution.steps) * length(bsResolution.steps)
    counterCurrentIterations = 0

    # print start message
    msg = paste0("Estimation at: ", round((counterCurrentIterations/counterTotalIterations)*100 ), "%"  )
    if(!veryQuiet) print(msg)

    # calculate binding sites for each filter step and width
    scoreAllDf = lapply(geneResolution.steps, function(bsFilterStep){
        # bsFilterStep = 0.2

        # apply current gene-wise filter (has two modes)
        if (isTRUE(sensitive)) {
            # -> sensitive mode
            # define regions in which sites should be kept
            cRng = getRanges(redObj)
            cRngMerge = reduce(cRng, min.gapwidth = sensitive.size)
            cRngMerge = cRngMerge[width(cRngMerge) >= sensitive.minWidth]

            # currFilterObj = pureClipGeneWiseFilter(object = redObj, anno.genes = anno.genes, cutoff = bsFilterStep, quiet = quiet, ...)
            currFilterObj = pureClipGeneWiseFilter(object = redObj, anno.genes = anno.genes, cutoff = bsFilterStep, quiet = quiet)
            cRngFilter = getRanges(currFilterObj)
            currRng = IRanges::subsetByOverlaps(cRngFilter, cRngMerge)

        } else {
            # -> normal mode
            # currFilterObj = pureClipGeneWiseFilter(object = redObj, anno.genes = anno.genes, cutoff = bsFilterStep, quiet = quiet, ...)
            currFilterObj = pureClipGeneWiseFilter(object = redObj, anno.genes = anno.genes, cutoff = bsFilterStep, quiet = quiet)
            currRng = getRanges(currFilterObj)
        }

        # calculate binding sites for current bsWidth
        rngPerWidth = lapply(bsResolution.steps, function(bsWidthStep){
            if (bsResolution == "fine") {
                # calculate with full binding sites
                currBsObj = makeBindingSites(object = currFilterObj,
                                             bsSize = bsWidthStep,
                                             minWidth = est.minWidth,
                                             quiet = quiet)
                currBs = getRanges(currBsObj)
            }
            if (bsResolution == "medium") {
                # approximate binding sites by a single merge and extend round
                currBs = .approximateBindingSites_medium(rng = currRng,
                                                         sgn = sgnMerge,
                                                         bsSize = bsWidthStep,
                                                         minWidth = est.minWidth)
            }
            if (bsResolution == "coarse") {
                # approximate binding sites by reduced pureclip sites center
                currBs = .approximateBindingSites_coarse(rng = currRng,
                                                         bsSize = bsWidthStep,
                                                         minWidth = est.minWidth)
            }
            return(currBs)
        })
        rngPerWidth = unlist(GRangesList(rngPerWidth))

        # handle all plus ranges
        cRangePlus = subset(rngPerWidth, strand == "+")
        if (length(cRangePlus) > 0) {
            bsSumPlus = sum(sgnMerge$signalPlus[cRangePlus])
            extendedRangePlus = cRangePlus + cRangePlus$bsSize
            exSumPlus = sum(sgnMerge$signalPlus[extendedRangePlus])
            # extendedRangePlus$signalToFlankRatio = (bsSumPlus) / (((exSumPlus - bsSumPlus) / 2) + est.offset)
            extendedRangePlus$signalToFlankRatio = (bsSumPlus / exSumPlus)
        } else {
            extendedRangePlus = cRangePlus
        }

        # handle all minus ranges
        cRangeMinus = subset(rngPerWidth, strand == "-")
        if (length(cRangeMinus) > 0) {
            bsSumMinus = sum(sgnMerge$signalMinus[cRangeMinus])
            extendedRangeMinus = cRangeMinus + cRangeMinus$bsSize
            exSumMinus = sum(sgnMerge$signalMinus[extendedRangeMinus])
            # extendedRangeMinus$signalToFlankRatio = (bsSumMinus) / (((exSumMinus - bsSumMinus) / 2) + est.offset)
            extendedRangeMinus$signalToFlankRatio = (bsSumMinus / exSumMinus)
        } else {
            extendedRangeMinus = cRangeMinus
        }

        if (isTRUE(reportScoresPerBindingSite)) {
            # report the scores for each binding site
            # -> does not allow the usage of the standard plot
            df = rbind(as.data.frame(mcols(extendedRangeMinus)), as.data.frame(mcols(extendedRangePlus))) %>%
                group_by(bsSize) %>%
                mutate(geneWiseFilter = bsFilterStep)

        } else {
            # combine both and make results dataframe
            df = rbind(as.data.frame(mcols(extendedRangeMinus)), as.data.frame(mcols(extendedRangePlus))) %>%
                group_by(bsSize) %>%
                summarize(signalToFlankRatio = median(signalToFlankRatio), .groups = "keep") %>%
                mutate(geneWiseFilter = bsFilterStep)
        }

        # update chunk counter
        counterCurrentIterations <<- counterCurrentIterations + length(bsResolution.steps)
        msg = paste0("Estimation at: ", round((counterCurrentIterations/counterTotalIterations)*100 ), "%"  )
        if(!veryQuiet) print(msg)
        return(df)
    })

    df = do.call("rbind", scoreAllDf)

    # ---
    # Store results for plotting
    if (isTRUE(reportScoresPerBindingSite)){
        object@plotData$estimateBsWidth$dataEnhanced = df
    } else {
        object@plotData$estimateBsWidth$data = df
    }

    resEstimate = .findFirstMaximum(df, est.minimumStepGain = est.minimumStepGain)
    est.bsSize = resEstimate$est.bsSize
    est.geneFilter = resEstimate$est.geneFilter
    est.option = resEstimate$est.option

    if (est.option == "local") {
        msg = paste0("No global maximum found with est.maxBsWidth=", est.maxBsWidth, ". Using local maximum instead.\n")
        if(!veryQuiet) warning(msg)
    }
    if (est.option == "fallback") {
        msg = paste0("No global or local maximum found with est.maxBsWidth=", est.maxBsWidth, ". Falling back to simple maximum instead.\n")
        if(!veryQuiet) warning(msg)
    }
    if (est.option == "error") {
        msg = paste0("No maximum found at all. Raise 'est.minimumStepGain' argument/, change bsResolution and/or geneResolution arguments/ raise est.maxBsWidth and check visually. \n")
        stop(msg)
    }

    # Store estimate type in options (global/ local)
    optstr$option = est.option
    object@params$estimateBsWidth = optstr

    # Store estimated parameters
    object@params$bsSize = est.bsSize
    object@params$geneFilter = est.geneFilter

    # ---
    # Store for results
    rng = getRanges(object)
    resultLine = data.frame(
        funName = "estimateBsWidth()", class = "estimate",
        nIn = length(rng), nOut = length(rng),
        per = paste0(round(length(rng)/ length(rng), digits = 2)*100,"%"),
        options = paste0("bsResolution=", bsResolution, ", geneResolution=", geneResolution,
                         ", est.maxBsWidth=", est.maxBsWidth, ", est.maxSites=", est.maxSites, ", est.subsetChromosome=", paste(est.subsetChromosome, collapse = ","))
    )
    object@results = rbind(object@results, resultLine)

    return(object)
}



