


#' RBP binding site definition for iCLIP data
#'
#' This is the main function that performs the binding site definition analysis
#' through the following steps:
#' \enumerate{
#' \item Filter PureCLIP sites by their score distribution: \code{\link{pureClipGlobalFilter}}
#' \item Estimate the appropriate binding site width together with the optimal gene-wise filter level: \code{\link{estimateBsWidth}}
#' \item Filter PureCLIP sites by their score distribution per gene: \code{\link{pureClipGeneWiseFilter}}
#' \item Define equally sized binding sites: \code{\link{makeBindingSites}}
#' \item Perform replicate reproducibility filter: \code{\link{reproducibilityFilter}}
#' \item Assign binding sites to their hosting genes: \code{\link{assignToGenes}}
#' \item Assign binding sites to their hosting transcript regions: \code{\link{assignToTranscriptRegions}}
#' \item Re-assign PureCLIP scores to binding sites: \code{\link{annotateWithScore}}
#' }
#'
#' For complete details on each step, see the manual pages of the respective
#' functions. The \code{BSFind} function returns a \code{\link{BSFDataSet}}
#' with ranges merged into binding sites. A full flowchart for the entire process
#' can be visualized with \code{\link{processingStepsFlowChart}}. For each of the
#' individual steps dedicated diagnostic plots exists. Further information can be
#' found in our Bioconductor vignette:
#' \url{https://www.bioconductor.org/packages/release/bioc/html/BindingSiteFinder.html}
#'
#' If no binding site size is provided through \code{bsSize}, then
#' \code{\link{estimateBsWidth}} is called to estimate the optimal size for the
#' given data-set. The result of this estimation can be looked at with
#' \code{\link{estimateBsWidthPlot}} and arguments can be adjusted if needed.
#'
#'
#' @param object a \code{\link{BSFDataSet}} object with stored ranges
#' @param bsSize an odd integer value specifying the size of the output
#' binding sites
#'
#' @param est.bsResolution character; level of resolution of the binding site
#' width in function \code{\link{estimateBsWidth}}
#' @param est.geneResolution character; level of resolution of the gene-wise
#' filtering in function \code{\link{estimateBsWidth}}
#' @param est.maxBsWidth numeric; the largest binding site width which should
#' considered in the testing
#' @param est.minimumStepGain numeric; the minimum additional gain in the score
#' in percent the next binding site width has to have, to be selected as best option
#' @param est.maxSites numeric; maximum number of PureCLIP sites that are used
#' @param est.subsetChromosome character; define on which chromosome the
#' estimation should be done in function \code{\link{estimateBsWidth}}
#' @param est.minWidth the minimum size of regions that are subjected to the
#' iterative merging routine, after the initial region concatenation.
#' @param est.offset constant added to the flanking count in the signal-to-flank
#' ratio calculation to avoid division by Zero
#' @param est.sensitive logical; whether to enable sensitive pre-filtering before
#' binding site merging or not
#' @param est.sensitive.size numeric; the size (in nucleotides) of the merged
#' sensitive region
#' @param est.sensitive.minWidth numeric; the minimum size (in nucleoties) of the
#' merged sensitive region
#'
#' @param merge.minWidth the minimum size of regions that are subjected to the
#' iterative merging routine, after the initial region concatenation.
#' @param merge.minCrosslinks the minimal number of positions to overlap with at least
#' one crosslink event in the final binding sites
#' @param merge.minClSites the minimal number of crosslink sites that have to
#' overlap a final binding site
#' @param merge.CenterIsClSite logical, whether the center of a final binding
#' site must be covered by an initial crosslink site
#' @param merge.CenterIsSummit logical, whether the center of a final binding
#' site must exhibit the highest number of crosslink events
#'
#' @param cutoff.globalFilter numeric; defines the cutoff for which sites to
#' keep, the smallest step is 1\% (0.01) in function
#' \code{\link{pureClipGlobalFilter}}
#' @param cutoff.geneWiseFilter numeric; defines the cutoff for which sites to
#' remove in in function \code{\link{pureClipGeneWiseFilter}}. The smallest step
#' is 1\% (0.01). A cutoff of 5\% will remove the lowest 5\% sites, given their
#' score, on each gene, thus keeping the strongest 95\%.
#'
#' @param repro.cutoff numeric; percentage cutoff to be used for the
#' reproducibility quantile filtering
#' @param repro.nReps numeric; number of replicates that must meet the cutoff
#' defined in \code{repro.cutoff} for a binding site to be called reproducible.
#' Defaults to N-1.
#' @param repro.minCrosslinks numeric; minimal number of crosslinks a binding
#' site needs to have to be called reproducible. Acts as a lower boundary for
#' \code{repro.cutoff}. Defaults to 1.
#'
#' @param overlaps.geneWiseFilter character; how overlaps should be handled in
#'  \code{\link{pureClipGeneWiseFilter}}
#' @param overlaps.geneAssignment character; how overlaps should be handled in
#'  \code{\link{assignToGenes}}
#' @param overlaps.rule.geneAssignment character vector; a vector of gene types
#' that should be used to handle overlaps if option 'hierarchy' is selected
#' for \code{\link{assignToGenes}}. The order of the vector is the order of
#' the hierarchy.
#' @param overlaps.TranscriptRegions character; how overlaps should be handled in
#'  \code{\link{assignToTranscriptRegions}}
#' @param overlaps.rule.TranscriptRegions character vector; a vector of gene types
#' that should be used to handle overlaps if option 'hierarchy' is selected
#' for \code{\link{assignToTranscriptRegions}}. The order of the vector is the order of
#' the hierarchy.
#'
#' @param stf.flank character; how the flanking region shoule be set. Options are
#' 'bs', 'manual'
#' @param stf.flank.size numeric; if flank='manual' provide the desired flanking size
#'
#' @param match.score character; meta column name of the crosslink site
#' @param match.geneID character; meta column name of the genes
#' @param match.geneName character; meta column name of the gene name
#' @param match.geneType character; meta column name of the gene type
#' @param match.ranges.score a GRanges object, with numeric column for the score
#' to match in function \code{\link{annotateWithScore}}
#' @param match.option.score character; meta column name of the crosslink site
#' in function \code{\link{annotateWithScore}}
#'
#' @param anno.annoDB an object of class \code{OrganismDbi} that contains
#' the gene annotation.
#' @param anno.genes an object of class \code{\link{GenomicRanges}} that represents
#' the gene ranges directly
#' @param anno.transcriptRegionList an object of class \code{\link{CompressedGRangesList}}
#' that holds an ranges for each transcript region
#'
#' @param quiet logical; whether to print messages
#' @param veryQuiet logical; whether to suppress all messages
#' @param ... additional arguments passed to \code{\link{estimateBsWidth}},
#' \code{\link{makeBindingSites}} and \code{\link{reproducibilityFilter}}
#'
#'
#' @return an object of class \code{\link{BSFDataSet}} with ranges merged into
#' binding sites given the inputs.
#'
#' @seealso \code{\link{BSFDataSet}}, \code{\link{estimateBsWidth}},
#' \code{\link{pureClipGlobalFilter}}, \code{\link{pureClipGeneWiseFilter}},
#' \code{\link{assignToGenes}}, \code{\link{assignToTranscriptRegions}},
#' \code{\link{annotateWithScore}}, \code{\link{reproducibilityFilter}}
#'
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' # Load genes
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#' # load transcript regions
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[2])
#' BSFind(object = bds, bsSize = 9, anno.genes = gns,
#'  anno.transcriptRegionList = regions, est.subsetChromosome = "chr22")
#'
#' @export
BSFind <- function(
        # input object
        object,
        # binding site size
        bsSize = NULL,
        # binding site size estimation
        est.bsResolution = "medium",
        est.geneResolution = "medium",
        est.maxBsWidth = 13,
        est.minimumStepGain = 0.02,
        est.maxSites = Inf,
        est.subsetChromosome = "chr1",
        est.minWidth = 2,
        est.offset = 1,
        est.sensitive = FALSE,
        est.sensitive.size = 5,
        est.sensitive.minWidth = 2,
        # binding site merging
        merge.minWidth = 2,
        merge.minCrosslinks = 2,
        merge.minClSites = 1,
        merge.CenterIsClSite = TRUE,
        merge.CenterIsSummit = TRUE,
        # cutoffs
        cutoff.globalFilter = 0.01,
        cutoff.geneWiseFilter = NULL,
        # reproducibility
        repro.cutoff = NULL,
        repro.nReps = NULL,
        repro.minCrosslinks = 1,
        # overlap arguments
        overlaps.geneWiseFilter = "keepSingle",
        overlaps.geneAssignment = "frequency",
        overlaps.rule.geneAssignment = NULL,
        overlaps.TranscriptRegions = "frequency",
        overlaps.rule.TranscriptRegions = NULL ,
        # signal-to-flank
        stf.flank = "bs",
        stf.flank.size = NULL,
        # matching arguments
        match.score = "score",
        match.geneID = "gene_id",
        match.geneName = "gene_name",
        match.geneType = "gene_type",
        match.ranges.score = NULL,
        match.option.score = "max",
        # data sources
        anno.annoDB = NULL,
        anno.genes = NULL,
        anno.transcriptRegionList = NULL,
        # controls
        quiet = TRUE,
        veryQuiet = FALSE,
        ...
) {
    # INPUT CHECKS
    # --------------------------------------------------------------------------
    # General
    # ---
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is.logical(quiet))
    stopifnot(is.logical(veryQuiet))

    # check meta data
    this.meta = getMeta(object)
    if (length(levels(this.meta$condition)) > 1) {
        msg0 = paste0("Found ", length(levels(this.meta$condition)), " different conditions in the input object.\n")
        msg1 = paste0("BSFind can only be used on data from a single condition.\n")
        msg2 = paste0("Please run BSFind sparately for each condition, then combine both objects with combineBSF.\n ")
        stop(c(msg0,msg1,msg2))
    }

    # Check annotation source
    # ---
    if (is.null(anno.annoDB) & is.null(anno.genes) & is.null(anno.transcriptRegionList)) {
        msg = paste0("None of the required annotation sources anno.annoDB or anno.genes was specified. ")
        stop(msg)
    }
    if (!is.null(anno.annoDB) & !is.null(anno.genes) | !is.null(anno.annoDB) & !is.null(anno.transcriptRegionList)) {
        msg = paste0("Both of the required annotation sources anno.annoDB or anno.genes/ anno.transcriptRegionList are specified. Please provide only one source type. ")
        stop(msg)
    }
    if (!is.null(anno.annoDB) & is.null(anno.genes) & is.null(anno.transcriptRegionList)) {
        stopifnot(is(anno.annoDB, "OrganismDb"))
        if (!is.null(anno.genes)) {
            msg = paste0("Parameter anno.annoDB and anno.genes are specified at the same time. Use only one of them.")
            stop(msg)
        }
    }
    if (is.null(anno.annoDB) & !is.null(anno.genes) & !is.null(anno.transcriptRegionList)) {
        stopifnot(is(anno.genes, "GenomicRanges"))
        stopifnot(is(anno.transcriptRegionList, "GenomicRangesList"))
    }
    if (is.null(anno.annoDB) & !is.null(anno.genes) & is.null(anno.transcriptRegionList)) {
        msg = paste0("Parameter anno.transcriptRegionList is missing, while anno.annoDB not specified.")
        stop(msg)
    }
    if (is.null(anno.annoDB) & is.null(anno.genes) & !is.null(anno.transcriptRegionList)) {
        msg = paste0("Parameter anno.genes is missing, while anno.annoDB not specified.")
        stop(msg)
    }

    # Check values
    # ---
    obj_colname = colnames(mcols(getRanges(object)))
    if (! match.score %in% obj_colname) {
        msg = paste0("Column specified in `match.score` (", match.score, "), does not match any column in object score.")
        stop()
    }

    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    obj = object

    # Begin main workflow
    if(!veryQuiet) message("pureClipGlobalFilter...")
    obj = pureClipGlobalFilter(obj,
                               match.score = match.score,
                               cutoff = cutoff.globalFilter,
                               quiet = quiet)

    # Check if parameter bsSize and cutoff.geneWiseFilter should be estimated
    if (is.null(bsSize) | is.null(cutoff.geneWiseFilter)) {
        if(!veryQuiet) message("estimateBsWidth...")
        obj = estimateBsWidth(obj,
                              bsResolution = est.bsResolution,
                              geneResolution = est.geneResolution,
                              est.maxBsWidth = est.maxBsWidth,
                              est.minimumStepGain = est.minimumStepGain,
                              est.maxSites = est.maxSites,
                              est.subsetChromosome = est.subsetChromosome,
                              est.minWidth = est.minWidth,
                              est.offset = est.offset,
                              sensitive = est.sensitive,
                              sensitive.size = est.sensitive.size,
                              sensitive.minWidth = est.sensitive.minWidth,
                              anno.annoDB = anno.annoDB,
                              anno.genes = anno.genes,
                              quiet = quiet,
                              veryQuiet = veryQuiet)
        if (is.null(bsSize) & !is.null(cutoff.geneWiseFilter)) {
            obj@params$geneFilter = cutoff.geneWiseFilter
        }
        if (!is.null(bsSize) & is.null(cutoff.geneWiseFilter)) {
            obj@params$bsSize = bsSize
        }
    }

    # Continue main workflow
    if(!veryQuiet) message("pureClipGeneWiseFilter...")
    obj = pureClipGeneWiseFilter(obj,
                                 cutoff = cutoff.geneWiseFilter,
                                 overlaps = overlaps.geneWiseFilter,
                                 match.score = match.score,
                                 match.geneID = match.geneID,
                                 anno.annoDB = anno.annoDB,
                                 anno.genes = anno.genes,
                                 quiet = quiet)

    if(!veryQuiet) message("makeBindingSites...")
    obj = makeBindingSites(obj,
                           bsSize = bsSize,
                           minWidth = merge.minWidth,
                           minCrosslinks = merge.minCrosslinks,
                           minClSites = merge.minClSites,
                           centerIsClSite = merge.CenterIsClSite,
                           centerIsSummit = merge.CenterIsSummit,
                           quiet = quiet)

    if(!veryQuiet) message("reproducibilityFilter...")
    obj = reproducibilityFilter(obj,
                                returnType = "BSFDataSet",
                                cutoff = repro.cutoff,
                                nReps = repro.nReps,
                                minCrosslinks = repro.minCrosslinks,
                                quiet = quiet)

    if(!veryQuiet) message("assignToGenes...")
    obj = assignToGenes(obj,
                        overlaps = overlaps.geneAssignment,
                        overlaps.rule = overlaps.rule.geneAssignment,
                        anno.annoDB = anno.annoDB,
                        anno.genes = anno.genes,
                        match.geneID = match.geneID,
                        match.geneName = match.geneName,
                        match.geneType = match.geneType,
                        quiet = quiet)

    if(!veryQuiet) message("assignToTranscriptRegions...")
    obj = assignToTranscriptRegions(obj,
                                    overlaps = overlaps.TranscriptRegions,
                                    overlaps.rule = overlaps.rule.TranscriptRegions,
                                    anno.annoDB = anno.annoDB,
                                    anno.transcriptRegionList = anno.transcriptRegionList,
                                    quiet = quiet)

    if(!veryQuiet) message("annotateWithScore...")
    if (is.null(match.ranges.score)) {
        match.ranges.score = getRanges(object)
        match.score = match.score
    }
    obj = annotateWithScore(obj,
                            match.ranges = match.ranges.score,
                            match.score = match.score,
                            match.option = match.option.score,
                            quiet = quiet)

    if(!veryQuiet) message("calculateSignalToFlankScore...")
    obj = calculateSignalToFlankScore(obj,
                                      flank = stf.flank,
                                      flank.size = stf.flank.size,
                                      quiet = quiet)
    return(obj)
}
