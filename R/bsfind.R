


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

#' @param est.geneResolution character; level of resolution of the gene-wise
#' filtering in function \code{\link{estimateBsWidth}}
#' @param est.bsResolution character; level of resolution of the binding site
#' width in function \code{\link{estimateBsWidth}}
#' @param est.subsetChromosome character; define on which chromosome the
#' estimation should be done in function \code{\link{estimateBsWidth}}
#'
#' @param cutoff.globalFilter numeric; defines the cutoff for which sites to
#' keep, the smallest step is 1\% (0.01) in function
#' \code{\link{pureClipGlobalFilter}}
#' @param cutoff.geneWiseFilter numeric; defines the cutoff for which sites to
#' remove in in function \code{\link{pureClipGeneWiseFilter}}. The smallest step
#' is 1\% (0.01). A cutoff of 5\% will remove the lowest 5\% sites, given their
#' score, on each gene, thus keeping the strongest 95\%.
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
#' @param ... additional arguments passed to \code{\link{makeBindingSites}} and
#' \code{\link{reproducibilityFilter}}
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
        est.geneResolution = "medium",
        est.bsResolution = "medium",
        est.subsetChromosome = "chr1",
        # cutoffs
        cutoff.globalFilter = 0.01,
        cutoff.geneWiseFilter = NULL,
        # overlap arguments
        overlaps.geneWiseFilter = "keepSingle",
        overlaps.geneAssignment = "frequency",
        overlaps.rule.geneAssignment = NULL,
        overlaps.TranscriptRegions = "frequency",
        overlaps.rule.TranscriptRegions = NULL ,
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
    if(!veryQuiet) message("estimateBsWidth...")
    if (is.null(bsSize) | is.null(cutoff.geneWiseFilter)) {
        obj = estimateBsWidth(obj,
                              geneResolution = est.geneResolution,
                              bsResolution = est.bsResolution,
                              est.subsetChromosome = est.subsetChromosome,
                              anno.annoDB = anno.annoDB,
                              anno.genes = anno.genes, ...)
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
    obj = makeBindingSites(obj, bsSize = bsSize, quiet = quiet, ...)

    if(!veryQuiet) message("reproducibilityFilter...")
    obj = reproducibilityFilter(obj, returnType = "BSFDataSet", quiet = quiet, ...)

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

    return(obj)
}
