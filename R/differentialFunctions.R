

#' Compute background coverage for binding sites per gene
#'
#' This function computes the background coverage used for the differential
#' binding analysis to correct for transcript level changes. Essentially,
#' the crosslink signal on each gene is split into crosslinks that can be
#' attributed to the binding sites and all other signal that can be attributed
#' to the background.
#'
#' To avoid that crosslinks from binding sites contaminate the background counts
#' a protective region around each binding sites can be spanned with
#' \code{use.offset} the default width of the offset region is half of the
#' binding site width, but can also be changed with the \code{ranges.offset}
#' parameter.
#'
#' Additional region that one wants to exclude from contributing to the
#' background signal can be incorporated as \code{GRanges} objects through
#' the \code{blacklist} option.
#'
#' It is expected that binding sites are assigned to hosting genes prior to
#' running this funciton (see \code{\link{BSFind}}). This means a unique gene ID
#' is present in the meta columns of each binding site ranges. If this is not the
#' case one can invoce the binding site to gene assignment with
#' \code{generate.geneID.bs}. The same is true for the blacklist regions with
#' option \code{generate.geneID.blacklist}.
#'
#' It is expected that all binding sites are of the same size
#' (See \code{\link{BSFind}} on how to achieve this). If this is however not
#' the case and one wants to keep binding sites of different with then option
#' \code{force.unequalSites} can be used.
#'
#' @param object a \code{\link{BSFDataSet}} object with two conditions
#' @param anno.annoDB an object of class \code{OrganismDbi} that contains
#' the gene annotation.
#' @param anno.genes an object of class \code{\link{GenomicRanges}} that represents
#' the gene ranges directly
#' @param blacklist GRanges; genomic ranges where the signal should be
#' excluded from the background
#' @param use.offset logical; if an offset region around the binding sites should
#' be used on which the signal is excluded from the background
#' @param ranges.offset numeric; number of nucleotides the offset window around
#' each binding site should be wide (defaults to 1/2 binding site width - NULL)
#' @param match.geneID.gene character; the name of the column with the gene ID
#' in the genes meta columns used for matching binding sites to genes
#' @param match.geneID.bs character; the name of the column with the gene ID
#' in the binding sites meta columns used for matching binding sites to genes
#' @param match.geneID.blacklist character; the name of the column with the gene
#' ID in the blacklist meta columns used for matching the blacklist regions
#' with the genes
#' @param generate.geneID.bs logical; if the binding site to gene matching
#' should be performed if no matching gene ID is provided
#' @param generate.geneID.blacklist logical; if the blacklist to gene matching
#' should be performed if no matching gene ID is provided
#' @param uniqueID.gene character; column name of a unique ID for all genes
#' @param uniqueID.bs character; column name of a unique ID for all binding sites
#' @param uniqueID.blacklist character; column name of a unique ID for all
#' blacklist regions
#' @param force.unequalSites logical; if binding sites of not identical width
#' should be allowed or not
#' @param quiet logical; whether to print messages or not
#' @param veryQuiet logical; whether to print messages or not
#' @param ... additional arguments; passed to \code{\link{assignToGenes}}
#'
#' @return an object of class \code{\link{BSFDataSet}} with counts for binding
#' sites, background and total gene added to the meta column of the ranges
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#'
#' # make binding sites
#' bds = makeBindingSites(bds, bsSize = 7)
#' bds = assignToGenes(bds, anno.genes = gns)
#'
#' # change meta data
#' m = getMeta(bds)
#' m$condition = factor(c("WT", "WT", "KO", "KO"), levels = c("WT", "KO"))
#' bds = setMeta(bds, m)
#'
#' # change signal
#' s = getSignal(bds)
#' names(s$signalPlus) = paste0(m$id, "_", m$condition)
#' names(s$signalMinus) = paste0(m$id, "_", m$condition)
#' bds = setSignal(bds, s)
#'
#' # make example blacklist region
#' myBlacklist = getRanges(bds)
#'set.seed(1234)
#' myBlacklist = sample(myBlacklist, size = 500) + 4
#'
#' # make background
#' bds.b1 = calculateBsBackground(bds, anno.genes = gns)
#'
#' # make background - no offset
#' bds.b2 = calculateBsBackground(bds, anno.genes = gns, use.offset = FALSE)
#'
#' # make background - use blacklist
#' bds.b3 = calculateBsBackground(bds, anno.genes = gns, blacklist = myBlacklist)
#'
#' @export
calculateBsBackground <- function(object,
                                  anno.annoDB = NULL,
                                  anno.genes = NULL,
                                  blacklist = NULL,
                                  use.offset = TRUE,
                                  ranges.offset = NULL,
                                  match.geneID.gene = "gene_id",
                                  match.geneID.bs = "geneID",
                                  match.geneID.blacklist = "geneID",
                                  generate.geneID.bs = FALSE,
                                  generate.geneID.blacklist = FALSE,
                                  uniqueID.gene = "gene_id",
                                  uniqueID.bs = "bsID",
                                  uniqueID.blacklist = "bsID",
                                  force.unequalSites = FALSE,
                                  quiet = FALSE,
                                  veryQuiet = TRUE,
                                  ...
){
    # bind local variables
    geneID <- . <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is.logical(quiet))
    stopifnot(is.logical(veryQuiet))
    stopifnot(is.logical(use.offset))
    stopifnot(is.logical(generate.geneID.bs))
    stopifnot(is.logical(generate.geneID.blacklist))
    stopifnot(is.logical(force.unequalSites))

    # check meta data
    this.meta = getMeta(object)
    if (length(levels(this.meta$condition)) == 1) {
        msg0 = paste0("Found only one conditions in the input object.\n")
        msg1 = paste0("Make sure to combine objects from different conditions in order to compare them.\n")
        stop(c(msg0,msg1))
        # TODO this could also be just a warning
    }

    # Check gene annotation source
    # -> Check if none is specified
    if (is.null(anno.annoDB) & is.null(anno.genes)) {
        msg = paste0("None of the required annotation sources anno.annoDB or anno.genes was specified. ")
        stop(msg)
    }
    # -> Check if both are specified
    if (!is.null(anno.annoDB) & !is.null(anno.genes)) {
        msg = paste0("Both of the required annotation sources anno.annoDB or anno.genes are specified. Please provide only one of the two. ")
        stop(msg)
    }
    # -> Checks if anno.annoDB should be used
    if (!is.null(anno.annoDB) & is.null(anno.genes)) {
        stopifnot(is(anno.annoDB, "OrganismDb"))
        datasource = "anno.annoDB"
        # extract relevant annotation
        anno.genes = genes(anno.annoDB, columns = c("GENEID"))
        # Create matching vectors for columns from input annotation
        selectID = as.character(anno.genes$GENEID)
    }
    # -> Checks if anno.genes should be used
    if (is.null(anno.annoDB) & !is.null(anno.genes)) {
        stopifnot(is(anno.genes, "GenomicRanges"))

        datasource = "anno.genes"
        # extract relevant annotation
        anno.genes = anno.genes
        # check correct annotation columns
        annoColNames = colnames(mcols(anno.genes))
        if (!match.geneID.gene %in% annoColNames) {
            msg0 = paste0("The column '", match.geneID.gene, "', for match.geneID.gene is not present.\n")
            msg1 = paste0("Provide a column for gene id matching. \n")
            stop(c(msg0,msg1))
        }
        # Create matching vectors for columns from input annotation
        selectID.gene = mcols(anno.genes)[match(match.geneID.gene, colnames(mcols(anno.genes)))][[1]]
    }
    # test if unique id is present
    if (!uniqueID.gene %in% colnames(mcols(anno.genes))) {
        # id is not present
        # -> throw error
        msg0 = paste0("No unique ID present for genes with name: ", uniqueID.gene, " for option 'uniqueID.gene'. \n")
        msg1 = paste0("Provide a unique ID for each range. \n")
        stop(c(msg0,msg1))
    } else {
        # id is present
        # -> test if it is unique
        this.uniqueID.gene = (mcols(anno.genes))[colnames(mcols(anno.genes)) == uniqueID.gene][[1]]
        if (any(duplicated(this.uniqueID.gene))) {
            msg0 = paste0("Unique IDs are not duplicated for IDs: ",
                          paste(this.uniqueID.gene[duplicated(this.uniqueID.gene)], collapse = ","))
        }
    }

    # Check if binding site to gene assignment was already done
    # -> if not run gene assignment first
    if (!match.geneID.bs %in% colnames(mcols(getRanges(object)))) {
        # no gene ID found in binding site meta columns
        msg0 = paste0("No matching column for gene ID found with name: '", match.geneID.bs, "' in binding site meta data.\n")
        if (!isTRUE(generate.geneID.bs)){
            # binding site to gene matching should not be done
            # -> stop with error
            msg1 = paste0("Binding site to gene matching not actiavted with: generate.geneID.bs=", generate.geneID.bs)
            stop(c(msg0,msg1))
        } else {
            # -> run binding site to gene matching
            msg1 = paste0("Trying 'assignToGenes()' to assing binding sites to hosting genes. \n")
            if(!quiet) warning(msg0, msg1)
            object = assignToGenes(object, anno.genes = anno.genes, quiet = quiet, ...)
        }
    }
    if (!uniqueID.bs %in% colnames(mcols(getRanges(object)))) {
        # id is not present
        # -> throw error
        msg0 = paste0("No unique ID present for binding sites with name: ", uniqueID.bs, " for option 'uniqueID.bs'. \n")
        msg1 = paste0("Provide a unique ID for each range. \n")
        stop(c(msg0,msg1))
    } else {
        # id is present
        # -> test if it is unique
        this.uniqueID.bs = mcols(getRanges(object))[colnames(mcols(getRanges(object))) == uniqueID.bs][[1]]
        if (any(duplicated(this.uniqueID.bs))) {
            msg0 = paste0("Unique IDs are not duplicated for IDs: ",
                          paste(this.uniqueID.bs[duplicated(this.uniqueID.bs)], collapse = ","))
        }
    }

    # Check if blacklist regions have a gene ID that can be used
    # -> if not run gene assignment first
    blacklist.inUse = FALSE
    if (!is.null(blacklist)) {
        # blacklist regions are present
        blacklist.inUse = TRUE
        # -> check gene ID
        if (!match.geneID.blacklist %in% colnames(mcols(blacklist))) {
            # no geneID found in binding site meta columns
            msg0 = paste0("No matching column for gene ID found with name: '", match.geneID.blacklist, "' in blacklist region meta data.\n")
            if (!isTRUE(generate.geneID.blacklist)){
                # binding site to gene matching should not be done
                # -> stop with error
                msg1 = paste0("Blacklist region to gene matching not actiavted with: generate.geneID.blacklist=", generate.geneID.blacklist)
                stop(c(msg0,msg1))
            } else {
                # -> run binding site to gene matching
                msg1 = paste0("Trying 'assignToGenes()' to assing blacklist regions to hosting genes. \n")
                if(!quiet) warning(msg0, msg1)
                # TODO
                # -> think about which parameters for assignToGenes need to be available to the user
                object.blacklist = setRanges(object, blacklist, quiet = quiet)
                object.blacklist = assignToGenes(object.blacklist, anno.genes = anno.genes, quiet = quiet, ...)
                blacklist = getRanges(object.blacklist)
            }
        }
        if (!uniqueID.blacklist %in% colnames(mcols(blacklist))) {
            # id is not present
            # -> throw error
            msg0 = paste0("No unique ID present for blacklist regions with name: ", uniqueID.blacklist, " for option 'uniqueID.blacklist'. \n")
            msg1 = paste0("Provide a unique ID for each range. \n")
            stop(c(msg0,msg1))
        } else {
            # id is present
            # -> test if it is unique
            this.uniqueID.blacklist = mcols(blacklist)[colnames(mcols(blacklist)) == uniqueID.blacklist][[1]]
            if (any(duplicated(this.uniqueID.blacklist))) {
                msg0 = paste0("Unique IDs are not duplicated for IDs: ",
                              paste(this.uniqueID.blacklist[duplicated(this.uniqueID.blacklist)], collapse = ","))
            }
        }
    }

    # get all ranges
    this.ranges.initial = getRanges(object)
    this.ranges = getRanges(object)

    # check ranges width
    if (length(unique(width(this.ranges))) > 1) {
        # binding sites are not all of the same size
        msg0 = paste0("Binding sites are not all of the same size. \n")
        msg1 = paste0("Found sizes: ", paste(unique(width(this.ranges)), collapse = ","))
        msg2 = paste0("Use 'combineBSF()' to unify binding site sizes. \n")
        # check if user wants to force continue despite binding sites being unequal
        if(!isTRUE(force.unequalSites)){
            stop(c(msg0, msg1, msg2))
        } else {
            if (!quiet) warning(c(msg0, msg1, msg2))
        }
    }

    # check for binding sites with geneID = NA and remove
    all.geneIDs = mcols(this.ranges)[match(match.geneID.bs, colnames(mcols(this.ranges)))][[1]]
    if (any(is.na(all.geneIDs))) {
        msg0 = paste0("Found ", length(which(is.na(all.geneIDs))), " binding sites with 'NA' as geneID. \n")
        msg1 = paste0("Removing those keeps ", length(which(!is.na(all.geneIDs))), " binding sites remaining. \n")
        if(!quiet) warning(c(msg0,msg1))
        this.ranges = this.ranges[which(!is.na(all.geneIDs))]
    }

    # final binding site check
    if (length(this.ranges) == 0) {
        # no more binding site ranges are left
        msg0 = paste0("No more binding sites left to work with.\n")
        msg1 = paste0("Check input data. \n")
        stop(c(msg0,msg1))
    }

    # MAIN COMPUTE
    # --------------------------------------------------------------------------

    # calculate coverage per binding site
    # ------------------------------------
    cov.bs = .coverageForDifferential(object, ranges = this.ranges,
                                      match.rangeID = uniqueID.bs,
                                      prefix = "counts.bs.", quiet = veryQuiet)

    # calculate coverage of hosting genes
    # ----------------------------------------
    this.matchIDs.gene = mcols(anno.genes)[,which(match.geneID.gene == colnames(mcols(anno.genes)))]
    this.matchIDs.bs = mcols(this.ranges)[,which(match.geneID.bs == colnames(mcols(this.ranges)))]
    hosting.genes = anno.genes[this.matchIDs.gene %in% this.matchIDs.bs]
    # calculate gene-wise coverage
    cov.hosting.genes = .coverageForDifferential(object, ranges = hosting.genes,
                                                 match.rangeID = uniqueID.gene,
                                                 prefix = "counts.gene.", quiet = veryQuiet)

    # group counts by gene
    cov.bs.merge = as.data.frame(mcols(cov.bs)) %>%
        select(geneID, contains("counts")) %>%
        group_by(geneID) %>%
        dplyr::summarise_all(sum)
    # match binding site counts per gene
    match.gene = mcols(cov.hosting.genes)[colnames(mcols(cov.hosting.genes)) == match.geneID.gene][[1]]

    match.bs = cov.bs.merge[colnames(cov.bs.merge) == match.geneID.bs][[1]]
    idx.bs = match(match.gene, match.bs)
    mcols(cov.hosting.genes) = cbind(mcols(cov.hosting.genes),
                                     select(cov.bs.merge, starts_with("counts"))[idx.bs,])


    # add the offset ranges per gene
    # ----------------------------------------
    if (isTRUE(use.offset)) {
        # offset should be used
        # -> check if the default or manual value should be used
        if(!is.null(ranges.offset)){
            # ranges.offset is set manually
            stopifnot(is(ranges.offset, "numeric"))
        } else {
            ranges.offset = floor(width(this.ranges)/2)
        }

        # turn offset size into offset ranges
        range.offset = this.ranges + ranges.offset

        #
        # # turn offset size into offset ranges
        # range.before = flank(this.ranges, width = ranges.offset, start = TRUE, both = FALSE)
        # range.after = flank(this.ranges, width = ranges.offset, start = FALSE, both = FALSE)
        #
        # # merge overlapping offsets
        # range.offset = c(range.before, range.after)


        range.offset.merge = reduce(range.offset, with.revmap = TRUE)
        idx = unlist(lapply(range.offset.merge$revmap, `[[`, 1))
        range.offset.merge$geneID = mcols(range.offset)[colnames(mcols(range.offset)) == match.geneID.bs][[1]][idx]
        range.offset.merge$id = seq_along(range.offset.merge)

        # calulate coverage on offset regions
        cov.offset = .coverageForDifferential(object, ranges = range.offset.merge,
                                              match.rangeID = "id",
                                              prefix = "counts.offset.", quiet = veryQuiet)

        # group counts by gene
        cov.offset.gene = as.data.frame(mcols(cov.offset)) %>%
            select(all_of(match.geneID.bs), contains("counts")) %>%
            dplyr::group_by_if(., is.character) %>%
            dplyr::summarise_all(sum)

        # match counts per gene
        match.gene = mcols(cov.hosting.genes)[colnames(mcols(cov.hosting.genes)) == match.geneID.gene][[1]]
        match.offset = cov.offset.gene[colnames(cov.offset.gene) == match.geneID.bs][[1]]
        idx.offset = match(match.gene, match.offset)

        # add counts to gene
        mcols(cov.hosting.genes) = cbind(mcols(cov.hosting.genes),
                                         select(cov.offset.gene, starts_with("counts"))[idx.offset,])
    }
    # add additional blacklist regions per gene
    # ----------------------------------------

    if (!is.null(blacklist)) {
        # blacklist region is present
        # -> calculate coverage
        cov.blacklist = .coverageForDifferential(object, ranges = blacklist,
                                                 match.rangeID = uniqueID.blacklist,
                                                 prefix = "counts.black.", quiet = veryQuiet)

        # group counts by gene
        # cov.blacklist = as.data.frame(mcols(cov.blacklist)) %>%
        #     select(geneID, contains("counts")) %>%
        #     group_by(geneID) %>%
        #     dplyr::summarise_all(sum)

        cov.blacklist = as.data.frame(mcols(cov.blacklist)) %>%
            select(all_of(match.geneID.blacklist), contains("counts")) %>%
            dplyr::group_by_if(., is.character) %>%
            dplyr::summarise_all(sum)

        # match counts per gene
        match.gene = mcols(cov.hosting.genes)[colnames(mcols(cov.hosting.genes)) == match.geneID.gene][[1]]

        match.blacklist = cov.blacklist[colnames(cov.blacklist) == match.geneID.bs][[1]]
        idx.blacklist = match(match.gene, match.blacklist)
        mcols(cov.hosting.genes) = cbind(mcols(cov.hosting.genes),
                                         select(cov.blacklist, starts_with("counts"))[idx.blacklist,])
    }

    # manage returns
    # ----------------------------------------
    countMatr = mcols(cov.hosting.genes) %>% as.data.frame() %>%
        select(starts_with("counts"), ) %>%
        replace(is.na(.), 0)

    if (isTRUE(use.offset) | !is.null(blacklist)) {
        # one additional resource to subtract from gene counts is present
        if (isTRUE(use.offset) & !is.null(blacklist)) {
            # offset and blacklist are used
            # -> subtract both with bs from gene counts
            sel.countsToRemove = countMatr %>% select(contains(c("offset", "black")))

            # combine counts to remove
            sample.ids = paste0(this.meta$id, "_", this.meta$condition)
            sel.countsToRemove = as.matrix(sel.countsToRemove)
            countsToRemove = lapply(sample.ids, function(x) {
                rowSums(sel.countsToRemove[, grepl(x, colnames(sel.countsToRemove))])
            })
            countsToRemove = do.call(cbind, countsToRemove)
        }
        if (isTRUE(use.offset) & is.null(blacklist)) {
            # offset is used and blacklist is not used
            # -> subtract offset with bs from gene counts
            sel.countsToRemove = countMatr %>% select(contains(c("offset")))
            countsToRemove = sel.countsToRemove
        }
        if (!isTRUE(use.offset) & !is.null(blacklist)) {
            # offset is not used and blacklist is used
            # -> subtract blacklist with bs from gene counts
            sel.countsToRemove = countMatr %>% select(contains(c("black")))
            countsToRemove = sel.countsToRemove
        }
    } else {
        # only one type signal is used for background
        if (!isTRUE(use.offset) & is.null(blacklist)) {
            # none of the options are used (no blacklist, not offset)
            # -> only subtract bs count from genes
            sel.countsToRemove = countMatr %>% select(contains(c("bs")))
        }
        countsToRemove = sel.countsToRemove
    }

    # make background counts
    sample.ids = paste0(this.meta$id, "_", this.meta$condition)
    countsBg = countMatr %>% select(contains("gene"))
    countsBg = countsBg - countsToRemove
    colnames(countsBg) = paste0("counts.bg.",sample.ids)

    # select background and gene coverage for each hosting gene and replicate
    # mcols(cov.hosting.genes) = cbind(select(as.data.frame(mcols(cov.hosting.genes)), -(starts_with("counts"))),
    #                                  select(as.data.frame(mcols(cov.hosting.genes)), starts_with("counts.gene")),
    #                                  countsBg)
    # select background coverage for each hosting gene and replicate
    mcols(cov.hosting.genes) = cbind(select(as.data.frame(mcols(cov.hosting.genes)), -(starts_with("counts"))),
                                     countsBg)

    # match hosting gene counts with binding site counts
    this.matchIDs.gene = mcols(cov.hosting.genes)[,which(match.geneID.gene == colnames(mcols(cov.hosting.genes)))]
    idx = match(this.matchIDs.bs, this.matchIDs.gene)
    mcols(cov.bs) = cbind(mcols(cov.bs),
                          select(as.data.frame(mcols(cov.hosting.genes)), starts_with("counts"))[idx,])

    # set object and return
    object = setRanges(object, cov.bs, quiet = veryQuiet)

    # ---
    # Store function parameters in list
    optstr = list(source = datasource,
                  blacklist = blacklist.inUse,
                  use.offset = use.offset,
                  ranges.offset = paste(unique(ranges.offset), sep = ","),
                  match.geneID.gene = match.geneID.gene,
                  match.geneID.bs = match.geneID.bs,
                  match.geneID.blacklist = match.geneID.blacklist,
                  generate.geneID.blacklist = generate.geneID.blacklist,
                  generate.geneID.bs = generate.geneID.bs,
                  uniqueID.blacklist = uniqueID.blacklist,
                  uniqueID.bs = uniqueID.bs,
                  uniqueID.gene = uniqueID.gene,
                  force.unequalSites = force.unequalSites)
    object@params$calculateBsBackground = optstr

    # ---
    # Store for results
    resultLine = data.frame(
        funName = "calculateBsBackground()", class = "transform",
        nIn = length(this.ranges.initial), nOut = length(cov.bs),
        per = paste0(round(length(cov.bs)/ length(this.ranges.initial), digits = 2)*100,"%"),
        options = paste0("Ranges.offset=", optstr$ranges.offset,
                         ", source=", optstr$datasource,
                         ", blacklist=", optstr$blacklist,
                         ", use.offset=", optstr$use.offset,
                         ", match.geneID.gene=", optstr$match.geneID.gene,
                         ", match.geneID.bs=", optstr$match.geneID.bs,
                         ", match.geneID.blacklist=", optstr$match.geneID.blacklist,
                         ", generate.geneID.blacklist=", optstr$generate.geneID.blacklist,
                         ", generate.geneID.bs=", optstr$generate.geneID.bs,
                         ", uniqueID.blacklist=", optstr$uniqueID.blacklist,
                         ", uniqueID.bs=", optstr$uniqueID.bs,
                         ", uniqueID.gene=", optstr$uniqueID.gene,
                         ", force.unequalSites=", optstr$force.unequalSites)
    )

    object@results = rbind(object@results, resultLine)

    return(object)
}


#' Filter for genes not suitable for differential testing
#'
#' This function removes genes where the differential testing protocol can
#' not be applied to, using count coverage information on the binding sites
#' and background regions per gene, through the following steps:
#' \enumerate{
#' \item Remove genes with overall not enough crosslinks: \code{minCounts}
#' \item Remove genes with a disproportion of counts in binding sites vs. the
#' background: \code{balanceBackground}
#' \item Remove genes where the expression between conditions is too much off
#' balance: \code{balanceCondition}
#' }
#'
#' To remove genes with overall not enough crosslinks (\code{minCounts})
#' all counts are summed up per gene across all samples and compared to the
#' minimal count threshold (\code{minCounts.cutoff}).
#'
#' To remove genes with a count disproportion between binding sites and
#' background regions crosslinks are summed up for binding sites and background
#' per gene. These sums are combined in a ratio. Genes where eg. 50\% of all
#' counts are within binding sites would be removed
#' (see \code{balanceBackground.cutoff.bs} and \code{balanceBackground.cutoff.bg}).
#'
#' To remove genes with very large expression differences between conditions,
#' crosslinks are summed up per gene for each condition. If now eg. the total
#' number of crosslinks is for 98\% in one condition and only 2\% of the
#' combined signal is in the second condition, expression levels are too
#' different for a reliable comparisson (see \code{balanceCondition.cutoff}).
#'
#' @param object a \code{\link{BSFDataSet}} object with computed count data
#' for binding sites and background regions
#' @param minCounts logical; whether to use the minimum count filter
#' @param minCounts.cutoff numeric; the minimal number of crosslink
#' per gene over all samples for the gene to be retained (default = 100)
#' @param balanceBackground logical; whether to use the counts balancing filter
#' between binding sites and background
#' @param balanceBackground.cutoff.bs numeric; the maximum fraction of the total signal
#' per gene that can be within binding sites (default = 0.2)
#' @param balanceBackground.cutoff.bg numeric; the minimum fraction of the total signal
#' per gene that can be within the background (default = 0.8)
#' @param balanceCondition logical; whether to use the counts balancing filter
#' between conditions
#' @param balanceCondition.cutoff numeric; the maximum fraction of the total signal
#' that can be attributed to only one condition
#' @param match.geneID character; the name of the column with the gene ID
#' in the binding sites meta columns used for matching binding sites to genes
#' @param flag logical; whether to remove or flag binding sites from genes that do not
#' pass any of the filters
#' @param quiet logical; whether to print messages or not
#' @param veryQuiet logical; whether to print messages or not
#'
#' @return an object of class \code{\link{BSFDataSet}} with biniding sites filtered
#' or flagged by the above filter options
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#'
#' # make testset
#' bds = makeBindingSites(bds, bsSize = 7)
#' bds = assignToGenes(bds, anno.genes = gns)
#' bds = imputeBsDifferencesForTestdata(bds)
#' bds = calculateBsBackground(bds, anno.genes = gns, use.offset = FALSE)
#'
#' # use all filters and remove binding sites that fail (default settings)
#' f0 = filterBsBackground(bds)
#'
#' # do not use the condition balancing filter
#' f1 = filterBsBackground(bds, balanceCondition = FALSE)
#'
#' # use only the minimum count filter and flag binding sites instead of
#' # removing them
#' f3 = filterBsBackground(bds, flag = TRUE, balanceCondition = FALSE,
#'  balanceBackground = FALSE)
#'
#' @export
filterBsBackground <- function(object,
                               minCounts = TRUE,
                               minCounts.cutoff = 100,
                               balanceBackground = TRUE,
                               balanceBackground.cutoff.bs = 0.3,
                               balanceBackground.cutoff.bg = 0.7,
                               balanceCondition = TRUE,
                               balanceCondition.cutoff = 0.05,
                               match.geneID = "geneID",
                               flag = FALSE,
                               quiet = FALSE,
                               veryQuiet = FALSE
){
    # Initialize local variables
    . <- s <- name <- value <- type <- bg <- gene.approx <- bs <- ratio.bg <- NULL
    ratio.bs <- f.bg <- f.bs <-condition <- total <- both <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is.logical(quiet))
    stopifnot(is.logical(veryQuiet))
    stopifnot(is.logical(minCounts))
    stopifnot(is.logical(balanceBackground))
    stopifnot(is.logical(balanceCondition))
    stopifnot(is.logical(flag))

    # get ranges
    this.ranges.initial = getRanges(object)
    this.ranges = getRanges(object)

    # check if ranges have geneID present
    if (!match.geneID %in% colnames(mcols(this.ranges))) {
        msg0 = paste0("No matching column for gene ID found with name: '", match.geneID, "' in binding site meta data.\n")
        msg1 = paste0("Please run calculateBsBackground() first. \n")
        stop(c(msg0, msg1))
    }

    # check if ranges have counts present
    if (!any(grepl("counts", colnames(as.data.frame(mcols(this.ranges)))))) {
        # meta data has no columns named counts
        # -> stop here; point user to calculateBsBackground
        msg0 = paste0("No associated fields with crosslink counts found for present binding sites. \n")
        msg1 = paste0("Make sure to run 'calculateBsBackground()' first. \n")
        stop(msg0,msg1)
    }
    # check if ranges have at least two samples per condition in counts
    this.counts = as.data.frame(mcols(this.ranges)) %>% select(contains("counts"))
    # check binding sites
    this.counts.bs = this.counts %>% select(contains("bs"))
    if (length(this.counts.bs) == 0) {
        msg0 = paste0("No counts for binding sites present. \n")
        msg1 = paste0("Make sure to run 'calculateBsBackground()' first. \n")
        stop(c(msg0, msg1))
    }
    # check background
    this.counts.bg = this.counts %>% select(contains("bg"))
    if (length(this.counts.bg) == 0) {
        msg0 = paste0("No counts for background present. \n")
        msg1 = paste0("Make sure to run 'calculateBsBackground()' first. \n")
        stop(c(msg0, msg1))
    }

    this.meta = getMeta(object)
    this.condition.reference = levels(this.meta$condition)[1]
    this.condition.comp = levels(this.meta$condition)[2]

    # MAIN
    # --------------------------------------------------------------------------

    # filter for genes with less counts than cutoff
    if (isTRUE(minCounts)) {
        # the gene-wise count filter should be used
        # -> group counts by gene and sample and filter
        all.counts = as.data.frame(mcols(this.ranges)) %>%
            select(all_of(match.geneID), contains("counts")) %>%
            dplyr::group_by_at(match.geneID) %>%
            summarise(dplyr::across(everything(), sum)) %>%
            mutate(s = select(., contains("counts")) %>% rowSums(.)) %>%
            select(all_of(match.geneID), s)
        idx = which(all.counts < minCounts.cutoff, arr.ind = TRUE)

        # ---
        # Store for plotting
        object@plotData$filterBsBackground$data.minCounts = all.counts


        if (length(idx) > 0) {
            if (!veryQuiet) message("Count filter ")
            # found positions with counts less than cutoff
            # remove these counts
            msg0 = paste0("Found ", format(nrow(idx), big.mark = ",", decimal.mark = "."),
                          " genes with less total counts than ", minCounts.cutoff, ".\n")
            # remove genes below threshold
            genesToRemove = all.counts[idx[,1],] %>% pull(match.geneID)
            allGenes = mcols(this.ranges)[colnames(mcols(this.ranges)) == match.geneID][[1]]
            idxToRemove = allGenes %in% genesToRemove

            if (isTRUE(flag)) {
                # flag genes with low counts
                # -> add additional meta column
                msg1 = paste0("Flagging a total of ",
                              format(length(genesToRemove), big.mark = ",", decimal.mark = "."),
                              " genes and ",
                              format(sum(idxToRemove), big.mark = ",", decimal.mark = "."),
                              " binding sites.\n")
                mcols(this.ranges)$flag.minCount = idxToRemove
            }
            if (!isTRUE(flag)) {
                # remove genes with low counts
                msg1 = paste0("Removing a total of ",
                              format(length(genesToRemove), big.mark = ",", decimal.mark = "."),
                              " genes and ",
                              format(sum(idxToRemove), big.mark = ",", decimal.mark = "."),
                              " binding sites.\n")
                this.ranges = this.ranges[!idxToRemove]
            }
            # inform user
            if (!quiet) warning(c(msg0, msg1))
        }
    }

    # filter for genes with binding site vs background ratio being off
    if (isTRUE(balanceBackground)) {
        if (!veryQuiet) message("Ratio filter ")
        # binding site to background ratio filter is used
        # -> calculate ratios per gene over all samples
        df.ratio = as.data.frame(mcols(this.ranges)) %>%
            select(all_of(match.geneID), starts_with("counts")) %>%
            pivot_longer(-all_of(match.geneID)) %>%
            mutate(type = sub(".*\\.(.*?)\\..*", "\\1", name)) %>%
            select(-name) %>%
            dplyr::group_by_if(., is.character) %>%
            summarise(sum = sum(value), .groups = "keep") %>%
            pivot_wider(values_from = sum, names_from = type) %>%
            mutate(gene.approx = bg + bs) %>%
            mutate(ratio.bg = bg / (gene.approx + 1)) %>%
            mutate(ratio.bs = bs / (gene.approx + 1)) %>%
            select(all_of(match.geneID), ratio.bg, ratio.bs) %>%
            pivot_longer(-all_of(match.geneID))

        # ---
        # Store for plotting
        object@plotData$filterBsBackground$data.balanceBackground = df.ratio

        # apply filter cutoffs
        df.filter = df.ratio %>% ungroup() %>%
            pivot_wider(names_from = name, values_from = value) %>%
            mutate(f.bg = ifelse(ratio.bg < balanceBackground.cutoff.bs, TRUE, FALSE),
                   f.bs = ifelse(ratio.bs > balanceBackground.cutoff.bs, TRUE, FALSE))

        msg0 = paste0("Found ", format(sum(df.filter$f.bg | df.filter$f.bs), big.mark = ".", decimal.mark = ","),
                      " genes with bs to bg ratio not meeting thresholds.\n")

        genesToRemove = df.filter %>% filter(f.bg | f.bs) %>% pull(match.geneID)
        allGenes = mcols(this.ranges)[colnames(mcols(this.ranges)) == match.geneID][[1]]
        idxToRemove = allGenes %in% genesToRemove

        if (isTRUE(flag)) {
            # flag genes with ratio above thresholds
            # -> add additional meta column
            msg1 = paste0("Flagging a total of ",
                          format(length(genesToRemove), big.mark = ",", decimal.mark = "."),
                          " genes and ",
                          format(sum(idxToRemove), big.mark = ",", decimal.mark = "."),
                          " binding sites.\n")
            mcols(this.ranges)$flag.ratio = idxToRemove
        }
        if (!isTRUE(flag)) {
            # remove genes with low counts
            msg1 = paste0("Removing a total of ",
                          format(length(genesToRemove), big.mark = ",", decimal.mark = "."),
                          " genes and ",
                          format(sum(idxToRemove), big.mark = ",", decimal.mark = "."),
                          " binding sites.\n")
            this.ranges = this.ranges[!idxToRemove]
        }
        # inform user
        if (!quiet) warning(c(msg0, msg1))
    }

    # filter for genes with large count imbalances between both conditions
    if (isTRUE(balanceCondition)) {
        if (!veryQuiet) message("Balance filter ")
        # genes should be filtered for count differences between conditions
        # -> calculate sum of counts per gene and condition
        df.ratio = as.data.frame(mcols(this.ranges)) %>%
            select(all_of(match.geneID), starts_with("counts")) %>%
            pivot_longer(-all_of(match.geneID)) %>%
            mutate(type = sub(".*\\.(.*?)\\..*", "\\1", name)) %>%
            mutate(condition = sub(".*_", "", name)) %>%
            filter(type == "bg") %>%
            select(-c(name, type)) %>%
            dplyr::group_by_if(., is.character) %>%
            summarise(sum = sum(value), .groups = "keep") %>%
            pivot_wider(values_from = sum, names_from = condition) %>%
            mutate(total = rowSums(dplyr::across(dplyr::where(is.numeric)))) %>%
            ungroup() %>%
            dplyr::relocate(all_of(match.geneID), total, all_of(this.condition.reference), all_of(this.condition.comp))

        # ---
        # Store for plotting
        object@plotData$filterBsBackground$data.balanceCondition = df.ratio

        df.filter = data.frame(
            geneID = df.ratio %>% select(all_of(match.geneID)),
            this.condition.reference = ifelse(c(df.ratio[,3][[1]] / df.ratio[,2][[1]]) < balanceCondition.cutoff, TRUE, FALSE),
            his.condition.comp = ifelse(c(df.ratio[,4][[1]] / df.ratio[,2][[1]]) < balanceCondition.cutoff, TRUE, FALSE)
        )
        df.filter$both = df.filter[,2] | df.filter[,3]

        msg0 = paste0("Found ", format(sum(df.filter$both, na.rm = TRUE), big.mark = ".", decimal.mark = ","),
                      " genes with condition balance ratio not meeting thresholds.\n")

        genesToRemove = df.filter %>% filter(both) %>% pull(match.geneID)
        allGenes = mcols(this.ranges)[colnames(mcols(this.ranges)) == match.geneID][[1]]
        idxToRemove = allGenes %in% genesToRemove

        if (isTRUE(flag)) {
            # flag genes with balance ratio above thresholds
            # -> add additional meta column
            msg1 = paste0("Flagging a total of ",
                          format(length(genesToRemove), big.mark = ",", decimal.mark = "."),
                          " genes and ",
                          format(sum(idxToRemove), big.mark = ",", decimal.mark = "."),
                          " binding sites.\n")
            mcols(this.ranges)$flag.balance = idxToRemove
        }
        if (!isTRUE(flag)) {
            # remove genes with low counts
            msg1 = paste0("Removing a total of ",
                          format(length(genesToRemove), big.mark = ",", decimal.mark = "."),
                          " genes and ",
                          format(sum(idxToRemove), big.mark = ",", decimal.mark = "."),
                          " binding sites.\n")
            this.ranges = this.ranges[!idxToRemove]
        }
        # inform user
        if (!quiet) warning(c(msg0, msg1))
    }

    object = setRanges(object, this.ranges, quiet = quiet)

    # ---
    # Store function parameters in list
    optstr = list(
        minCounts = minCounts,
        minCounts.cutoff = minCounts.cutoff,
        balanceBackground = balanceBackground,
        balanceBackground.cutoff.bs = balanceBackground.cutoff.bs,
        balanceBackground.cutoff.bg = balanceBackground.cutoff.bg,
        balanceCondition = balanceCondition,
        balanceCondition.cutoff = balanceCondition.cutoff,
        match.geneID = match.geneID,
        flag = flag)
    object@params$filterBsBackground = optstr

    # ---
    # Store for results
    resultLine = data.frame(
        funName = "filterBsBackground()", class = "transform",
        nIn = length(this.ranges.initial), nOut = length(this.ranges),
        per = paste0(round(length(this.ranges)/ length(this.ranges.initial), digits = 2)*100,"%"),
        options = paste0("minCounts=", optstr$minCounts,
                         ", source=", optstr$datasource,
                         ", minCounts.cutoff=", optstr$minCounts.cutoff,
                         ", balanceBackground=", optstr$balanceBackground,
                         ", balanceBackground.cutoff.bs=", optstr$balanceBackground.cutoff.bs,
                         ", balanceBackground.cutoff.bg=", optstr$balanceBackground.cutoff.bg,
                         ", balanceCondition=", optstr$balanceCondition,
                         ", balanceCondition.cutoff=", optstr$balanceCondition.cutoff,
                         ", match.geneID=", optstr$generate.genmatch.geneIDeID.bs,
                         ", flag=", optstr$flag)
    )
    object@results = rbind(object@results, resultLine)

    return(object)
}






#' Compute fold-changes per binding site
#'
#' Given count data for binding sites and background regions this function will
#' compute fold-changes between two condition for each binding site. Computation
#' is based on \code{\link[DESeq2]{DESeq}} using the Likelihood ratio test to
#' disentangle transcription level changes from binding site level changes.
#'
#' Fold-changes per binding sites are corrected for transcript level changes.
#' Essentially, background counts are used to model transcript level changes,
#' which are then used to compute fold-changes per binding site, which are
#' corrected for the observed transcript level changes. This is done by using a
#' Likelihood ratio test comparing the full model
#' (~condition + type + condition:type) to the reduced model (~condition + type).
#'
#' Fold-changes for the transcript level changes are modeled explicitly in a
#' second round of the \code{\link[DESeq2]{DESeq}} workflow, using the default
#' Wald test to compare changes between the conditions (~condition). Counts
#' attributed to binding sites are removed from the gene level counts.
#'
#' Results from both calculation rounds can be filtered and further manipulated
#' with parameters given in the DESeq2 framework
#' (see \code{\link[DESeq2]{results}}, \code{\link[DESeq2]{lfcShrink}}).
#'
#' @param object a \code{\link{BSFDataSet}} object
#' @param fitType either "parametric", "local", "mean", or "glmGamPoi" for the
#' type of fitting of dispersions to the mean intensity.
#' See \code{\link[DESeq2]{DESeq}} for more details.
#' @param sfType either "ratio", "poscounts", or "iterate" for the type of size
#' factor estimation. See \code{\link[DESeq2]{DESeq}} for more details.
#' @param minReplicatesForReplace the minimum number of replicates required in
#' order to use replaceOutliers on a sample. See \code{\link[DESeq2]{DESeq}} for
#' more details.
#' @param independentFiltering logical, whether independent filtering should be
#' applied automatically. See \code{\link[DESeq2]{results}} for more details.
#' @param alpha the significance cutoff used for optimizing the independent
#' filtering. See \code{\link[DESeq2]{results}} for more details.
#' @param pAdjustMethod he method to use for adjusting p-values. See
#' \code{\link[DESeq2]{results}} for more details.
#' @param minmu lower bound on the estimated count. See
#' \code{\link[DESeq2]{results}} for more details.
#' @param filterFun an optional custom function for performing independent
#' filtering and p-value adjustment. See \code{\link[DESeq2]{results}} for more
#' details.
#' @param use.lfc.shrinkage logical; whether to compute shrunken log2 fold
#' changes for the DESeq results. See \code{\link[DESeq2]{lfcShrink}} for more
#' details.
#' @param type if 'ashr', 'apeglml' or 'normal' should be used for fold change
#' shrinkage. See \code{\link[DESeq2]{lfcShrink}} for more details.
#' @param svalue logical, should p-values and adjusted p-values be replaced with
#' s-values when using apeglm or ashr. See \code{\link[DESeq2]{lfcShrink}}
#' for more details.
#' @param apeAdapt logical, should apeglm use the MLE estimates of LFC to adapt
#' the prior. See \code{\link[DESeq2]{lfcShrink}} for more details.
#' @param apeMethod what method to run apeglm, which can differ in terms of
#' speed. See \code{\link[DESeq2]{lfcShrink}} for more details.
#' @param match.geneID character; the name of the column with the gene ID
#' in the binding sites meta columns used for matching binding sites to genes
#' @param quiet logical; whether to print messages or not
#' @param veryQuiet logical; whether to print messages or not
#' @param replaceNegative logical; force negative counts to be replaces by 0.
#' Be careful when using this, having negative counts can point towards problems
#' with the gene annotation in use.
#' @param removeNA logical; force binding sites with any NA value to be removed.
#' Be careful when using this, having negative counts can point towards problems
#' with the gene annotation in use.
#'
#' @return a \code{\link{BSFDataSet}} object, with results from the
#' \code{\link[DESeq2]{DESeq}} analysis added to the meta columns of the
#' binding site ranges.
#'
#' @seealso \code{\link{calculateBsBackground}} \code{\link[DESeq2]{DESeq}}
#'
#' @importFrom stats relevel
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#'
#' # make example dataset
#' bds = makeBindingSites(bds, bsSize = 7)
#' bds = assignToGenes(bds, anno.genes = gns)
#' bds = imputeBsDifferencesForTestdata(bds)
#' bds = calculateBsBackground(bds, anno.genes = gns)
#'
#' # calculate fold changes - no shrinkage
#' bds = calculateBsFoldChange(bds)
#'
#' # calculate fold changes - with shrinkage
#' bds = calculateBsFoldChange(bds, use.lfc.shrinkage = TRUE)
#'
#' @export
calculateBsFoldChange <- function(object,
                                  # for DESeq
                                  fitType = "local",
                                  sfType = "ratio",
                                  minReplicatesForReplace = 10,
                                  # for results
                                  independentFiltering = TRUE,
                                  alpha = 0.05,
                                  pAdjustMethod = "BH",
                                  minmu = 0.5,
                                  filterFun = NULL,
                                  # for lfcShrink
                                  use.lfc.shrinkage = FALSE,
                                  type = c("ashr", "apeglm", "normal"),
                                  svalue = FALSE,
                                  apeAdapt = TRUE,
                                  apeMethod = "nbinomCR",
                                  # general
                                  match.geneID = "geneID",
                                  quiet = TRUE,
                                  veryQuiet = FALSE,
                                  replaceNegative = FALSE,
                                  removeNA = FALSE
                                  # forceThis = FALSE
){
    # Bind locale variables
    geneID <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is.logical(quiet))
    stopifnot(is.logical(veryQuiet))
    stopifnot(is.logical(use.lfc.shrinkage))

    # check if shrinkage type is set
    if (isTRUE(use.lfc.shrinkage)) {
        type = match.arg(type, choices = c("ashr", "apeglm", "normal"))
    }

    this.meta = getMeta(object)
    # check for correct pairwise setting in meta data
    if (length(levels(this.meta$condition)) != 2) {
        msg0 = paste0("Fold-changes can only be computed for pairwise comparisons.\n")
        if (length(levels(this.meta$condition)) < 2) {
            # only one condition found
            msg1 = paste0("Found that only one condition is present. There is nothing to compare.\n")
            msg2 = paste0("Please provide two conditions that can be compared. \n")
        }
        if (length(levels(this.meta$condition)) > 2) {
            msg1 = paste0("Found more than 2 conditions. Too many comparisons.\n")
            msg2 = paste0("Please split mulitple comparisons into sets of pairwise comparisions.\n")
        }
        stop(c(msg0,msg1,msg2))
    }

    this.ranges = getRanges(object)

    # check if ranges have geneID present
    if (!match.geneID %in% colnames(mcols(this.ranges))) {
        msg0 = paste0("No matching column for gene ID found with name: '", match.geneID, "' in binding site meta data.\n")
        msg1 = paste0("Please run calculateBsBackground() first. \n")
        stop(c(msg0, msg1))
    }

    # check if ranges have counts present
    if (!any(grepl("counts", colnames(as.data.frame(mcols(this.ranges)))))) {
        # meta data has no columns named counts
        # -> stop here; point user to calculateBsBackground
        msg0 = paste0("No associated fields with crosslink counts found for present binding sites. \n")
        msg1 = paste0("Make sure to run 'calculateBsBackground()' first. \n")
        stop(msg0,msg1)
    }

    # check if ranges have at least two samples per condition in counts
    this.counts = as.data.frame(mcols(this.ranges)) %>% select(contains("counts"))
    # check binding sites
    this.counts.bs = this.counts %>% select(contains("bs"))
    if (length(this.counts.bs) == 0) {
        msg0 = paste0("No counts for binding sites present. \n")
        msg1 = paste0("Make sure to run 'calculateBsBackground()' first. \n")
        stop(c(msg0,msg1))
    }
    this.counts.bs.samples = as.data.frame(table(sub(".*_", "", colnames(this.counts.bs))))
    if (!all(this.counts.bs.samples[,2] > 1)) {
        # found less than two replicates per condition for bs
        msg0 = paste0("Less than two replicates are present for each condition. \n")
        msg1 = paste0("Found conditions: ", paste(this.counts.bs.samples[,1], collapse = ", "),
                      " with counts: ", paste(this.counts.bs.samples[,2], collapse = ", "), "\n")
        stop(c(msg0, msg1))
    }
    # check background
    this.counts.bg = this.counts %>% select(contains("bg"))
    if (length(this.counts.bg) == 0) {
        msg0 = paste0("No counts for background present. \n")
        msg1 = paste0("Make sure to run 'calculateBsBackground()' first. \n")
        stop(c(msg0,msg1))
    }
    this.counts.bg.samples = as.data.frame(table(sub(".*_", "", colnames(this.counts.bg))))
    if (!all(this.counts.bg.samples[,2] > 1)) {
        # found less than two replicates per condition for bs
        msg0 = paste0("Less than two replicates are present for each condition. \n")
        msg1 = paste0("Found conditions: ", paste(this.counts.bg.samples[,1], collapse = ", "),
                      " with counts: ", paste(this.counts.bg.samples[,2], collapse = ", "), "\n")
        stop(c(msg0, msg1))
    }

    # TODO
    # ----------------------------------------------------------------
    # check for NA values
    idx = which(is.na(this.counts), arr.ind = TRUE)
    if (nrow(idx) > 0) {
        # found rows with NA values as counts in matrix
        # -> prompt user towards removing or checking
        msg0 = paste0(format(nrow(idx), big.mark = ",", decimal.mark = "."),
                      " ranges found with NA values.\n")
        if (isTRUE(removeNA)) {
            # force remove NA cases
            this.counts = this.counts[-idx[,1],]
            this.ranges = this.ranges[-idx[,1]]
            msg1 = paste0("Removed all cases.\n")
            if (!veryQuiet) warning(c(msg0, msg1))
        } else {
            msg1 = paste0("This can point towards problems with the gene annoation.\n")
            msg2 = paste0("You can either
                      (1) find problematic genes and fix the cause, or
                      (2) set 'removeNA=TRUE' to force remove those cases.\n")
            stop(c(msg0,msg1,msg2))
        }
    }

    # check for negative values
    idx = which(this.counts < 0, arr.ind = TRUE)
    if (nrow(idx) > 0) {
        # found rows with negative values as counts in matrix
        # -> prompt user towards removing or checking
        msg0 = paste0(format(nrow(idx), big.mark = ",", decimal.mark = "."),
                      " ranges found with negative values.\n")
        if (isTRUE(replaceNegative)) {
            # force replace negative counts
            this.counts[idx] = 0
            msg1 = paste0("Replaced all cases with 0.\n")
            if (!veryQuiet) warning(c(msg0, msg1))
        } else {
            msg1 = paste0("This can point towards problems with the gene annoation.\n")
            msg2 = paste0("You can either
                      (1) find problematic genes and fix the cause, or
                      (2) set 'replaceNegative=TRUE' to force them to 0.\n")
            stop(c(msg0,msg1,msg2))
        }
    }

    # if (isTRUE(forceThis)) {
    #     # TODO check if there are negative counts
    #     # -> temporary hack to replace negative counts with 0
    #     # -> this happens in rare instances where offset ranges of neighboring
    #     #    binding sites overlap, which cause crosslinks to be counted twice.
    #     #    On genes with low counts this could result in negative values for the
    #     #    background.
    #     # -> Optimal solution is to use a combination of `reduce + with.revmap` and
    #     #    `disjoin` to split overlapping offset ranges, which avoids the problem
    #     idx = which(this.counts < 0, arr.ind = TRUE)
    #     if (nrow(idx) > 0) {
    #         this.counts[idx] = 0
    #         msg0 = paste0(format(nrow(idx), big.mark = ",", decimal.mark = "."),
    #                       " ranges found with negative counts, replacing them with 0.\n")
    #         if (!veryQuiet) warning(c(msg0))
    #     }
    #
    #     # TODO
    #     # quick fix for NA values
    #     idx = which(is.na(this.counts), arr.ind = TRUE)
    #     if (nrow(idx) > 0) {
    #         this.counts = this.counts[-idx[,1],]
    #         this.ranges = this.ranges[-idx[,1]]
    #
    #         if (!veryQuiet) warning(c(msg0))
    #     }
    # }

    # ----------------------------------------------------------------
    # TODO



    # --------------------------------------------------------------------------
    # MAIN 1) - Testing binding sites
    # --------------------------------------------------------------------------
    if (!veryQuiet) message("1) Testing binding sites ")

    # construct coldata
    # -----------------
    coldata = rbind.data.frame(this.meta, this.meta)
    coldata$clPlus = NULL
    coldata$clMinus = NULL
    coldata$type = as.factor(c(rep("bs", nrow(this.meta)), rep("bg", nrow(this.meta))))

    # make counts matrix and pass on to deseq
    # ---------------------------------------
    count.matr = this.counts %>% select(contains(c("bs", "bg"))) %>% as.matrix()
    se = SummarizedExperiment::SummarizedExperiment(assays = list(counts = count.matr),
                                                    rowRanges = granges(this.ranges),
                                                    colData = coldata)
    mode(SummarizedExperiment::assay(se)) = "integer"
    dds = DESeq2::DESeqDataSet(se, design = ~condition + type + condition:type)

    # manage model parametrisation through factor levels
    # -> reference level is taken from meta data
    # -> is the first level from condition
    reference.level = levels(this.meta$condition)[1] # WT
    dds$condition = stats::relevel(dds$condition, reference.level)
    dds$type = stats::relevel(dds$type, "bg")

    # running the deseq2 test
    # -----------------------
    dds = DESeq2::DESeq(dds, test="LRT", reduced = ~ condition + type,
                        fitType = fitType,
                        sfType = sfType,
                        quiet = quiet,
                        minReplicatesForReplace = minReplicatesForReplace,
                        modelMatrixType = "standard",
                        useT = FALSE)

    # pulling deseq2 results
    # ----------------------
    comp.level = levels(this.meta$condition)[2] # KO
    this.contrast = paste0("condition", comp.level, ".typebs")

    if (!this.contrast %in% DESeq2::resultsNames(dds)) {
        msg0 = paste0("Contrast named: ", this.contrast, "not found in available names: ", paste(DESeq2::resultsNames(dds), collapse = ", "), ".\n")
        msg1 = paste0("Please check the factor levels of the 'condition' column in your meta data.\n")
        stop(c(msg0, msg1))
    }

    if (!is.null(filterFun)) {
        res = DESeq2::results(dds, name = this.contrast,
                              lfcThreshold = 0,
                              altHypothesis = "greaterAbs",
                              independentFiltering = independentFiltering,
                              alpha = alpha, pAdjustMethod = pAdjustMethod,
                              filterFun = filterFun, format = "DataFrame",
                              minmu = minmu)
    } else {
        res = DESeq2::results(dds, name = this.contrast,
                              lfcThreshold = 0,
                              altHypothesis = "greaterAbs",
                              independentFiltering = independentFiltering,
                              alpha = alpha, pAdjustMethod = pAdjustMethod,
                              format = "DataFrame", minmu = minmu)
    }

    # apply fold change shrinkage if needed
    # -> requiers additional packages to be installed
    if (isTRUE(use.lfc.shrinkage)) {
        res = DESeq2::lfcShrink(dds = dds, res = res, coef = this.contrast,
                                type = type, lfcThreshold = 0,
                                svalue = svalue, returnList = FALSE,
                                format = "DataFrame", saveCols = NULL,
                                apeAdapt = apeAdapt, apeMethod = apeMethod,
                                quiet = quiet)
    }

    # match with ranges
    colnames(res) = paste0('bs.', colnames(res))
    mcols(this.ranges) = cbind.data.frame(mcols(this.ranges), res)


    # --------------------------------------------------------------------------
    # MAIN 2) - Testing the background
    # --------------------------------------------------------------------------
    if (!veryQuiet) message("2) Testing genes ")

    # construct coldata
    # -----------------
    coldata.bg = this.meta
    coldata.bg$clPlus = NULL
    coldata.bg$clMinus = NULL

    # make counts matrix and pass on to deseq
    # ---------------------------------------
    count.matr.bg = as.data.frame(count.matr)
    count.matr.bg$geneID = mcols(this.ranges)[colnames(mcols(this.ranges)) == match.geneID][[1]]
    count.matr.bg = count.matr.bg %>% select(geneID, starts_with("counts.bg"))
    count.matr.bg = count.matr.bg[!duplicated(count.matr.bg$geneID),]
    rownames(count.matr.bg) = count.matr.bg$geneID
    count.matr.bg$geneID = NULL

    se.bg = SummarizedExperiment::SummarizedExperiment(assays = list(counts = as.matrix(count.matr.bg)),
                                                       colData = coldata.bg)
    mode(SummarizedExperiment::assay(se.bg)) = "integer"
    dds.bg = DESeq2::DESeqDataSet(se.bg, design = ~ condition)
    dds.bg$condition = stats::relevel(dds.bg$condition, reference.level)

    # running the deseq2 test
    # -----------------------
    dds.bg = DESeq2::DESeq(dds.bg, test = "Wald",
                           fitType = fitType,
                           sfType = sfType,
                           quiet = quiet,
                           minReplicatesForReplace = minReplicatesForReplace,
                           modelMatrixType = "standard",
                           useT = FALSE)


    # pulling deseq2 results
    # ----------------------
    this.contrast.bg = paste0("condition_", comp.level, "_vs_", reference.level)

    if (!is.null(filterFun)) {
        res.bg = DESeq2::results(dds.bg, name = this.contrast.bg,
                                 lfcThreshold = 0,
                                 altHypothesis = "greaterAbs",
                                 independentFiltering = independentFiltering,
                                 alpha = alpha, pAdjustMethod = pAdjustMethod,
                                 filterFun = filterFun, format = "DataFrame",
                                 minmu = minmu)
    } else {
        res.bg = DESeq2::results(dds.bg, name = this.contrast.bg,
                                 lfcThreshold = 0,
                                 altHypothesis = "greaterAbs",
                                 independentFiltering = independentFiltering,
                                 alpha = alpha, pAdjustMethod = pAdjustMethod,
                                 format = "DataFrame", minmu = minmu)
    }


    # apply fold change shrinkage if needed
    # -> requiers additional packages to be installed
    if (isTRUE(use.lfc.shrinkage)) {
        res.bg = DESeq2::lfcShrink(dds = dds.bg, res = res.bg, coef = this.contrast.bg,
                                   type = type, lfcThreshold = 0,
                                   svalue = svalue, returnList = FALSE,
                                   format = "DataFrame", saveCols = NULL,
                                   apeAdapt = apeAdapt, apeMethod = apeMethod,
                                   quiet = quiet)
    }

    # match with binding sites
    match.bs.ranges = mcols(this.ranges)[colnames(mcols(this.ranges)) == match.geneID][[1]]
    colnames(res.bg) = paste0('bg.', colnames(res.bg))
    idx = match(match.bs.ranges, rownames(res.bg))
    mcols(this.ranges) = cbind.data.frame(mcols(this.ranges), res.bg[idx,])

    object = setRanges(object, this.ranges, quiet = quiet)

    # ---
    # Store function parameters in list
    optstr = list(fitType = fitType, sfType = sfType,
                  minReplicatesForReplace = minReplicatesForReplace,
                  lfcThreshold = 0,
                  independentFiltering = independentFiltering,
                  alpha = alpha, pAdjustMethod = pAdjustMethod,
                  minmu = minmu,
                  use.lfc.shrinkage = use.lfc.shrinkage,
                  type = type,
                  svalue = svalue, apeAdapt = apeAdapt,
                  apeMethod = apeMethod, match.geneID = match.geneID)
    object@params$calculateBsFoldChange = optstr

    # ---
    # Store for results
    resultLine = data.frame(
        funName = "calculateBsFoldChange()", class = "transform",
        nIn = length(this.ranges), nOut = length(this.ranges),
        per = paste0(round(length(this.ranges)/ length(this.ranges), digits = 2)*100,"%"),
        options = paste0("Alpha=", optstr$alpha,
                         ", fitType=", optstr$fitType,
                         ", sfType=", optstr$sfType,
                         ", minReplicatesForReplace=", optstr$minReplicatesForReplace,
                         ", lfcThreshold=", optstr$lfcThreshold,
                         ", independentFiltering=", optstr$independentFiltering,
                         ", minmu=", optstr$minmu,
                         ", use.lfc.shrinkage=", optstr$use.lfc.shrinkage,
                         ifelse(length(optstr$type) == 1, paste0(", type=", optstr$type), ""),
                         ", svalue=", optstr$svalue,
                         ", apeAdapt=", optstr$apeAdapt,
                         ", apeMethod=", optstr$apeMethod,
                         ", match.geneID=", optstr$match.geneID)
    )
    object@results = rbind(object@results, resultLine)

    return(object)
}


# ------------------------------------------------------------------------------
# Exported helper
# ------------------------------------------------------------------------------


#' Impute artificial differences in the example data set
#'
#' A function that works only on the test data set provided with the package. It
#' is used for internal testing and the making of examples to showcase the
#' differential binding functions.
#'
#' Differences between samples are artificially introduced by removing the
#' signal on a random set of binding sites of the input.
#'
#' @param object a \code{\link{BSFDataSet}} object; explicitly the test data set
#' from the extdata folder
#' @param size numeric; the number of positions on which signal should be
#' deleted, counting from the start
#' @param change.per numeric; the percentage of ranges that should be effected
#' by the change.
#'
#' @return object a \code{\link{BSFDataSet}} object with the signal slot adapted
#' to reflect changes in binding between two artificial conditions
#'
#' @examples
#' # load clip data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' load(list.files(files, pattern = ".rds$", full.names = TRUE)[1])
#'
#' bds = makeBindingSites(bds, bsSize = 7)
#' bds = assignToGenes(bds, anno.genes = gns)
#'
#' bds = imputeBsDifferencesForTestdata(bds)
#'
#' @export
imputeBsDifferencesForTestdata <- function(object,
                                           size = 5,
                                           change.per = 0.1
){
    # manipulate meta data sample names
    this.meta = getMeta(object)
    this.meta$condition = factor(c("WT", "WT", "KO", "KO"), levels = c("WT", "KO"))

    # manipulate signal names
    s = getSignal(object)
    names(s$signalPlus) = paste0(this.meta$id, "_", this.meta$condition)
    names(s$signalMinus) = paste0(this.meta$id, "_", this.meta$condition)

    # manipulate ranges
    this.ranges = getRanges(object)
    set.seed(1234)
    diff.ranges = sample(this.ranges, size = length(this.ranges)* change.per)
    idx = sample(c(TRUE,FALSE), size = length(diff.ranges), replace = TRUE)
    diff.ranges.up = diff.ranges[idx]
    diff.ranges.down = diff.ranges[!idx]

    # manipulate signal
    vToKill.up = lapply(1:size, function(x){
        start(ranges(diff.ranges.up)) + x
    })
    vToKill.up = do.call(c, vToKill.up)
    s$signalPlus$`1_WT`$chr22[vToKill.up] = 0
    s$signalPlus$`2_WT`$chr22[vToKill.up] = 0
    s$signalMinus$`1_WT`$chr22[vToKill.up] = 0
    s$signalMinus$`2_WT`$chr22[vToKill.up] = 0

    vToKill.down = lapply(1:size, function(x){
        start(ranges(diff.ranges.down)) + x
    })
    vToKill.down = do.call(c, vToKill.down)
    s$signalPlus$`3_KO`$chr22[vToKill.down] = 0
    s$signalPlus$`4_KO`$chr22[vToKill.down] = 0
    s$signalMinus$`3_KO`$chr22[vToKill.down] = 0
    s$signalMinus$`4_KO`$chr22[vToKill.down] = 0

    object = setMeta(object, this.meta)
    object = setSignal(object, s)
    return(object)
}

