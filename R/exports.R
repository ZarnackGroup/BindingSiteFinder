#' Function to export sorted RBP target genes
#'
#' Genes with binding sites are target genes of the RBP. They can be exported
#' as 'csv' or 'xlsx' file. Genes can be sorted by the sum of the individual
#' binding sites score, or by the number of binding sites per gene.
#'
#' As output option, one can either output all genes in a single file, or split
#' by either gene-type or transcript-region. This options requires that either
#' \code{\link{BSFind}} or the individual functions \code{\link{assignToGenes}},
#'  and \code{\link{assignToTranscriptRegions}} were run.
#'
#' @param object a \code{BSFDataSet} object with stored ranges
#' @param path A path to where the output should be stored
#' @param format output file format
#' @param sort sorting rule for genes
#' @param split if and how the output file should be split
#'
#' @return a file of the type specified in \code{\link{format}}
#'
#' @importFrom dplyr filter group_by mutate arrange desc ungroup
#' @importFrom utils write.csv
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' \dontrun{
#' # export
#' # exportTargetGenes(bds)
#' }
#'
#' @export
exportTargetGenes <- function(object,
                              path = "./", # where to export
                              format = c("csv", "xslx"), # CSV, Excel
                              sort = c("score", "bs"), # by what to sort (score, number of bs, number of crosslinks),
                              split = c("none", "geneType", "transcriptRegion")
){
    # bind varaibles locally
    geneID <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    # type checks
    stopifnot(is(object, "BSFDataSet"))

    # capture export options
    format = match.arg(format, choices = c("csv", "xslx"))
    sort = match.arg(sort, choices = c("score", "bs"))
    split = match.arg(split, choices = c("none", "geneType", "transcriptRegion"))


    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    expRng = getRanges(object)

    # sort output
    if (sort == "bs") {
        df = getRanges(object) %>% as.data.frame(row.names = NULL) %>% filter(!is.na(geneID)) %>%
            group_by(geneID) %>% mutate(n = n()) %>% arrange(desc(n)) %>% ungroup()
    }
    if (sort == "score") {
        if (! "score" %in% colnames(mcols(getRanges(object)))) {
            stop("Sorting by score not possible, no column found with name 'score'.")
        }
        df = getRanges(object) %>% as.data.frame(row.names = NULL) %>% filter(!is.na(geneID)) %>%
            group_by(geneID) %>% mutate(n = sum(score)) %>% arrange(desc(n)) %>% ungroup()
    }

    # manage format
    if (format == "csv") {
        # manage split
        if (split == "none") {
            write.csv(x = df, file = paste0(path, "results.csv"))
        }
        if (split == "geneType") {
            if (! "geneType" %in% colnames(df)) {
                msg = paste0("Column geneType not found. Make sure to run assignToGenes() or BSFind() to use this option.")
                stop(msg)
            } else {
                sf = split(df, df$geneType)
                for (i in 1:length(sf)) {
                    cName = names(sf[i])
                    write.csv(x = sf[[i]], file = paste0(path = path, cName, ".csv"))
                }
            }
        }
        if (split == "transcriptRegion") {
            if (! "transcriptRegion" %in% colnames(df)) {
                msg = paste0("Column transcriptRegion not found. Make sure to run assignToTranscriptRegion() or BSFind() to use this option.")
                stop(msg)
            } else {
                sf = split(df, df$transcriptRegion)
                for (i in 1:length(sf)) {
                    cName = names(sf[i])
                    write.csv(x = sf[[i]], file = paste0(path = path, cName, ".csv"))
                }
            }
        }
    }

    if (format == "xslx") {
        if (!requireNamespace("xlsx", quietly=TRUE)) {
            stop("format='xlsx' requires installing the package 'xlsx'")
        }
        if (split == "none") {
            xlsx::write.xlsx(x = df, file = paste0(path, "results.xlsx"), sheetName = "All")
        }
        if (split == "geneType") {
            if (! "geneType" %in% colnames(df)) {
                msg = paste0("Column geneType not found. Make sure to run assignToGenes() or BSFind() to use this option.")
                stop(msg)
            }  else {
                sf = split(df, df$geneType)
                for (i in 1:length(sf)) {
                    cName = names(sf[i])
                    xlsx::write.xlsx(x = sf[[i]], file = paste0(path, "results.xlsx"), sheetName = cName, append = TRUE)
                }
            }
        }
        if (split == "transcriptRegion") {
            if (! "transcriptRegion" %in% colnames(df)) {
                msg = paste0("Column transcriptRegion not found. Make sure to run assignToGenes() or BSFind() to use this option.")
                stop(msg)
            }  else {
                sf = split(df, df$transcriptRegion)
                for (i in 1:length(sf)) {
                    cName = names(sf[i])
                    xlsx::write.xlsx(x = sf[[i]], file = paste0(path, "results.xlsx"), sheetName = cName, append = TRUE)
                }
            }
        }

    }

}


#' Wrapper function to export binding sites as BED files
#'
#' Function that serves as a wrapper function for \code{rtracklayer::export}.
#'
#' @param object a \code{\link{BSFDataSet}} object with stored ranges
#' @param con A path or URL
#'
#' @return a .bed file
#'
#' @import rtracklayer
#' @importFrom utils tail
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' \dontrun{
#' # export
#' # exportToBED(bds, con = "./myfile.bed")
#' }
#'
#' @export
exportToBED <- function(object, con) {
    # initialize locale variables
    funName <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))

    # range to export
    expRng = getRanges(object)

    # do a basic export if nothing was computed
    if(length(object@params) == 0){
        mcols(expRng)$name = paste0("IDX", ":", seq_along(expRng))
        mcols(expRng)$score = NULL
        rtracklayer::export(object = expRng, con = con, format = "BED",
                            trackline = new("BasicTrackLine", name = "CLS_unprocessed", description = "Unprocessed", useScore = FALSE))
    }

    # check for the highest level function that was computed and export accordingly
    lastFun = tail(object@results, 1) %>% pull(funName)
    lastFun = substr(lastFun, start = 1, stop = nchar(lastFun)-2)
    optionStr = tail(object@results, 1) %>% pull(options)

    allFunctions = c("pureClipGlobalFilter", "estimateBsWidth",
                     "pureClipGeneWiseFilter", "makeBindingSites",
                     "reproducibilityFilter", "assignToGenes",
                     "assignToTranscriptRegions", "annotateWithScore",
                     "calculateSignalToFlankScore")

    if(lastFun == allFunctions[1]) {
        mcols(expRng)$name = paste0("IDX", ":", seq_along(expRng))
        mcols(expRng)$score = NULL
        rtracklayer::export(object = expRng, con = con, format = "BED",
                            trackLine = new("BasicTrackLine", name = "CLS_pureClipGlobalfilter",
                                            description = paste0("pureClipGlobalFilter; options: ", optionStr), useScore = FALSE))
    }
    if(lastFun == allFunctions[2]) {
        mcols(expRng)$name = paste0("IDX", ":", seq_along(expRng))
        mcols(expRng)$score = NULL
        rtracklayer::export(object = expRng, con = con, format = "BED",
                            trackLine = new("BasicTrackLine", name = "CLS_estimateBsWidth",
                                            description = paste0("estimateBsWidth; options: ", optionStr), useScore = FALSE))
    }
    if(lastFun == allFunctions[3]) {
        mcols(expRng)$name = paste0("IDX", ":", seq_along(expRng), "_", expRng$geneID)
        mcols(expRng)$score = NULL
        rtracklayer::export(object = expRng, con = con, format = "BED",
                            trackLine = new("BasicTrackLine", name = "CLS_pureClipGeneWiseFilter",
                                            description = paste0("pureClipGeneWiseFilter; options: ", optionStr), useScore = FALSE))
    }
    if(lastFun == allFunctions[4]){
        mcols(expRng)$name = paste0(expRng$bsID)
        mcols(expRng)$score = NULL
        rtracklayer::export(object = expRng, con = con, format = "BED",
                            trackLine = new("BasicTrackLine", name = "BS_makeBindingsites",
                                            description = paste0("makeBindingSites; options: ", optionStr), useScore = FALSE))
    }
    if(lastFun == allFunctions[5]){
        mcols(expRng)$name = paste0(expRng$bsID)
        mcols(expRng)$score = NULL
        rtracklayer::export(object = expRng, con = con, format = "BED",
                            trackLine = new("BasicTrackLine", name = "BS_reproducibilityFilter",
                                            description = paste0("reproducibilityFilter; options: ", optionStr), useScore = FALSE))
    }
    if(lastFun == allFunctions[6]){
        mcols(expRng)$name = paste0(expRng$bsID, "_", expRng$geneID, "_", expRng$geneType)
        mcols(expRng)$score = NULL
        rtracklayer::export(object = expRng, con = con, format = "BED",
                            trackLine = new("BasicTrackLine", name = "BS_assignToGenes",
                                            description = paste0("assignToGenes; options: ", optionStr), useScore = FALSE))
    }
    if(lastFun == allFunctions[7]){
        mcols(expRng)$name = paste0(expRng$bsID, "_", expRng$geneID, "_", expRng$geneType, "_", expRng$transcriptRegion)
        mcols(expRng)$score = NULL
        rtracklayer::export(object = expRng, con = con, format = "BED",
                            trackLine = new("BasicTrackLine", name = "BS_assignToTranscriptRegions",
                                            description = paste0("assignToTranscriptRegions; options: ", optionStr), useScore = FALSE))
    }
    if(lastFun == allFunctions[8]){
        mcols(expRng)$name = paste0("BsID", ":", expRng$bsID, "_", expRng$geneID, "_", expRng$geneType, "_", expRng$transcriptRegion, "_score:", expRng$score)
        mcols(expRng)$score = NULL
        rtracklayer::export(object = expRng, con = con, format = "BED",
                            trackLine = new("BasicTrackLine", name = "BS_annotateWithScore",
                                            description = paste0("annotateWithScore; options: ", optionStr), useScore = FALSE))
    }
    if(lastFun == allFunctions[9]){
        mcols(expRng)$name = paste0("BsID", ":", expRng$bsID, "_", expRng$geneID, "_", expRng$geneType, "_", expRng$transcriptRegion, "_score:", expRng$score, "_flankScore:", expRng$signalToFlankRatio)
        mcols(expRng)$score = NULL
        rtracklayer::export(object = expRng, con = con, format = "BED",
                            trackLine = new("BasicTrackLine", name = "BS_calculateSignalToFlankScore",
                                            description = paste0("annotateWithScore; options: ", optionStr), useScore = FALSE))
    }
    lastFun = "assignToGenaes"
    if(!lastFun %in% allFunctions) {
        stop("Function names not matching. Please check your data or export manually with `rtracklayer::export(getRanges(myObj))`.\n")
    }

}


#' Create a table of all workflow steps for reporting
#'
#' Function that creates a printable table with all steps and numbers for each
#' of the workflow steps that were carried out.
#'
#' If \code{\link{option}} is set to `reduced`, only the most necessary information
#' are collected. Option `full` contains a full list of all options and parameters
#' that were set in any of the workflow functions. Option `extended` contains
#' extra information about the binding site merging step.
#'
#' @param object a \code{\link{BSFDataSet}} object with stored ranges
#' @param option character; how detailed the table should be
#'
#' @return a \code{kableExtra} table
#'
#' @importFrom kableExtra kable
#'
#'
#' @examples
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#' # apply 5% filter
#' bds = pureClipGlobalFilter(object = bds, cutoff = 0.05)
#' processingStepsTable(bds)
#'
#' @export
processingStepsTable <- function(object, # bindingSiteFinder
                                 option = c("reduced", "full", "extended") # option to which level of detail the results should be shown
) {

    nIn <- nOut <- NULL

    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))

    # handle options
    option = match.arg(option, choices = c("reduced", "full", "extended"))
    res = object@results

    res = format(res, big.mark = ",", decimal.mark = ".")
    colnames(res) = c("Step", "Type", "#N In", "#N Out", "%", "Options")
    cap = paste0("Processing step overview (", option, " view).")

    if (option == "full") {
        k = kableExtra::kable(res, caption = cap)
    }
    if (option == "reduced") {
        k = kableExtra::kable(res[c(1,3,4,5)], caption = cap)
    }
    if (option == "extended") {
        s = getSummary(object)
        ss = data.frame(step = s[[1]][2:6], nIn = s[[2]][1:5], nOut = s[[2]][2:6]) %>%
            mutate(per = paste0(round(nOut/ nIn, digits = 2)*100, "%")) %>%
            mutate(options = c(paste0("bsSize=",object@params$makeBindingSites$bsSize,
                                      ", minWidth=", object@params$makeBindingSites$minWidth,
                                      ", sub.chr=", object@params$makeBindingSites$sub.chr),
                               paste0("minCrosslinks=", object@params$makeBindingSites$minCrosslinks),
                               paste0("minClSites=", object@params$makeBindingSites$minClSites),
                               paste0("centerIsClSite=", object@params$makeBindingSites$centerIsClSite),
                               paste0("centerIsSummit=", object@params$makeBindingSites$centerIsSummit)))
        resDetail = data.frame(Step = paste0("makeBindingSites()-",ss$step),
                               Type = "binding sites", nIn = ss$nIn, nOut = ss$nOut, per = ss$per,
                               Options = ss$options)
        colnames(resDetail) = c("Step", "Type", "#N In", "#N Out", "%", "Options")

        resNewPre = res[c(1:(which(res$Step == "makeBindingSites()"))),]
        resNewPost = res[c((which(res$Step == "makeBindingSites()")+1):nrow(res)),]
        resNew = rbind.data.frame(resNewPre, resDetail, resNewPost)

        k = kableExtra::kable(resNew, caption = cap)
    }

    return(k)
}


#' Collapse signal from replicates
#'
#' Collapses all replicates merges all samples from a \link{BSFDataSet} object
#' into a single signal stream, only split by minus and plus strand.
#'
#' @param object a \code{BSFDataSet} object
#' @param collapseAll TRUE/FALSE, if all samples should be collapsed (TRUE), or
#' if samples should be kept separate by condition (FALSE)
#' @return object of type \code{\link{BSFDataSet}} with updated signal
#'
#' @seealso \code{\link{BSFDataSet}}
#'
#' @examples
#'
#' # load data
#' files <- system.file("extdata", package="BindingSiteFinder")
#' load(list.files(files, pattern = ".rda$", full.names = TRUE))
#'
#' bdsNew = collapseReplicates(bds)
#'
#' @export
collapseReplicates <- function(object, collapseAll = FALSE) {
    # stop if object is not of class BSF
    stopifnot(is(object, "BSFDataSet"))
    # stop if collapseAll parameter is not logical
    if (!is.logical(collapseAll)) {
        stop("Option collapseAll not logical, set to TRUE/FALSE")
    }

    sgn = getSignal(object)
    mta = getMeta(object)

    # collapse per condition in meta table
    if (!isTRUE(collapseAll)) {
        # get conditions to split by
        cond = levels(mta$condition)
        # handle plus strand
        plus = sgn$signalPlus
        idx = lapply(cond, function(currCond){
            grep(currCond, names(plus))
        })
        plusSum = lapply(seq_along(idx), function(x){
            currSampleIdx = idx[[x]]
            p = sgn$signalPlus[currSampleIdx]
            pSum = 0
            for (i in seq_along(p)) {
                pSum = pSum + p[[i]]
            }
            names(pSum) = names(p[[1]])
            return(pSum)
        })
        names(plusSum) = cond
        # handle minus strand
        minus = sgn$signalMinus
        idx = lapply(cond, function(currCond){
            grep(currCond, names(minus))
        })
        minusSum = lapply(seq_along(idx), function(x){
            currSampleIdx = idx[[x]]
            m = sgn$signalMinus[currSampleIdx]
            mSum = 0
            for (i in seq_along(m)) {
                mSum = mSum + m[[i]]
            }
            names(mSum) = names(m[[1]])
            return(mSum)
        })
        names(minusSum) = cond
        mrgSgn = list(signalPlus = (plusSum), signalMinus = (minusSum))
        newObj = setSignal(object, mrgSgn)
    }

    # collapse all replicates regardless of the conditions
    if (isTRUE(collapseAll)) {

        p = sgn$signalPlus
        m = sgn$signalMinus

        pSum = 0
        for (i in seq_along(p)) {
            pSum = pSum + p[[i]]
        }
        names(pSum) = names(p[[1]])
        mSum = 0
        for (i in seq_along(m)) {
            mSum = mSum + m[[i]]
        }
        names(mSum) = names(m[[1]])

        # set condition name
        lp = list(pSum)
        names(lp) = "All"
        lm = list(mSum)
        names(lm) = "All"

        mrgSgn = list(signalPlus = (lp), signalMinus = (lm))
        newObj = setSignal(object, mrgSgn)
    }

    return(newObj)
}
