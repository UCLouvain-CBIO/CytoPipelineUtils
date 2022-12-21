# CytoPipelineUtils - Copyright (C) <2022>
# <Université catholique de Louvain (UCLouvain), Belgique>
#
#   Description and complete License: see LICENSE file.
#
# This program (CytoPipelineUtils) is free software:
#   you can redistribute it and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details (<http://www.gnu.org/licenses/>).

### FUNCTIONS FOR MARGIN EVENTS REMOVAL ###

#' @title remove margin events, using flowAI
#' @description wrapper around flowAI::flow_auto_qc(), with inputs directed 
#' towards specifically remove the margin events. 
#' In the current implementation, all the signal channels, i.e. both scatter 
#' and fluo channels are scanned.
#' @param x a flowCore::flowFrame or a flowCore::flowSet
#' @param ... additional parameters passed to flowAI::flow_auto_qc(), apart
#' from the following ones : remove_from, output, ChExcludeFM, html_report,
#' mini_report, fcs_QC, fcs_highQ, fcs_lowQ, folder_results
#'
#' @return a flowCore::flowFrame with removed doublets events from the input
#' @export
#'
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' ff_m <- removeMarginsFlowAI(x = fsRaw[[2]])
#'     
removeMarginsFlowAI <- function(x, ...) {
  if (!inherits(x, "flowFrame") && !inherits(x, "flowSet")) {
    stop("x should be a flowCore::flowFrame or a flowCore::flowSet")
  }
  channel2Exclude <-
    flowCore::colnames(x)[!CytoPipeline::areSignalCols(x)]
  
  
  xOut <-
    flowAI::flow_auto_qc(x,
                         remove_from = "FM",
                         output = 1,
                         ChExcludeFM = channel2Exclude,
                         html_report = FALSE,
                         mini_report = FALSE,
                         fcs_QC = FALSE,
                         fcs_highQ = FALSE,
                         fcs_lowQ = FALSE,
                         folder_results = FALSE,
                         ...)
  return(xOut)
  
}


### FUNCTIONS FOR DOUBLETS REMOVAL ###

#' @title remove doublets from a flowFrame, using PeacoQC
#' @description wrapper around PeacoQC::RemoveDoublets().
#' Can apply the PeacoQC function subsequently on several channel pairs,
#' e.g. (FSC-A, FSC-H) and (SSC-A, SSC-H)
#' @param ff a flowCore::flowFrame
#' @param areaChannels a character vector containing the name of the 'area type'
#' channels one wants to use
#' @param heightChannels a character vector containing the name of the
#' 'height type' channels one wants to use
#' @param nmads a numeric vector with the bandwidth above the ratio allowed, per
#' channels pair (cells are kept if the ratio between -A channel\[i\] and
#' -H channel\[i\] is smaller than the median ratio + nmad\[i\] times the median
#' absolute deviation of the ratios). Default is 4, for all channel pairs.
#' @param verbose If set to TRUE, the median ratio and width will be printed.
#' @param ... additional parameters passed to PeacoQC::RemoveDoublets()
#'
#' @return a flowCore::flowFrame with removed doublets events from the input
#' @importFrom CytoPipeline appendCellID
#' @export
#'
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#'     
#' ff_c <-
#'     compensateFromMatrix(ff_m,
#'                          matrixSource = "fcs")        
#'
#' ff_s <-
#'     removeDoubletsPeacoQC(ff_c,
#'                           areaChannels = c("FSC-A", "SSC-A"),
#'                           heightChannels = c("FSC-H", "SSC-H"),
#'                           nmads = c(3, 5))                            
removeDoubletsPeacoQC <- function(ff,
                                  areaChannels,
                                  heightChannels,
                                  nmads = rep(4, length(areaChannels)),
                                  verbose = TRUE,
                                  ...) {

    # if not present already, add a column with Cell ID
    ff <- appendCellID(ff)

    # validate common scatter channel parameters
    nScatterFilters <- length(areaChannels)
    if (nScatterFilters < 1 || nScatterFilters > 2) {
        stop(
            "nb of scatter channels for doublets removal ",
            "should be either 1 or 2!"
        )
    }
    if (length(heightChannels) != nScatterFilters) {
        stop(
            "inconsistency between length of area and ",
            "height channel vectors!"
        )
    }

    if (length(nmads) != nScatterFilters) {
        stop("inconsistency between length of area channel and nMAD vectors!")
    }
    for (i in seq_len(nScatterFilters)) {
        ff <-
            PeacoQC::RemoveDoublets(ff,
                channel1 = areaChannels[i],
                channel2 = heightChannels[i],
                nmad = nmads[i],
                verbose = verbose
            )
    }

    return(ff)
}



#' @title remove doublets from a flowFrame, using flowStats
#' @description Wrapper around flowStats::singletGate().
#' Can apply the flowStats function subsequently on several channel pairs,
#' e.g. (FSC-A, FSC-H) and (SSC-A, SSC-H)
#' @param ff a flowCore::flowFrame
#' @param areaChannels a character vector containing the name of the 'area type'
#' channels one wants to use
#' @param heightChannels a character vector containing the name of the
#' 'height type' channels one wants to use
#' @param widerGate a boolean as wider_gate parameter to
#' flowStats::singletGate()
#' @param ... additional parameters passed to flowStats::singletGate()
#'
#' @return a flowCore::flowFrame with removed doublets events from the input
#' @importFrom CytoPipeline appendCellID
#' @export
#'
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#'     
#' ff_c <-
#'     compensateFromMatrix(ff_m,
#'                          matrixSource = "fcs")        
#' 
#' ff_s <-
#'     removeDoubletsFlowStats(ff_c,
#'                             areaChannels = c("FSC-A", "SSC-A"),
#'                             heightChannels = c("FSC-H", "SSC-H"),
#'                             widerGate = TRUE)
#'                             
removeDoubletsFlowStats <- function(ff,
                                    areaChannels,
                                    heightChannels,
                                    widerGate = FALSE,
                                    ...) {
    # if not present already, add a column with Cell ID
    ff <- appendCellID(ff)

    # validate common scatter channel parameters
    nScatterFilters <- length(areaChannels)
    if (nScatterFilters < 1 || nScatterFilters > 2) {
        stop(
            "nb of scatter channels for doublets removal ",
            "should be either 1 or 2!"
        )
    }
    if (length(heightChannels) != nScatterFilters) {
        stop(
            "inconsistency between length of area ",
            "and height channel vectors!"
        )
    }

    for (i in seq_len(nScatterFilters)) {
        currentSingletGate <-
            flowStats::singletGate(ff,
                filterId = paste0(
                    "Singlets_",
                    areaChannels[i]
                ),
                area = areaChannels[i],
                height = heightChannels[i],
                wider_gate = widerGate
            )

        if (i == 1) {
            singletGateCombined <- currentSingletGate
        } else {
            singletGateCombined <- singletGateCombined & currentSingletGate
        }
    }

    fltSinglet <- flowCore::filter(ff, singletGateCombined)

    ff <- flowCore::Subset(ff, fltSinglet)
    #ff <- ff[fltSinglet@subSet, ]

    return(ff)
}

#' @title remove debris from a flowFrame, using flowClust
#' @description this function removes debris from a flowFrame,
#' using clustering capabilities of flowClust::tmixFilter(). The idea is to
#' pre-select a number of clusters to be found in the (FSC,SSC) 2D view, and
#' eliminate the cluster that is the closest to the origin.

#' @param ff a flowCore::flowFrame
#' @param FSCChannel the name of the FSC channel
#' @param SSCChannel the name of the SSC channel
#' @param nClust number of clusters to identify
#' @param ... additional parameters passed to flowClust::tmixFilter()
#'
#' @return a flowCore::flowFrame with removed debris events from the input
#' @importFrom CytoPipeline appendCellID
#' @export
#' 
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#'     
#' ff_c <-
#'     compensateFromMatrix(ff_m,
#'                          matrixSource = "fcs")        
#' 
#' 
#' ff_cells <-
#'     removeDebrisFlowClustTmix(ff_c,
#'                               FSCChannel = "FSC-A",
#'                               SSCChannel = "SSC-A",
#'                               nClust = 3,
#'                               level = 0.97,
#'                               B = 100)
#' 
removeDebrisFlowClustTmix <- function(ff,
                                      FSCChannel,
                                      SSCChannel,
                                      nClust,
                                      ...) {

    # if not present already, add a column with Cell ID
    ff <- appendCellID(ff)

    # handle ellipsis arguments, as 'tmixFilter' does not accept unknown args
    passedEllipsisArgs <- list(...)
    newEllipsisArgs <- list()

    argNames <-
        c(
            "expName", "K", "B", "tol", "nu", "lambda", "nu.est", "trans",
            "min.count", "max.count", "min", "max", "level", "u.cutoff",
            "z.cutoff", "randomStart", "B.init", "tol.init", "seed", "criterion"
        )
    for (argN in argNames) {
        if (!is.null(passedEllipsisArgs[[argN]])) {
            newEllipsisArgs[[argN]] <- passedEllipsisArgs[[argN]]
        }
    }

    cellsFilter <-
        do.call(flowClust::tmixFilter,
            args = c(
                list(
                    filterId = "tmixFilter",
                    parameters = c(FSCChannel, SSCChannel),
                    K = nClust
                ),
                newEllipsisArgs
            )
        )


    resCellsFilter <- flowCore::filter(ff, cellsFilter)

    FSCMedians <- vapply(
        X = seq_len(nClust),
        FUN.VALUE = double(1),
        FUN = function(x, ff, flt) {
            resCellsFltr <- flt[[x]]
            # stats::median(flowCore::exprs(ff)[
            #     resCellsFltr@subSet, FSCChannel
            # ],
            # na.rm = TRUE
            # )
            stats::median(
                flowCore::exprs(
                    flowCore::Subset(ff, resCellsFltr))[, FSCChannel],
                na.rm = TRUE)
        },
        ff = ff, flt = resCellsFilter
    )

    debrisIndex <- which.min(FSCMedians)
    keptClustersIndexes <- setdiff(seq_len(nClust), debrisIndex)
    tokeepFilter <- resCellsFilter[[keptClustersIndexes[1]]]
    if (nClust > 2) {
        for (i in keptClustersIndexes[-1]) {
            tokeepFilter <- tokeepFilter | resCellsFilter[[i]]
        }
    }
    selectedCells <- flowCore::filter(ff, tokeepFilter)
    ff <- flowCore::Subset(ff, selectedCells)
    #ff <- ff[selectedCells@subSet, ]

    return(ff)
}

### FUNCTIONS FOR DEAD CELLS REMOVAL ###


# @title remove dead cells from a flowFrame
# @description this function removes dead cells from a flowFrame, using a
# specific '(a)live/dead' channel, and the openCyto::gate_tail() gating
# function (see doc of the openCyto package)

# @param ff a flowCore::flowFrame
# @param preTransform if TRUE, apply the transList scale transform prior to
# running the gating algorithm
# @param transList applied in conjunction with preTransform == TRUE
# @param LDMarker a character containing the exact name of the marker
# corresponding to Live/Dead channel, or the Live/Dead channel name itself
# @param ... additional parameters passed to openCyto::gate_tail()
#
# @return a flowCore::flowFrame with removed dead cells from the input
# @export
# 
# @importFrom CytoPipeline appendCellID
#
# @examples
#
# rawDataDir <-
#     paste0(system.file("extdata", package = "CytoPipeline"), "/")
# sampleFiles <-
#     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
# 
# truncateMaxRange <- FALSE
# minLimit <- NULL
# 
# # create flowCore::flowSet with all samples of a dataset
# fsRaw <- readSampleFiles(
#     sampleFiles = sampleFiles,
#     whichSamples = "all",
#     truncate_max_range = truncateMaxRange,
#     min.limit = minLimit)
# 
# suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#     
# ff_c <-
#     compensateFromMatrix(ff_m,
#                          matrixSource = "fcs")        
#
# transList <- 
#     estimateScaleTransforms(        
#         ff = ff_c,
#         fluoMethod = "estimateLogicle",
#         scatterMethod = "linear",
#         scatterRefMarker = "BV785 - CD3")
# 
# ff_lcells <-
#     removeDeadCellsGateTail(ff_c,
#                             preTransform = TRUE,
#                             transList = transList,
#                             LDMarker = "L/D Aqua - Viability",
#                             num_peaks = 2,
#                             ref_peak = 2,
#                             strict = FALSE,
#                             positive = FALSE)
#                             
# removeDeadCellsGateTail <- function(ff,
#                                     preTransform = FALSE,
#                                     transList = NULL,
#                                     LDMarker,
#                                     ...) {
#     
#     # handle ellipsis arguments, as 'openCyto::gate_tail'
#     # does not accept unknown args
#     passedEllipsisArgs <- list(...)
#     
#     newArgs <- list()
#     
#     argNames <-
#         c(
#             "num_peaks", "ref_peak", "strict", "tol", "side", "min", "max",
#             "bias", "positive", "deriv", "bandwidth", "adjust", "num_points",
#             "range.x", "binned", "se", "w"
#         )
#     for (argN in argNames) {
#         if (!is.null(passedEllipsisArgs[[argN]])) {
#             newArgs[[argN]] <- passedEllipsisArgs[[argN]]
#         }
#     }
#     
#     # will be needed by openCyto::gate_tail to support 
#     # BiocParallel::SnowParams
#     #requireNamespace("flowCore")
#     #require(flowCore, include.only = c("exprs", "rectangleGate"))
#     
#     message(paste0("Removing Dead Cells Gate Tail events from file : ", 
#                    flowCore::identifier(ff)))
#     
#     # if not present already, add a column with Cell ID
#     ff <- appendCellID(ff)
#     
#     if (preTransform) {
#         if (is.null(transList)) {
#             stop(
#                 "tranformation list needs to be provided ",
#                 "if preTransform = TRUE!"
#             )
#         }
#         ffIn <- flowCore::transform(ff, transList)
#     } else {
#         ffIn <- ff
#     }
#     
#     if (LDMarker %in% flowCore::colnames(ff)) {
#         LDChannel <- LDMarker
#     } else {
#         LDChannel <- getChannelNamesFromMarkers(ffIn, markers = LDMarker)
#     }
#     
#     
#     
#     liveGate <-
#         do.call(openCyto::gate_tail,
#                 args = c(
#                     list(ffIn,
#                          channel = LDChannel,
#                          filterId = "Live_Cells"
#                     ),
#                     newArgs
#                 )
#         )
#     
#     selectedLive <- flowCore::filter(ffIn, liveGate)
#     
#     # note we take ff and not ffIn (no transfo)
#     ff <- flowCore::Subset(ff, selectedLive)
#     
#     return(ff)
#         
#     
#     
#     
#     
# 
# }

#' @title remove dead cells from a flowFrame
#' @description this function removes dead cells from a flowFrame, using a
#' specific '(a)live/dead' channel, and the flowDensity::deGate() gating
#' function (see doc of the flowDensity package)
#'
#' @param ff a flowCore::flowFrame
#' @param preTransform if TRUE, apply the transList scale transform prior to
#' running the gating algorithm
#' @param transList applied in conjunction with preTransform == TRUE
#' @param LDMarker a character containing the exact name of the marker
#' corresponding to Live/Dead channel, or the Live/Dead channel name itself
#' @param keepPositivePop logical flag stating whether we want to keep, 
#' after gating, the population that is positive for `LDMarker` (TRUE) 
#' or negative for `LDmarker` (FALSE)  
#' @param ... additional parameters passed to flowDensity::deGate()
#'
#' @return a flowCore::flowFrame with removed dead cells from the input
#' @export
#' 
#' @importFrom CytoPipeline appendCellID
#'
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#'     
#' ff_c <-
#'     compensateFromMatrix(ff_m,
#'                          matrixSource = "fcs")        
#'
#' transList <- 
#'     estimateScaleTransforms(        
#'         ff = ff_c,
#'         fluoMethod = "estimateLogicle",
#'         scatterMethod = "linear",
#'         scatterRefMarker = "BV785 - CD3")
#' 
#' ff_lcells <-
#'     removeDeadCellsDeGate(ff_c,
#'                           preTransform = TRUE,
#'                           transList = transList,
#'                           LDMarker = "L/D Aqua - Viability",
#'                           keepPositivePop = FALSE)
#'                             
removeDeadCellsDeGate <- function(ff,
                                  preTransform = FALSE,
                                  transList = NULL,
                                  LDMarker, 
                                  keepPositivePop = FALSE,
                                  ...) {

    
    message("Removing Dead Cells events (using flowDensity::deGate())", 
            " from file : ", flowCore::identifier(ff))

    # if not present already, add a column with Cell ID
    ff <- appendCellID(ff)

    if (preTransform) {
        if (is.null(transList)) {
            stop(
                "tranformation list needs to be provided ",
                "if preTransform = TRUE!"
            )
        }
        ffIn <- flowCore::transform(ff, transList)
    } else {
        ffIn <- ff
    }

    if (LDMarker %in% flowCore::colnames(ff)) {
        LDChannel <- LDMarker
    } else {
        LDChannel <- 
          CytoPipeline::getChannelNamesFromMarkers(ffIn, markers = LDMarker)
    }
    
    threshold <- 
        flowDensity::deGate(ffIn,
                            channel = LDChannel,
                            ...)
    
    if (keepPositivePop) {
        selectedLive <- flowCore::exprs(ffIn)[, LDChannel] >= threshold    
    } else {
        selectedLive <- flowCore::exprs(ffIn)[, LDChannel] <= threshold
    }



    # note we take ff and not ffIn (no transfo)
    ff <- ff[selectedLive, ]

    return(ff)

}

### FUNCTIONS for Quality Control ###


#' @title perform QC with flowCut
#' @description this function is a wrapper around flowCut::flowCut()
#' function.
#' It also pre-selects the channels to be handled (=> all signal channels)
#' @param ff a flowCore::flowFrame
#' @param preTransform if TRUE, apply the transList scale transform prior to
#' running the gating algorithm
#' @param transList applied in conjunction with preTransform
#' @param outputDiagnostic if TRUE, stores diagnostic files generated by
#' flowCut in outputDir directory
#' @param outputDir used in conjunction with outputDiagnostic
#' @param verbose if TRUE messages comments on the QC process
#' @param ... additional parameters passed to flowCut::flowCut()
#'
#' @return a flowCore::flowFrame with removed low quality events from the input
#' @importFrom CytoPipeline appendCellID
#' @export
#'
#'
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#'
#' ff_QualityControl <-
#'     qualityControlFlowCut(
#'         fsRaw[[2]],
#'         MaxContin = 0.1,
#'         MeanOfMeans = 0.13,
#'         MaxOfMeans = 0.15,
#'         MaxValleyHgt = 0.1,
#'         MaxPercCut = 0.3,
#'         LowDensityRemoval = 0.1,
#'         RemoveMultiSD = 7,
#'         AlwaysClean = FALSE,
#'         IgnoreMonotonic = FALSE,
#'         MonotonicFix = NULL,
#'         Measures = c(1:8))
#'       
qualityControlFlowCut <- function(ff,
                                  preTransform = FALSE,
                                  transList = NULL,
                                  outputDiagnostic = FALSE,
                                  outputDir = NULL,
                                  verbose = TRUE,
                                  ...) {
    # if not present already, add a column with Cell ID
    ff <- appendCellID(ff)

    if (preTransform) {
        if (is.null(transList)) {
            stop(
                "tranformation list needs to be provided ",
                "if preTransform = TRUE!"
            )
        }
        ffIn <- flowCore::transform(ff, transList)
    } else {
        ffIn <- ff
    }

    # qualityControl with flowCut
    message("Applying flowCut method...")
    channelsIndices <- which(CytoPipeline::areSignalCols(ffIn))
    if (outputDiagnostic) {
        Plot <- "All"
        if (!is.null(outputDir)) {
            Directory <- outputDir
        } else {
            filePrefixWithDir <- "resultsQC"
        }
    } else {
        Directory <- NULL # not used
        Plot <- "None"
    }

    res <-
        flowCut::flowCut(
            f = ffIn,
            Channels = channelsIndices,
            Directory = Directory,
            FileID = NULL,
            Plot = Plot,
            Verbose = verbose,
            ...
        )
    # browser()
    badEventIDs <- res$ind

    goodEvents <- !(seq_len(flowCore::nrow(ffIn)) %in% badEventIDs)
    ff <- ff[goodEvents, ] # note we take ff and not ffIn (no transfo)

    return(ff)
}

#' @title perform QC with flowClean
#' @description this function is a wrapper around flowClean::clean()
#' function.
#' It also pre-selects the channels to be handled (=> all signal channels)
#' @param ff a flowCore::flowFrame
#' @param preTransform if TRUE, apply the transList scale transform prior to
#' running the gating algorithm
#' @param transList applied in conjunction with preTransform
#' @param outputDiagnostic if TRUE, stores diagnostic files generated by
#' flowClean in outputDir directory
#' @param outputDir used in conjunction with outputDiagnostic
#' @param verbose if TRUE messages comments on the QC process
#' @param ... additional parameters passed to flowClean::clean()
#' 
#' @return a flowCore::flowFrame with removed low quality events from the input
#' @importFrom CytoPipeline appendCellID
#' 
#' @export
#' 
#' @examples
#' 
#' rawDataDir <-
#' paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#' paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#' sampleFiles = sampleFiles,
#' whichSamples = "all",
#' truncate_max_range = truncateMaxRange,
#' min.limit = minLimit)
#' 
#' ff_QualityControl <- suppressWarnings(
#' qualityControlFlowClean(fsRaw[[2]],
#' binSize = 0.01, # default
#' nCellCutoff = 500, # default
#' cutoff = "median", # default
#' fcMax = 1.3, # default
#' nstable = 5))
#' 
qualityControlFlowClean <- function(ff,
                                    preTransform = FALSE,
                                    transList = NULL,
                                    outputDiagnostic = FALSE,
                                    outputDir = NULL,
                                    verbose = TRUE,
                                    ...) {
    #browser()
    # if not present already, add a column with Cell ID
    ff <- appendCellID(ff)

    if (preTransform) {
        if (is.null(transList)) {
            stop(
                "tranformation list needs to be provided ",
                "if preTransform = TRUE!"
            )
        }
        ffIn <- flowCore::transform(ff, transList)
    } else {
        ffIn <- ff
    }

    # qualityControl with flowClean
    message("Applying flowClean method...")
    vectMarkers <- which(CytoPipeline::areSignalCols(ffIn))

    if (outputDiagnostic) {
        diagnostic <- TRUE
        if (!is.null(outputDir)) {
            filePrefixWithDir <- outputDir
        } else {
            filePrefixWithDir <- "resultsQC"
        }

        # add original fcs file name in prefix,
        # as flowClean is designed to work for one fcs at the time
        filename <- basename(flowCore::keyword(ffIn, "FILENAME")$FILENAME)
        # removing extension
        filename <- sub("([^.]+)\\.[[:alnum:]]+$", "\\1", filename)
        filePrefixWithDir <- paste0(filePrefixWithDir, filename)
    } else {
        filePrefixWithDir <- NULL # not used
        diagnostic <- FALSE
    }

    # note we call it using returnVector = FALSE and extract
    # the goodVsBadVector later to work around a flowClean bug (to be corrected)
    tempDF <-
        flowClean::clean(
            fF = ffIn,
            vectMarkers = vectMarkers,
            filePrefixWithDir = filePrefixWithDir,
            ext = ".fcs", # not used
            diagnostic = diagnostic,
            announce = verbose,
            returnVector = FALSE,
            ...
        )

    areGoodEvents <- tempDF[, "GoodVsBad"] < 10000
    ff <- ffIn[areGoodEvents, ]

    return(ff)
}

