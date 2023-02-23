

#' @title Process a flowjo workspace file
#' @description Reads a flowjo workspace file using the flowWorkspace library
#' and returns a vector with a label for each cell of a flowFrame
#' from a set of specified gates. If the flowFrame was generated from several
#' files (contains a 'File' column containing a file index), then the groups
#' should contain the group name to apply for each of these files/samples
#' @param ff a flowCore::flowFrame
#' @param wspFile a flowjo workspace
#' @param groups vector of flowjo groups to parse from the flowjo workspace
#' (i.e. one for each sample/file used to generate the flowFrame).
#' Default c("All Samples")
#' @param sampleInGroups vector of flowjo sample index in group to use from
#' the flowjo workspace (i.e. one for each sample/file used to generate the
#' flowFrame). If not provided, will start from 1 in each group and increment.
#' Default NULL.
#' @param cellTypes vector of flowJo gate nodes to parse (should be common to
#' all groups!). Default NULL then uses either all leaf node from first group,
#' or all gates from first group, depending on 'defaultCellTypes'.
#' @param defaultCellTypes if cellTypes are not specified, than triggers the
#' choice of either :
#' - all leaf nodes ('leaves') (the default)
#' - all nodes ('all')
#' @param specialCharsInChannels vector of special char litterals that are
#' replaced in FlowJo workspace
#' @param FlowJoCharsinChannels vector of new char litterals to be applied
#' in FlowJo workspace
#' @param withFJv10TimeCorrection if TRUE (default), applies a time correction
#' to cope with a bug in flowJo 10.0. Should be kept TRUE with FloJo v>10.0 if
#' time channel is used in gating.
#' 
#' Since FlowJo v10, time channel is read differently from fcs files (see 
#' https://docs.flowjo.com/flowjo/experiment-based-platforms/kinetics/
#' how-flowjo-handles-time/). The resulting scale transformation on the time 
#' channel is normally stored in the corresponding tranformation gain parameter 
#' in the FlowJo workspace.
#' However, somehow this gain parameter can be sometimes wrong in the .wsp file, 
#' and so the workaround consists in correcting the time channel values of the
#' cytoset prior to gating, in order to have the gate labels calculated 
#' correctly.
#' @param ... Extra arguments passed to `CytoML::flowjo_to_gatingset()`
#'
#' @return a list with :
#' - the first element ("matrix") is a matrix containing filtering results for
#' each specified gate
#' - the second element ("labels") is a vector which assigns one label to each
#' cell. If no cell type correspond to a specific cell, than the keyword
#' 'unlabeled' is assigned to this cell. If the cell belongs to several gates
#' (meaning that the gates ar not disjoints), than this cell is assigned to the
#' gate with the less matching cells.
#' @export
#'
getFlowJoLabels <- function(ff,
                            wspFile,
                            groups = "All Samples",
                            sampleInGroups = NULL,
                            cellTypes = NULL,
                            defaultCellTypes = c("leaves", "all"),
                            specialCharsInChannels = c("/"),
                            FlowJoCharsinChannels = c("_"),
                            withFJv10TimeCorrection = TRUE,
                            ...)
{
  
  if (!inherits(ff, "flowFrame"))
    stop("ff should be a flowframe")
  
  if (length(groups) == 0 || !is.character(groups)) {
    stop("groups should contain at least one group and groups should be of ",
         " character type!")
  }
  
  # replace special characters by flowJo ones in channel names
  if (length(specialCharsInChannels) != length(FlowJoCharsinChannels)) {
    msg <- "inconsistency in lengthes of 'specialCharsInChannels' vs"
    msg <- paste0(msg, "'FlowJoCharsinChannels'")
    stop(msg)
  }
  
  cnames <- flowCore::colnames(ff)
  for (i in seq_along(specialCharsInChannels)) {
    for (j in seq_along(cnames)) {
      cnames[j] <- gsub(specialCharsInChannels[i],
                        FlowJoCharsinChannels[i],
                        cnames[j])
    }
  }
  flowCore::colnames(ff) <- cnames
  
  # check consistence between groups and number of original files in ff
  # and store array of cells per file
  whichCellsPerSample <- list()
  
  if ("File" %in% flowCore::colnames(ff)) {
    sampleIdsInFile <- unique(flowCore::exprs(ff)[,"File"])
    nSamples <- length(sampleIdsInFile)
    if (length(groups) != nSamples)
      stop("groups should contain a group per sample ",
           "from which the flow frame was generated")
    for(i in 1:nSamples){
      whichCellsPerSample[[i]] <-
        which(flowCore::exprs(ff)[,"File"] == sampleIdsInFile[i])
    }
  } else {
    nSamples <- 1
    if (length(groups) != nSamples)
      stop("only one sample was detected in flow frame ",
           "=> groups should contain only one element")
    whichCellsPerSample[[1]] <- seq(nrow(flowCore::exprs(ff)))
  }
  
  ws <- CytoML::open_flowjo_xml(wspFile)
  
  # store gating set per distinct group
  
  #browser()
  distinctGroups <- unique(groups)
  GSPerGroup <- list()
  
  timeCh <- CytoPipeline::findTimeChannel(ff)
  timeStep <- 0.01
  correctTimeGain <- 0.01
  
  if (!is.null(timeCh)) {
    btim_kwd <- flowCore::keyword(ff, "$BTIM")
    etim_kwd <- flowCore::keyword(ff, "$ETIM")
    
    # work-around to cope with different handling of time channel
    # since FlowJo v10 (not yet integrated into flowWorkspace)
    if (withFJv10TimeCorrection) {
      if (!is.null(flowCore::keyword(ff, "$TIMESTEP"))) {
        timeStep <- as.numeric(flowCore::keyword(ff, "$TIMESTEP"))
      } else {
        message("$TIMESTEP keywork not found, applying 0.01 by default")
      }
      if (!is.null(btim_kwd) && !is.null(etim_kwd)) {
        btim <- as.POSIXct(flowCore::keyword(ff, "$BTIM")[["$BTIM"]],
                           format="%H:%M:%S")
        etim <- as.POSIXct(flowCore::keyword(ff, "$ETIM")[["$ETIM"]],
                           format="%H:%M:%S")
        
        targetDiffTime <- as.numeric(etim - btim, units = "secs")
        timeOrder <- order(flowCore::exprs(ff)[,timeCh])
        actualDiffTime <- 
          (flowCore::exprs(ff)[timeOrder[flowCore::nrow(ff)], timeCh] - 
             flowCore::exprs(ff)[timeOrder[1], timeCh]) 
        correctTimeGain <- targetDiffTime / actualDiffTime
      } else {
        correctTimeGain <- timeStep
        message("$BTIM and/or $ETIM keyword not found, unable to calculate",
                "correct gain in FJ10 workaround => applying timestep value")
      }
    }
  }
  
  timeChannelCFs <- list()
  for (gr in distinctGroups) {
    # extend_val set to -Inf to discard any gate automatic extension
    # for all vals <0 to -4000
    GSPerGroup[[gr]] <- CytoML::flowjo_to_gatingset(ws,
                                                    name = gr,
                                                    execute = FALSE,
                                                    transform = TRUE,
                                                    extend_val = -Inf,
                                                    ...)
    
    if (!is.null(timeCh) && withFJv10TimeCorrection) {
      # getting actual time gain parameter retrieved from wsp file
      transfoFunction <- 
        flowWorkspace::gh_get_transformations(GSPerGroup[[gr]][[1]],
                                              channel = timeCh)
      if (is.null(transfoFunction)) {
        actualTimeGain <- 1.0
        message("No transformation found for time channel, assuming actual ",
                "gain in FJ10 workaround is 1")
      } else {
        actualTimeGain <- transfoFunction(1.0)
      }
      
      timeChannelCFs[[gr]] <- correctTimeGain / actualTimeGain
    }
  }
  
  # if cellTypes not provided => take all leaf node of first gating set
  
  if (is.null(cellTypes)) {
    defaultCellTypes <- match.arg(defaultCellTypes)
    gr <- groups[1]
    if (defaultCellTypes == "leaves") {
      cellTypes <- flowWorkspace::gs_get_leaf_nodes(GSPerGroup[[gr]],
                                                    path = "auto")
    } else {
      cellTypes <- flowWorkspace::gs_get_pop_paths(GSPerGroup[[gr]],
                                                   path = "auto")
    }
  }
  
  # per cellType, store array of corresponding cell
  cellMap <- list()
  foundCellType <- rep(FALSE, length(cellTypes))
  
  for (i in 1:nSamples) {
    gr <- groups[i]
    if (is.null(sampleInGroups)) {
      whichSampleInGroup <- sum((groups == groups[i])[1:i])
    } else {
      whichSampleInGroup <- sampleInGroups[i]
    }
    
    #browser()
    
    cs <- flowWorkspace::cytoset()
    flowWorkspace::cs_add_cytoframe(
      cs, "unique_sample",
      flowWorkspace::flowFrame_to_cytoframe(ff[whichCellsPerSample[[i]]]))
    
    # applying work-around for time channel if needed
    if (!is.null(timeCh) && withFJv10TimeCorrection) {
      flowCore::exprs(cs[[1]])[,timeCh] <- 
        flowCore::exprs(cs[[1]])[,timeCh] * timeChannelCFs[[gr]]
    }
    
    gs <-
      suppressMessages(
        flowWorkspace::gh_apply_to_cs(
          x = GSPerGroup[[gr]][[whichSampleInGroup]],
          cs = cs,
          compensation_source = "none"))
    
    nodeListFull <- flowWorkspace::gs_get_pop_paths(gs, path = "full")
    nodeList <- flowWorkspace::gs_get_pop_paths(gs, path = "auto")
    
    for (j in seq_along(cellTypes)) {
      ct <- cellTypes[j]
      nodeIndex <- which(nodeList == ct)
      if (length(nodeIndex) == 0) {
        next
      } else if (length(nodeIndex) > 1) {
        msg <- paste0("Ambiguity with cell type [", ct, "], found several ")
        msg <- paste0(msg, "times in sample ", i)
        warning(msg)
        nodeIndex <- nodeIndex[1]
      }
      
      foundCellType[j] <- TRUE
      node <- nodeListFull[nodeIndex]
      
      cellIndices <- flowWorkspace::gh_pop_get_indices(gs[[1]],
                                                       y = node)
      cellMap[[ct]] <- c(cellMap[[ct]], whichCellsPerSample[[i]][cellIndices])
    }
  }
  
  for (j in seq_along(cellTypes)) {
    if (!foundCellType[j]) {
      ct <- cellTypes[j]
      stop("Cell type [", ct, "] not found in any sample of any group !")
    }
  }
  
  # prepare outputs
  # 1. 'matrix with gates in columns and cells in rows, data = TRUE/FALSE'
  res <- list()
  nEvents <- flowCore::nrow(ff)
  nCellTypes <- length(cellTypes)
  res$matrix <- matrix(data = rep(FALSE, nEvents*(nCellTypes+1)),
                       nrow = nEvents,
                       ncol = nCellTypes+1,
                       dimnames = list(NULL ,c(cellTypes, "unlabeled")))
  
  for (ct in cellTypes){
    res$matrix[cellMap[[ct]], ct] <- TRUE
  }
  
  res$matrix[, nCellTypes+1] <- (rowSums(res$matrix) == 0)
  
  
  # 2. 'labels' with one label attached to each cell. If conflict, take the
  # label of the most rare population
  res$labels <- rep("unlabeled", nEvents)
  if (nCellTypes > 1) {
    nCellsByType <- colSums(res$matrix[,1:nCellTypes])
  } else {
    nCellsByType <- sum(res$matrix[,1])
  }
  
  for (ct in cellTypes[order(nCellsByType, decreasing = TRUE)]) {
    res$labels[cellMap[[ct]]] <- ct
  }
  
  return(res)
}