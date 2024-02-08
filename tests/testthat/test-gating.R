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

data(OMIP021Samples)
wspFile <- system.file("extdata",
                       "OMIP021_samples_FlowJo.wsp",
                       package = "CytoPipelineUtils")
compMatrix <- flowCore::spillover(OMIP021Samples[[1]])$SPILL


fs_c <- CytoPipeline::runCompensation(OMIP021Samples, 
                                      spillover = compMatrix)

test_that("getFlowJoLabels works", {
    group <- "Donors"
    cellTypes <- c("Cells","Debris")
    
    labels <- getFlowJoLabels(fs_c[[1]],
                              wspFile = wspFile,
                              groups = group,
                              cellTypes = cellTypes)$labels
    
    expect_equal(sum(labels == "Cells"), 3687)
    expect_equal(sum(labels == "Debris"), 457)
    expect_equal(sum(labels == "unlabeled"), 856)
    
    labels <- getFlowJoLabels(fs_c[[2]],
                              wspFile = wspFile,
                              groups = group,
                              sampleInGroups = 2,
                              cellTypes = cellTypes)$labels
    
    expect_equal(sum(labels == "Cells"), 4383)
    expect_equal(sum(labels == "Debris"), 151)
    expect_equal(sum(labels == "unlabeled"), 466)
    
    # now with gate having parent gate and using fluorochrome markers
    # (check impact of compensation and transformation)
    # for compensation, one needs to apply exactly the same matrix as in flowJo
    # transformation should not play a role as soon as getFlowJoLabels() is
    # called on un-transformed data
    
    group <- c("Donors")
    cellTypes <- c("CD4+")
    
    labels1 <- getFlowJoLabels(fs_c[[1]],
                               wspFile = wspFile,
                               groups = group,
                               sampleInGroups = 1,
                               cellTypes = cellTypes)$labels
    
    expect_equal(sum(labels1 == "CD4+"), 1213)
    expect_equal(sum(labels1 == "unlabeled"), 3787)
    
    labels2 <- getFlowJoLabels(fs_c[[2]],
                               wspFile = wspFile,
                               groups = group,
                               sampleInGroups = 2,
                               cellTypes = cellTypes)$labels
    
    expect_equal(sum(labels2 == "CD4+"), 1339)
    expect_equal(sum(labels2 == "unlabeled"), 3661)
    
    # same with CellTypes == NULL => leaf nodes
    labels <- getFlowJoLabels(fs_c[[2]],
                              wspFile = wspFile,
                              groups = group,
                              sampleInGroups = 2)$labels
    
    expect_equal(sum(labels == "CD4+"), 1339)
    expect_equal(sum(labels == "Debris"), 151)
    expect_equal(sum(labels == "unlabeled"), 3510)
    
    # same with CellTypes == NULL and defaultCellTytpes == 'all'
    labels <- getFlowJoLabels(fs_c[[2]],
                              wspFile = wspFile,
                              defaultCellTypes = "all",
                              groups = group,
                              sampleInGroups = 2)$labels
    
    expect_equal(sum(labels == "Cells"), 3044)
    expect_equal(sum(labels == "CD4+"), 1339)
    expect_equal(sum(labels == "Debris"), 151)
    expect_equal(sum(labels == "root"), 466)
    
    # with aggregate file from samples one and two
    agg <- 
        CytoPipeline::aggregateAndSample(fs_c,
                                         nTotalEvents = 10000,
                                         seed = 1,
                                         writeOutput = FALSE)
    
    
    groups <- c("Donors","Donors")
    sampleInGroups <- c(1,2)
    cellTypes <- c("Cells", "Debris")
    
    labels <- getFlowJoLabels(agg,
                              wspFile = wspFile,
                              groups = groups,
                              cellTypes = cellTypes)$labels
    
    expect_equal(sum(labels == "Cells"), 8070)
    expect_equal(sum(labels == "Debris"), 608)
    expect_equal(sum(labels == "unlabeled"), 1322)
    
    
    # now take non disjunct gates => conflicts to solve
    groups <- c("Donors","Donors")
    sampleInGroups <- c(1,2)
    cellTypes <- c("Cells", "Debris", "CD4+")
    
    res <- getFlowJoLabels(agg,
                           wspFile = wspFile,
                           groups = groups,
                           cellTypes = cellTypes)
    labels <- res$labels
    
    expect_equal(sum(labels == "Cells"), 5518)
    expect_equal(sum(labels == "CD4+"), 2552)
    expect_equal(sum(labels == "Debris"), 608)
    expect_equal(sum(labels == "unlabeled"), 1322)
    
    expect_equal(sum(res$matrix[,"Cells"]), 8070)
    
    # what happens when cell type not found ?
    cellTypes <- c("Cells", "Debrais", "CD4+")
    
    expect_warning(resW <- getFlowJoLabels(agg,
                                 wspFile = wspFile,
                                 groups = groups,
                                 cellTypes = cellTypes),
                   regexp = "not found in any sample of any group")
    
    labels <- resW$labels
    
    expect_equal(sum(labels == "Cells"), 5518)
    expect_equal(sum(labels == "CD4+"), 2552)
    expect_equal(sum(labels == "Debrais"), 0)
    expect_equal(sum(labels == "unlabeled"), 1930)
    
    expect_equal(sum(resW$matrix[,"Cells"]), 8070)
    expect_true(all(is.na(resW$matrix[,"Debrais"])))
    
    # with subset file from sample two only (special case where file 1 is not
    # represented in aggregate)
    
    aggSample2 <- agg[flowCore::exprs(agg)[,"File"] == 2,]
    
    groups <- "Donors"
    sampleInGroups <- c(2)
    cellTypes <- c("Cells", "Debris")
    
    labels <- getFlowJoLabels(aggSample2,
                              wspFile = wspFile,
                              groups = groups,
                              sampleInGroups = sampleInGroups,
                              cellTypes = cellTypes)$labels
    
    expect_equal(sum(labels == "Cells"), 4383)
    expect_equal(sum(labels == "Debris"), 151)
    expect_equal(sum(labels == "unlabeled"), 466)
  
})

test_that("getEventNbFromFJGates works", {
    
    #1. with specified gates
    
    cellTypes <- c("Cells","Debris")
    
    eventNbDF <- getEventNbFromFJGates(
        fs_c[[1]],
        wspFile = wspFile,
        gates = cellTypes
    )
    
    expGateNames <- c("Cells", "Debris", "unlabeled")
    expGateEvents <- c(3687, 457, 856)
    
    expect_equal(eventNbDF$gate, expGateNames)
    expect_equal(eventNbDF$eventNb, expGateEvents)
    
    #2. with all gates
    
    eventNbDF <- getEventNbFromFJGates(
        fs_c[[1]],
        wspFile = wspFile,
        defaultGates = "all")
    
    expGateNames <- c("root", "Cells", "CD4+", "Debris", "unlabeled")
    expGateEvents <- c(5000, 3687, 1213, 457, 0)
    
    expect_equal(eventNbDF$gate, expGateNames)
    expect_equal(eventNbDF$eventNb, expGateEvents)
    
    #3. with leave gates only 
    
    eventNbDF <- getEventNbFromFJGates(
        fs_c[[1]],
        wspFile = wspFile)
    
    expGateNames <- c("CD4+", "Debris", "unlabeled")
    expGateEvents <- c(1213, 457, 3330)
    
    expect_equal(eventNbDF$gate, expGateNames)
    expect_equal(eventNbDF$eventNb, expGateEvents)
        
})


