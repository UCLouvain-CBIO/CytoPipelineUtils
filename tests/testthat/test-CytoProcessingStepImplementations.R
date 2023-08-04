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

# obtain OMIP021UTSamples, light-weight version used specifically for these 
# unit tests
path <- system.file("scripts",
                    package = "CytoPipeline"
)

source(file.path(path,"MakeOMIP021UTSamples.R"))

test_that("removeMarginsflowAI works", {
    fs_raw <- OMIP021UTSamples

    ff_m <-
        suppressWarnings(removeMarginsFlowAI(x = fs_raw[[1]]))

    ref_ff_m <- readRDS(test_path("fixtures", "ff_m2.rds"))

    #saveRDS(ff_m, test_path("fixtures", "ff_m2.rds"))

    expect_equal(
        flowCore::exprs(ff_m),
        flowCore::exprs(ref_ff_m)
    )
})



test_that("removeDoubletsFlowStats works", {
    ref_ff_c <- readRDS(test_path("fixtures", "ff_c.rds"))


    ff_s <-
        removeDoubletsFlowStats(ref_ff_c,
            areaChannels = c("FSC-A", "SSC-A"),
            heightChannels = c("FSC-H", "SSC-H"),
            widerGate = TRUE
        )

    ref_ff_s <- readRDS(test_path("fixtures", "ff_s.rds"))

    # saveRDS(ff_s, test_path("fixtures", "ff_s.rds"))

    expect_equal(
        flowCore::exprs(ff_s),
        flowCore::exprs(ref_ff_s)
    )
})


test_that("removeDoubletsPeacoQC works", {
    ref_ff_c <- readRDS(test_path("fixtures", "ff_c.rds"))

    ff_s2 <-
        removeDoubletsPeacoQC(ref_ff_c,
            areaChannels = c("FSC-A", "SSC-A"),
            heightChannels = c("FSC-H", "SSC-H"),
            nmads = c(3, 5)
        )

    ref_ff_s2 <- readRDS(test_path("fixtures", "ff_s2.rds"))

    # saveRDS(ff_s2, test_path("fixtures", "ff_s2.rds"))

    expect_equal(
        flowCore::exprs(ff_s2),
        flowCore::exprs(ref_ff_s2)
    )
})

test_that("removeDebrisFlowClustTmix works", {
    ref_ff_s <- readRDS(test_path("fixtures", "ff_s.rds"))

    ff_cells <-
        removeDebrisFlowClustTmix(ref_ff_s,
            FSCChannel = "FSC-A",
            SSCChannel = "SSC-A",
            nClust = 3,
            level = 0.97,
            B = 100
        )

    ref_ff_cells <- readRDS(test_path("fixtures", "ff_cells.rds"))

    # saveRDS(ff_cells, test_path("fixtures", "ff_cells.rds"))

    expect_equal(
        flowCore::exprs(ff_cells),
        flowCore::exprs(ref_ff_cells)
    )
})


# test_that("removeDeadCellsGateTail works", {
#     ref_ff_cells <- readRDS(test_path("fixtures", "ff_cells.rds"))
# 
#     transListPath <- file.path(system.file("extdata", 
#                                            package = "CytoPipeline"),
#                             "/OMIP021_TransList.rds")
#     refTransList <- readRDS(transListPath)
# 
#     ff_lcells <-
#         removeDeadCellsGateTail(ref_ff_cells,
#             preTransform = TRUE,
#             transList = refTransList,
#             LDMarker = "L/D Aqua - Viability",
#             num_peaks = 2,
#             ref_peak = 2,
#             strict = FALSE,
#             positive = FALSE
#         )
# 
#     ref_ff_lcells <- readRDS(test_path("fixtures", "ff_lcells.rds"))
# 
#     # saveRDS(ff_lcells, test_path("fixtures", "ff_lcells.rds"))
# 
#     expect_equal(
#         flowCore::exprs(ff_lcells),
#         flowCore::exprs(ref_ff_lcells)
#     )
#     
#     # same with channel name instead of marker for L/D
#     ff_lcells2 <-
#         removeDeadCellsGateTail(ref_ff_cells,
#                                 preTransform = TRUE,
#                                 transList = refTransList,
#                                 LDMarker = "Comp-525/50Violet-A",
#                                 num_peaks = 2,
#                                 ref_peak = 2,
#                                 strict = FALSE,
#                                 positive = FALSE
#         )
#     
#     expect_equal(
#         flowCore::exprs(ff_lcells2),
#         flowCore::exprs(ref_ff_lcells)
#     )
# })

test_that("removeDeadCellsDeGate works", {
    ref_ff_cells <- readRDS(test_path("fixtures", "ff_cells.rds"))
    
    transListPath <- file.path(system.file("extdata", 
                                           package = "CytoPipeline"),
                            "OMIP021_TransList.rds")
    refTransList <- readRDS(transListPath)
    
    ff_lcells <-
        removeDeadCellsDeGate(ref_ff_cells,
                              preTransform = TRUE,
                              transList = refTransList,
                              LDMarker = "L/D Aqua - Viability")
    
    ref_ff_lcells <- readRDS(test_path("fixtures", "ff_lcells.rds"))
    
    # saveRDS(ff_lcells, test_path("fixtures", "ff_lcells.rds"))
    
    expect_equal(
        flowCore::exprs(ff_lcells),
        flowCore::exprs(ref_ff_lcells)
    )
    
    # same with channel name instead of marker for L/D
    ff_lcells2 <-
        removeDeadCellsDeGate(ref_ff_cells,
                              preTransform = TRUE,
                              transList = refTransList,
                              LDMarker = "Comp-525/50Violet-A"
        )
    
    expect_equal(
        flowCore::exprs(ff_lcells2),
        flowCore::exprs(ref_ff_lcells)
    )
})




test_that("qualityControlFlowCut works", {
    fs_raw <- OMIP021UTSamples

    ff_QualityControl <- suppressMessages(
        qualityControlFlowCut(fs_raw[[1]],
            MaxContin = 0.1,
            MeanOfMeans = 0.13,
            MaxOfMeans = 0.15,
            MaxValleyHgt = 0.1,
            MaxPercCut = 0.3,
            LowDensityRemoval = 0.1,
            RemoveMultiSD = 7,
            AlwaysClean = FALSE,
            IgnoreMonotonic = FALSE,
            MonotonicFix = NULL,
            Measures = c(1:8)
        )
    )


    ref_ff_qualityControl_flowCut <-
        readRDS(test_path("fixtures", "ff_QC_flowCut.rds"))

    # saveRDS(ff_QualityControl, test_path("fixtures", "ff_QC_flowCut.rds"))

    expect_equal(
        flowCore::exprs(ff_QualityControl),
        flowCore::exprs(ref_ff_qualityControl_flowCut)
    )
})

test_that("qualityControlFlowClean works", {
    fs_raw <- OMIP021UTSamples

    ff_QualityControl <- suppressWarnings(
        qualityControlFlowClean(fs_raw[[1]],
            binSize = 0.01, # default
            nCellCutoff = 500, # default
            cutoff = "median", # default
            fcMax = 1.3, # default
            nstable = 5
        )
    )

    ref_ff_qualityControl_flowClean <-
        readRDS(test_path("fixtures", "ff_QC_flowClean.rds"))

    # saveRDS(ff_QualityControl, test_path("fixtures", "ff_QC_flowClean.rds"))

    expect_equal(
        flowCore::exprs(ff_QualityControl),
        flowCore::exprs(ref_ff_qualityControl_flowClean)
    )
})

test_that("applyFlowJoGate works", {
  fs_raw <- OMIP021UTSamples
  
  # initial compensation first
  compMatrix <- flowCore::spillover(OMIP021UTSamples[[1]])$SPILL
  
  fs_c <- CytoPipeline::runCompensation(OMIP021UTSamples, 
                                        spillover = compMatrix)
  
  # flow jo workspace
  wspFile <- system.file("extdata",
                         "OMIP021_samples_FlowJo.wsp",
                         package = "CytoPipelineUtils")
  
  ref_ff_FJ_gated1 <-
    readRDS(test_path("fixtures", "ff_FJ_gated_cells.rds"))
  
  ff_FJ_gated1 <- applyFlowJoGate(
    fs_c[[1]],
    wspFile = wspFile,
    gateName = "Cells")
  
  #saveRDS(ff_FJ_gated1, test_path("fixtures", "ff_FJ_gated_cells.rds"))
  
  expect_equal(
    flowCore::exprs(ff_FJ_gated1),
    flowCore::exprs(ref_ff_FJ_gated1)
  )
  
  # ggplotFilterEvents(fs_c[[1]],
  #                    ff_FJ_gated1,
  #                    xChannel = "FSC-A",
  #                    yChannel = "SSC-A")
  
  ref_ff_FJ_gated2 <-
    readRDS(test_path("fixtures", "ff_FJ_gated_CD4.rds"))
  
  ff_FJ_gated2 <- applyFlowJoGate(
    ff_FJ_gated1,
    wspFile = wspFile,
    gateName = "CD4+")
  
  #saveRDS(ff_FJ_gated2, test_path("fixtures", "ff_FJ_gated_CD4.rds"))
  
  expect_equal(
    flowCore::exprs(ff_FJ_gated2),
    flowCore::exprs(ref_ff_FJ_gated2)
  )
  
  # ggplotFilterEvents(ff_FJ_gated1,
  #                    ff_FJ_gated2,
  #                    xChannel = "CD3",
  #                    yChannel = "CD4",
  #                    xScale = "logicle",
  #                    yScale = "logicle")
  
})

test_that("anonymizeMarkers works", {
    retFF <- anonymizeMarkers(OMIP021Samples[[1]],
                              oldMarkerNames = c("FSC-A","BV785 - CD3"),
                              newMarkerName = c("Fwd Scatter-A", "CD3"),
                              newExperimentName = "My experiment")
    
    checkMkName <- getChannelNamesFromMarkers(retFF, markers = "Fwd Scatter-A")
    expect_equal(checkMkName, "FSC-A")
    expect_equal(flowCore::keyword(retFF, "$P1S")[["$P1S"]], "Fwd Scatter-A")
    
    checkMkName <- getChannelNamesFromMarkers(retFF, markers = "CD3")
    expect_equal(checkMkName, "670/30Violet-A")
    expect_equal(flowCore::keyword(retFF, "$P10S")[["$P10S"]], "CD3")
    
    expect_equal(
        flowCore::keyword(
            retFF, "EXPERIMENT NAME")[["EXPERIMENT NAME"]], "My experiment")
})

