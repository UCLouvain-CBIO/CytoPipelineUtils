test_that("getFlowJoLabels works", {
  data(OMIP021Samples)
  wspFile <- system.file("extdata",
                         "OMIP021_samples_FlowJo.wsp",
                         package = "CytoPipelineUtils")
  compMatrix <- flowCore::spillover(OMIP021Samples[[1]])$SPILL
  
  group <- "Donors"
  cellTypes <- c("Cells","Debris")
  
  fs_c <- CytoPipeline::runCompensation(OMIP021Samples, 
                                        spillover = compMatrix)
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
  
  expect_error(getFlowJoLabels(agg,
                               wspFile = wspFile,
                               groups = groups,
                               cellTypes = cellTypes),
               regexp = "not found in any sample of any group")
  
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