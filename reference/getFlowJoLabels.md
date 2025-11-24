# Process a flowjo workspace file

Reads a flowjo workspace file using the flowWorkspace library and
returns a vector with a label for each cell of a flowFrame from a set of
specified gates. If the flowFrame was generated from several files
(contains a 'File' column containing a file index), then the groups
should contain the group name to apply for each of these files/samples

## Usage

``` r
getFlowJoLabels(
  ff,
  wspFile,
  groups = "All Samples",
  sampleInGroups = NULL,
  cellTypes = NULL,
  defaultCellTypes = c("leaves", "all"),
  specialCharsInChannels = c("/"),
  FlowJoCharsinChannels = c("_"),
  withFJv10TimeCorrection = TRUE,
  ...
)
```

## Arguments

- ff:

  a flowCore::flowFrame

- wspFile:

  a flowjo workspace

- groups:

  vector of flowjo groups to parse from the flowjo workspace (i.e. one
  for each sample/file used to generate the flowFrame). Default c("All
  Samples")

- sampleInGroups:

  vector of flowjo sample index in group to use from the flowjo
  workspace (i.e. one for each sample/file used to generate the
  flowFrame). If not provided, will start from 1 in each group and
  increment. Default NULL.

- cellTypes:

  vector of flowJo gate nodes to parse (should be common to all
  groups!). Default NULL then uses either all leaf node from first
  group, or all gates from first group, depending on 'defaultCellTypes'.

- defaultCellTypes:

  if cellTypes are not specified, than triggers the choice of either :

  - all leaf nodes ('leaves') (the default)

  - all nodes ('all')

- specialCharsInChannels:

  vector of special char litterals that are replaced in FlowJo workspace

- FlowJoCharsinChannels:

  vector of new char litterals to be applied in FlowJo workspace

- withFJv10TimeCorrection:

  if TRUE (default), applies a time correction to cope with a bug in
  flowJo 10.0. Should be kept TRUE with FloJo v\>10.0 if time channel is
  used in gating.

  Since FlowJo v10, time channel is read differently from fcs files (see
  https://docs.flowjo.com/flowjo/experiment-based-platforms/kinetics/
  how-flowjo-handles-time/). The resulting scale transformation on the
  time channel is normally stored in the corresponding tranformation
  gain parameter in the FlowJo workspace. However, somehow this gain
  parameter can be sometimes wrong in the .wsp file, and so the
  workaround consists in correcting the time channel values of the
  cytoset prior to gating, in order to have the gate labels calculated
  correctly.

- ...:

  Extra arguments passed to
  [`CytoML::flowjo_to_gatingset()`](https://rdrr.io/pkg/CytoML/man/flowjo_to_gatingset.html)

## Value

a list with :

- the first element ("matrix") is a matrix containing filtering results
  for each specified gate

- the second element ("labels") is a vector which assigns one label to
  each cell. If no cell type correspond to a specific cell, than the
  keyword 'unlabeled' is assigned to this cell. If the cell belongs to
  several gates (meaning that the gates ar not disjoints), than this
  cell is assigned to the gate with the less matching cells.
