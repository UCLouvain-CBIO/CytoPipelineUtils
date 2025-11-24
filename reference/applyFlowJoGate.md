# perform manual gating from a FlowJo gate

perform manual gating by: reading a flowjo workspace file using the
CytoML package, and apply one of the gates to the input flowFrame. Note
that not only the selected gate, but also all its hierarchy of gates,
will be applied to the input dataset.

## Usage

``` r
applyFlowJoGate(ff, gateName, wspFile, FJGroupName = NULL, ...)
```

## Arguments

- ff:

  a flowCore::flowFrame

- gateName:

  the name of the flowJo gate that will be applied

- wspFile:

  a flowjo workspace

- FJGroupName:

  optional flowjo group name in which to look for the sample file. This
  can be helpful in case of duplicate file names.

- ...:

  Extra arguments passed to
  [`getFlowJoLabels()`](https://uclouvain-cbio.github.io/CytoPipelineUtils/reference/getFlowJoLabels.md)

## Value

a flowFrame containing only the events inside the gate
