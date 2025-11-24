# perform manual gating from a FlowJo gate

perform manual gating by: reading a flowjo workspace file using the
flowWorkspace library, and apply the selected gates to the input
flowFrame. Then provides nb of remaining events for the selected gates.

## Usage

``` r
getEventNbFromFJGates(
  ff,
  wspFile,
  gates = NULL,
  defaultGates = c("leaves", "all"),
  ...
)
```

## Arguments

- ff:

  a flowCore::flowFrame

- wspFile:

  a flowjo workspace

- gates:

  vector of flowJo gate nodes to parse (should be common to all
  groups!). Default NULL then uses either all leaf node from first
  group, or all gates from first group, depending on 'defaultGates'.

- defaultGates:

  if gates are not specified, than triggers the choice of either :

  - all leaf nodes ('leaves') (the default)

  - all nodes ('all')

- ...:

  Extra arguments passed to
  [`getFlowJoLabels()`](https://uclouvain-cbio.github.io/CytoPipelineUtils/reference/getFlowJoLabels.md)

## Value

a dataframe with the nb of events per gate, including one additional
'unlabelled' gate for un-labelled events.
