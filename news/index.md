# Changelog

## CytoPipelineUtils 0.99

### CytoPipelineUtils 0.99.9

- [`applyFlowJoGate()`](https://uclouvain-cbio.github.io/CytoPipelineUtils/reference/applyFlowJoGate.md)
  now accept a Flow Jo group name to look for samples. This can help to
  discriminate files in the nasty case when some sample files have the
  same base name.

### CytoPipelineUtils 0.99.8

- migrated to GHA cache v4

### CytoPipelineUtils 0.99.7

- implemented `pkgdown` customization

### CytoPipelineUtils 0.99.6

- [`anonymizeMarkers()`](https://uclouvain-cbio.github.io/CytoPipelineUtils/reference/anonymizeMarkers.md)
  can now work with no marker anonymization (only keyword changes)
- implemented
  [`addFlowJoGatesInfo()`](https://uclouvain-cbio.github.io/CytoPipelineUtils/reference/addFlowJoGatesInfo.md)

### CytoPipelineUtils 0.99.5

- in
  [`applyFlowJoGate()`](https://uclouvain-cbio.github.io/CytoPipelineUtils/reference/applyFlowJoGate.md):
  try to avoid using ‘All Samples’ generic sample group to query flowjo
  workspace (in order to avoid CytoML warning when using ‘All Samples’)

### CytoPipelineUtils 0.99.4

- Added
  [`getEventNbFromFJGates()`](https://uclouvain-cbio.github.io/CytoPipelineUtils/reference/getEventNbFromFJGates.md)
- In
  [`getFlowJoLabels()`](https://uclouvain-cbio.github.io/CytoPipelineUtils/reference/getFlowJoLabels.md),
  if one of the `cellTypes` is not found in any sample of any group,
  issue warning instead of error, populates corresponding output with
  NA’s

### CytoPipelineUtils 0.99.3

- [`applyFlowJoGate()`](https://uclouvain-cbio.github.io/CytoPipelineUtils/reference/applyFlowJoGate.md)
  can now match sample file with sample names in wsp file, by using
  pattern matching, instead of full file name matching

### CytoPipelineUtils 0.99.2

- added 2nd version of
  [`removeDebrisFlowClust()`](https://uclouvain-cbio.github.io/CytoPipelineUtils/reference/removeDebrisFlowClust.md)
- [`applyFlowJoGate()`](https://uclouvain-cbio.github.io/CytoPipelineUtils/reference/applyFlowJoGate.md)
  can now use different flowJo wsp files, depending on `pData`

### CytoPipelineUtils 0.99.1

- added support for FlowJo manual gating as a CytoProcessingStep
- used file.path() instead of paste0() to construct file paths

### CytoPipelineUtils 0.99.0

- First version
