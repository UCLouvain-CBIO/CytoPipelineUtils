# CytoPipelineUtils 0.99

## CytoPipelineUtils 0.99.6
- implemented `pkgdown` customization

## CytoPipelineUtils 0.99.5
- in applyFlowJoGate(): try to avoid using 'All Samples' generic sample group
to query flowjo workspace (in order to avoid CytoML warning when using 
'All Samples')

## CytoPipelineUtils 0.99.4
- Added `getEventNbFromFJGates()`
- In `getFlowJoLabels()`, if one of the `cellTypes` is not found in any sample 
of any group, issue warning instead of error, populates corresponding output 
with NA's

## CytoPipelineUtils 0.99.3
- `applyFlowJoGate()` can now match sample file with sample names in wsp file,
by using pattern matching, instead of full file name matching

## CytoPipelineUtils 0.99.2
- added 2nd version of `removeDebrisFlowClust()`
- `applyFlowJoGate()` can now use different flowJo wsp files, 
depending on `pData`

## CytoPipelineUtils 0.99.1

- added support for FlowJo manual gating as a CytoProcessingStep
- used file.path() instead of paste0() to construct file paths 

## CytoPipelineUtils 0.99.0

- First version
