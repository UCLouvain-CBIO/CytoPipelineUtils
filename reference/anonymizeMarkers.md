# anonymizing some markers of a flowFrame

: in a flowCore::flowFrame, update some of the marker names. This
translated into :

- updating the desc field of the parameters pheno data dataframe

- updating the corresponding keyword value in the flowFrame.

## Usage

``` r
anonymizeMarkers(
  ff,
  oldMarkerNames,
  newMarkerNames,
  toUpdateKeywords = NULL,
  ...
)
```

## Arguments

- ff:

  a flowCore::flowFrame

- oldMarkerNames:

  the marker names to be amended, or alternatively, the channel names
  for which to update the marker names

- newMarkerNames:

  the new marker names to be given to the provided `oldMarkerNames`

- toUpdateKeywords:

  a list with new keyword values

- ...:

  other arguments (not used)

## Value

a new flowCore::flowFrame with the updated marker names (and possibly,
the new experiment name)

## Examples

``` r

data(OMIP021Samples)

retFF <- anonymizeMarkers(OMIP021Samples[[1]],
                          oldMarkerNames = c("FSC-A","BV785 - CD3"),
                          newMarkerName = c("Fwd Scatter-A", "CD3"),
                          toUpdateKeywords = list(
                              "EXPERIMENT NAME" = "MyExperiment"))
```
