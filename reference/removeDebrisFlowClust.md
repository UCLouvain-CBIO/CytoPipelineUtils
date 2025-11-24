# remove debris from a flowFrame, using flowClust

this function removes debris from a flowFrame, using clustering
capabilities of flowClust::tmixFilter(). The idea is to pre-select a
number of clusters to be found in the (FSC,SSC) 2D view, and select all
clusters behalve the one closest nearest the origin. Then we take all
events that are inside the `probaLevel` quantile curve for at least one
of the remaining clusters.

## Usage

``` r
removeDebrisFlowClust(
  ff,
  FSCChannel,
  SSCChannel,
  nClust,
  probaLevel = 0.9,
  ...
)
```

## Arguments

- ff:

  a flowCore::flowFrame

- FSCChannel:

  the name of the FSC channel

- SSCChannel:

  the name of the SSC channel

- nClust:

  number of clusters to identify

- probaLevel:

  the probability level

- ...:

  additional parameters passed to flowClust::flowClust()

## Value

a flowCore::flowFrame with removed debris events from the input

## Examples

``` r

rawDataDir <- system.file("extdata", package = "CytoPipeline")
sampleFiles <-
    file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))

truncateMaxRange <- FALSE
minLimit <- NULL

# create flowCore::flowSet with all samples of a dataset
fsRaw <- readSampleFiles(
    sampleFiles = sampleFiles,
    whichSamples = "all",
    truncate_max_range = truncateMaxRange,
    min.limit = minLimit)

suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#> Removing margins from file : Donor2.fcs
    
ff_c <-
    compensateFromMatrix(ff_m,
                         matrixSource = "fcs")        


ff_cells <-
    removeDebrisFlowClust(ff_c,
                          FSCChannel = "FSC-A",
                          SSCChannel = "SSC-A",
                          nClust = 3,
                          probaLevel = 0.9,
                          B = 100)
#> Using the serial version of flowClust
```
