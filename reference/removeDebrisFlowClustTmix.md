# remove debris from a flowFrame, using flowClust

this function removes debris from a flowFrame, using clustering
capabilities of flowClust::tmixFilter(). The idea is to pre-select a
number of clusters to be found in the (FSC,SSC) 2D view, and eliminate
the cluster that is the closest to the origin.

## Usage

``` r
removeDebrisFlowClustTmix(ff, FSCChannel, SSCChannel, nClust, ...)
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

- ...:

  additional parameters passed to flowClust::tmixFilter()

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
    removeDebrisFlowClustTmix(ff_c,
                              FSCChannel = "FSC-A",
                              SSCChannel = "SSC-A",
                              nClust = 3,
                              level = 0.97,
                              B = 100)
#> The prior specification has no effect when usePrior=no
#> Using the serial version of flowClust
```
