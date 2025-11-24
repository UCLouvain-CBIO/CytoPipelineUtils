# remove doublets from a flowFrame, using flowStats

Wrapper around flowStats::singletGate(). Can apply the flowStats
function subsequently on several channel pairs, e.g. (FSC-A, FSC-H) and
(SSC-A, SSC-H)

## Usage

``` r
removeDoubletsFlowStats(
  ff,
  areaChannels,
  heightChannels,
  widerGate = FALSE,
  ...
)
```

## Arguments

- ff:

  a flowCore::flowFrame

- areaChannels:

  a character vector containing the name of the 'area type' channels one
  wants to use

- heightChannels:

  a character vector containing the name of the 'height type' channels
  one wants to use

- widerGate:

  a boolean as wider_gate parameter to flowStats::singletGate()

- ...:

  additional parameters passed to flowStats::singletGate()

## Value

a flowCore::flowFrame with removed doublets events from the input

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

ff_s <-
    removeDoubletsFlowStats(ff_c,
                            areaChannels = c("FSC-A", "SSC-A"),
                            heightChannels = c("FSC-H", "SSC-H"),
                            widerGate = TRUE)
#> Warning: replacing previous import ‘flowViz::contour’ by ‘graphics::contour’ when loading ‘flowStats’
                            
```
