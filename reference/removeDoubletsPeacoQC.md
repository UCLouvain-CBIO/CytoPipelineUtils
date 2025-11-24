# remove doublets from a flowFrame, using PeacoQC

wrapper around PeacoQC::RemoveDoublets(). Can apply the PeacoQC function
subsequently on several channel pairs, e.g. (FSC-A, FSC-H) and (SSC-A,
SSC-H)

## Usage

``` r
removeDoubletsPeacoQC(
  ff,
  areaChannels,
  heightChannels,
  nmads = rep(4, length(areaChannels)),
  verbose = TRUE,
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

- nmads:

  a numeric vector with the bandwidth above the ratio allowed, per
  channels pair (cells are kept if the ratio between -A channel\[i\] and
  -H channel\[i\] is smaller than the median ratio + nmad\[i\] times the
  median absolute deviation of the ratios). Default is 4, for all
  channel pairs.

- verbose:

  If set to TRUE, the median ratio and width will be printed.

- ...:

  additional parameters passed to PeacoQC::RemoveDoublets()

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
    removeDoubletsPeacoQC(ff_c,
                          areaChannels = c("FSC-A", "SSC-A"),
                          heightChannels = c("FSC-H", "SSC-H"),
                          nmads = c(3, 5))                            
#> Median ratio: 1.11051333489904, width: 0.146677138600654
#> Median ratio: 1.07502856441465, width: 0.15333572575707
```
