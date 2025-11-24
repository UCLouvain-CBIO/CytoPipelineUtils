# remove margin events, using flowAI

wrapper around flowAI::flow_auto_qc(), with inputs directed towards
specifically remove the margin events. In the current implementation,
all the signal channels, i.e. both scatter and fluo channels are
scanned.

## Usage

``` r
removeMarginsFlowAI(x, ...)
```

## Arguments

- x:

  a flowCore::flowFrame or a flowCore::flowSet

- ...:

  additional parameters passed to flowAI::flow_auto_qc(), apart from the
  following ones : remove_from, output, ChExcludeFM, html_report,
  mini_report, fcs_QC, fcs_highQ, fcs_lowQ, folder_results

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

ff_m <- removeMarginsFlowAI(x = fsRaw[[2]])
#> Quality control for the file: Donor2
#> 10.44% of anomalous cells detected in the flow rate check. 
#> 0% of anomalous cells detected in signal acquisition check. 
#> 0.1% of anomalous cells detected in the dynamic range check. 
    
```
