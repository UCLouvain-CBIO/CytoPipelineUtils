# remove dead cells from a flowFrame

this function removes dead cells from a flowFrame, using a specific
'(a)live/dead' channel, and the flowDensity::deGate() gating function
(see doc of the flowDensity package)

## Usage

``` r
removeDeadCellsDeGate(
  ff,
  preTransform = FALSE,
  transList = NULL,
  LDMarker,
  keepPositivePop = FALSE,
  ...
)
```

## Arguments

- ff:

  a flowCore::flowFrame

- preTransform:

  if TRUE, apply the transList scale transform prior to running the
  gating algorithm

- transList:

  applied in conjunction with preTransform == TRUE

- LDMarker:

  a character containing the exact name of the marker corresponding to
  Live/Dead channel, or the Live/Dead channel name itself

- keepPositivePop:

  logical flag stating whether we want to keep, after gating, the
  population that is positive for `LDMarker` (TRUE) or negative for
  `LDmarker` (FALSE)

- ...:

  additional parameters passed to flowDensity::deGate()

## Value

a flowCore::flowFrame with removed dead cells from the input

## Details

rawDataDir \<- system.file("extdata", package = "CytoPipeline")
sampleFiles \<- file.path(rawDataDir, list.files(rawDataDir, pattern =
"Donor"))

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

transList <- 
    estimateScaleTransforms(        
        ff = ff_c,
        fluoMethod = "estimateLogicle",
        scatterMethod = "linear",
        scatterRefMarker = "BV785 - CD3")

ff_lcells <-
    removeDeadCellsDeGate(ff_c,
                          preTransform = TRUE,
                          transList = transList,
                          LDMarker = "L/D Aqua - Viability",
                          keepPositivePop = FALSE)
#> Removing Dead Cells events (using flowDensity::deGate()) from file : Donor2.fcs
                            
```
