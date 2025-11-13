# Function to subset and convert SummarizedExperiment data

Function to subset and convert SummarizedExperiment data

## Usage

``` r
convert_sumexp(data, select)
```

## Arguments

- data:

  a `SummarizedExperiment` object

- select:

  Name of OTU or taxa to select, must match taxa name in `data`

## Value

A `data.frame` object, with elements `W` as the observed counts, `M` as
the sequencing depth, and the sample data with their original names.
