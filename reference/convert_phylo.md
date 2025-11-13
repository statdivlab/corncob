# Function to subset and convert phyloseq data

Function to subset and convert phyloseq data

## Usage

``` r
convert_phylo(data, select)
```

## Arguments

- data:

  a `phyloseq` object

- select:

  Name of OTU or taxa to select, must match taxa name in `data`

## Value

A `data.frame` object, with elements `W` as the observed counts, `M` as
the sequencing depth, and the sample data with their original names.
