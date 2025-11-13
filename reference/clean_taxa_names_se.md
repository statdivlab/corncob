# Rename taxa

Renames taxa to have short human-readable names

## Usage

``` r
clean_taxa_names_se(x, name = "OTU")
```

## Arguments

- x:

  Object of class `SummarizedExperiment`

- name:

  Character, defaults to `"OTU"`. Optional. String to use in every taxa
  name.

## Value

Object of class `SummarizedExperiment`, with taxa renamed (defaults to
OTU1, OTU2, ...), with the original taxa names saved as an attribute.

## Details

The original taxa names are saved as the `original_names` attribute. See
the example for an example of how to access the original names.
