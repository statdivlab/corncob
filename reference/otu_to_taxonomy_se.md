# Transform OTUs to their taxonomic label

Transform OTUs to their taxonomic label

## Usage

``` r
otu_to_taxonomy_se(OTU, data, level = NULL)
```

## Arguments

- OTU:

  String vector. Names of OTU labels in `data`

- data:

  `phyloseq` object with a taxonomy table

- level:

  (Optional). Character vector. Desired taxonomic levels for output.

## Value

String vector. Names of taxonomic labels matching labels of `OTU`.
