# Small soil phylum data for examples, sample data as data frame combined with counts for OTU 1 and sequencing depth.

A small subset of
[`soil_phylo_sample`](https://statdivlab.github.io/corncob/reference/soil_phylo_sample.md)
used for examples. A data frame made from the \`phyloseq\` object with
only sample data and counts for OTU 1.

## Usage

``` r
soil_phylum_small_otu1
```

## Format

A phyloseq-class experiment-level object with sample data and OTU 1
counts.

- sam_data:

  sample data with the following covariates:

  - `Plants`, values `0` and `1`. Index for different plants

  - `Day`, values `0` (initial sampling point), `1` (12 days after
    treatment additions), and `2` (82 days after treatment additions).
    Index for different days of measurement

  - `Amdmt`, values `0` (no additions), `1` (biochar additions), and `2`
    (fresh biomass additions). Index for different soil additives.

  - `DayAmdmt`, values `00`, `01`, `02`, `10`, `11`, `12`, `20`, `21`,
    and `22`. A single index for the combination of `Day` and `Amdmt`
    with `Day` as the first digit and `Amdmt` as the second digit.

  - `ID`, values `A`, `B`, `C`, `D`, and `F`. Index for different soil
    plots.

  - `W`, counts for OTU1 in each sample. This OTU corresponds with the
    phylum *Proteobacteria*.

  - `M`, the sequencing depth for each sample.

## References

Whitman, T., Pepe-Ranney, C., Enders, A., Koechli, C., Campbell, A.,
Buckley, D. H., Lehmann, J. (2016). *Dynamics of microbial community
composi-tion and soil organic carbon mineralization in soil following
addition of pyrogenic andfresh organic matter*. The ISME journal,
10(12):2918. \<doi: 10.1038/ismej.2016.68\>.
