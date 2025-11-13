# Soil data, sample data.

A data frame made from a soil \`phyloseq\` object with only sample data.

## Usage

``` r
soil_phylo_sample
```

## Format

A phyloseq-class experiment-level object with sample data.

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

## References

Whitman, T., Pepe-Ranney, C., Enders, A., Koechli, C., Campbell, A.,
Buckley, D. H., Lehmann, J. (2016). *Dynamics of microbial community
composi-tion and soil organic carbon mineralization in soil following
addition of pyrogenic andfresh organic matter*. The ISME journal,
10(12):2918. \<doi: 10.1038/ismej.2016.68\>.
