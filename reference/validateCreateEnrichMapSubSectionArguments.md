# Validate common arguments passed to createEnrichMap() and createEnrichMapMultiBasic() functions

Validate the common arguments passed to createEnrichMap() and
createEnrichMapMultiBasic() functions.

## Usage

``` r
validateCreateEnrichMapSubSectionArguments(
  showCategory,
  categoryLabel,
  categoryNode,
  line
)
```

## Arguments

- showCategory:

  a positive `integer` or a `vector` of `characters` representing terms.
  If a `integer`, the first `n` terms will be displayed. If `vector` of
  terms, the selected terms will be displayed.

- categoryLabel:

  a positive `numeric` representing the amount by which plotting
  category nodes label size should be scaled relative to the default
  (1).

- categoryNode:

  a positive `numeric` representing he amount by which plotting category
  nodes should be scaled relative to the default (1).

- line:

  a non-negative `numeric` representing the scale of line width.

## Value

`TRUE` when all arguments are valid

## Author

Astrid Deschênes

## Examples

``` r

## Check that all arguments are valid
enrichViewNet:::validateCreateEnrichMapSubSectionArguments(
    showCategory=20, categoryLabel=1.1, categoryNode=1, 
    line=0.5)
#> [1] TRUE
```
