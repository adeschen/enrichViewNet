# Remove root term if present in the list of selected terms

Remove root term if present in the list of selected terms

## Usage

``` r
removeRootTerm(gostResult)
```

## Arguments

- gostResult:

  a `data.frame` containing the terms retained for the creation of the
  network.

## Value

a `data.frame` of selected terms without the root term.

## Author

Astrid Deschênes

## Examples

``` r

## Loading dataset containing result from an enrichment analysis done with
## gprofiler2
data(demoGOST)

## Only retained the WikiPathways results
results <- demoGOST$result[demoGOST$result$source == "WP", ]

## Remove WIKIPATHWAYS root term
enrichViewNet:::removeRootTerm(gostResult=results)
#>     query significant     p_value term_size query_size intersection_size
#> 1 query_1        TRUE 0.007158877        24         15                 2
#> 2 query_1        TRUE 0.007778592        25         15                 2
#> 3 query_1        TRUE 0.029002400       239         15                 3
#>   precision     recall   term_id source
#> 1 0.1333333 0.08333333 WP:WP4925     WP
#> 2 0.1333333 0.08000000 WP:WP3613     WP
#> 3 0.2000000 0.01255230 WP:WP2882     WP
#>                                                term_name effective_domain_size
#> 1                              Unfolded protein response                 24109
#> 2 Photodynamic therapy-induced unfolded protein response                 24109
#> 3                         Nuclear Receptors Meta-Pathway                 24109
#>   source_order   parents
#> 1          570 WP:000000
#> 2          233 WP:000000
#> 3          184 WP:000000

```
