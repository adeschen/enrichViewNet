# Change name description in data frame when more than one term as the same name

The function adds the ID of the term at the name of the name description
when terms with different ID share the same name. For example, "MAPK
pathway" would become "MAPK pathway (KEGG:04010)".

## Usage

``` r
manageNameDuplicationInEmap(clProfDF)
```

## Arguments

- clProfDF:

  a `data.frame` containing the enrichment terms that are going to be
  graphed as an enrichment map. The `data.frame` must contain at least
  one case of duplicated name descriptions for different term IDs.

## Value

a `data.frame` with the name descriptions modified for terms with
duplicated name descriptions but different term IDs.

## Author

Astrid Deschênes

## Examples

``` r

## Data frame with duplicated name descriptions for different IDs
clustData <- data.frame(Cluster=c("group 1" , "group 1", "group 2"), 
    ID=c("WP:WP4925", "WP:WP382", "KEGG:04010"),
    Description=c("Unfolded protein response", 
                    rep("MAPK signaling pathway", 2)),
    GeneRatio=c("4/157", "3/157", "3/157"),
    BgRatio=c("4/24022", "3/24022", "3/24022"),
    pvalues=c(1.55e-4, 8.13e-8, 4.33e-5),
    p.adjust=c(1e-3, 1e-3, 1.4e-3), qvalue=c(1e-3, 1e-3, 1.4e-3), 
    geneID=c("ENSG000107968/ENSG000120129/ENSG000123358/ENSG000158050",
        "ENSG000107968/ENSG000120129/ENSG000158050",
        "ENSG000107968/ENSG000120129/ENSG000158050"),
    Count=c(4, 3, 3))

## Change the name descriptions for the duplicated terms
enrichViewNet:::manageNameDuplicationInEmap(clProfDF=clustData)
#>   Cluster         ID                         Description GeneRatio BgRatio
#> 1 group 1  WP:WP4925           Unfolded protein response     4/157 4/24022
#> 2 group 1   WP:WP382   MAPK signaling pathway (WP:WP382)     3/157 3/24022
#> 3 group 2 KEGG:04010 MAPK signaling pathway (KEGG:04010)     3/157 3/24022
#>    pvalues p.adjust qvalue
#> 1 1.55e-04   0.0010 0.0010
#> 2 8.13e-08   0.0010 0.0010
#> 3 4.33e-05   0.0014 0.0014
#>                                                    geneID Count
#> 1 ENSG000107968/ENSG000120129/ENSG000123358/ENSG000158050     4
#> 2               ENSG000107968/ENSG000120129/ENSG000158050     3
#> 3               ENSG000107968/ENSG000120129/ENSG000158050     3
    
```
