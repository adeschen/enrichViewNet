# Change query name when more than one queries have the same name

The function adds an unique ID at the name of the query when more than
one query have the same name.

## Usage

``` r
manageQueryDuplicationInEmap(queryList)
```

## Arguments

- queryList:

  a `list` containing the query names.

## Value

a `list` with the query names modified for queries with duplicated
names.

## Author

Astrid Deschênes

## Examples

``` r

## List of query names with duplicated names
queryList <- list("parental_vs_DMSO", "rosa_vs_DMSO", "parental_vs_DMSO", 
    "rosa_vs_DMSO", "parental_vs_Control", "rosa_vs_DMSO")

## Change the query names for the duplicated names
enrichViewNet:::manageQueryDuplicationInEmap(queryList=queryList)
#> [[1]]
#> [1] "parental_vs_DMSO (1)"
#> 
#> [[2]]
#> [1] "rosa_vs_DMSO (1)"
#> 
#> [[3]]
#> [1] "parental_vs_DMSO (2)"
#> 
#> [[4]]
#> [1] "rosa_vs_DMSO (2)"
#> 
#> [[5]]
#> [1] "parental_vs_Control"
#> 
#> [[6]]
#> [1] "rosa_vs_DMSO (3)"
#> 
    
```
