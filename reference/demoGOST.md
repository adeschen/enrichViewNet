# The result of a functional enrichment analysis done with `gprofiler2` (<https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html>).

The object is a `list` with 2 entries. It contains the results of the
enrichment analysis as well as the metadata related to the analysis.

## Usage

``` r
data(demoGOST)
```

## Format

The `list` contains two entries. The `result` entry contains a
`data.frame` with the significant results obtained by an enrichment
analysis done with `gprofiler2`. The `meta` entry contains a named list
with all the metadata for the query.

## Value

A `list` containing two entries. The `result` entry contains a
`data.frame` with the significant results obtained by an enrichment
analysis done with `gprofiler2`. The `meta` entry contains a named list
with all the metadata for the query.

## Details

This dataset can be used to test the [`createNetwork`](createNetwork.md)
function.

## See also

- [createNetwork](createNetwork.md) for transforming functional
  enrichment results from gprofiler2 into a Cytoscape network

- [createEnrichMap](createEnrichMap.md) for transforming functional
  enrichment results from gprofiler2 into an enrichment map

## Examples

``` r

## Loading dataset containing result from an enrichment analysis done with
## gprofiler2
data(demoGOST)

## Create network for WikiPathways results 
## in Cytoscape (if the application is open)
## Otherwise, create a CX file in the temporary directory
## The file can be opened in Cytoscape
createNetwork(gostObject=demoGOST, source="WP", title="Wikipathways",
    fileName=file.path(tempdir(), "Wikipathways_Demo.cx"))
#> Unable to connect to Cytoscape. 
#> CX JSON file will be created.
#> Preparing information for generating CX JSON file.
#> CX JSON file "/tmp/RtmpKsInZl/Wikipathways_Demo.cx" has been created.
#> [1] TRUE

```
