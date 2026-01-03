# Extract information about nodes and edges to be used to create network when interaction column is missing

Create a list containing all node and edge information needed to create
the network

## Usage

``` r
extractInformationWhenNoIntersection(gostResults, gostObject)
```

## Arguments

- gostResults:

  a `data.frame` containing the terms retained for the creation of the
  network. The `data.frame` does not contain a column called
  "intersection".

- gostObject:

  a `list` created by gprofiler2 that contains the results and meta-data
  from an enrichment analysis.

## Value

`list` containing 2 entries:

- `"geneNodes"`: a `data.frame` containing the information about the
  nodes present in the network. The nodes are genes.

- `"termNodes"`: a `data.frame` containing the information about the
  nodes present in the network. The nodes are terms.

- `"edges"`: a `data.frame` containing the information about the edges
  present in the network. The edges connect one gene to one term.

## Author

Astrid Deschênes

## Examples

``` r

## Loading dataset containing result from an enrichment analysis done with
## gprofiler2
data(demoGOST)

## Only retained the WikiPathways results
results <- demoGOST$result[demoGOST$result$source == "WP", ]

information <- enrichViewNet:::extractInformationWhenNoIntersection(
                gostResults=results, gostObject=demoGOST)
```
