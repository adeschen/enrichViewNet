# Extract node and edge information from the enrichment results

Create a list containing all node and edge information needed to create
the network

## Usage

``` r
extractNodesAndEdgesInformation(gostResults, gostObject)
```

## Arguments

- gostResults:

  a `data.frame` containing the terms retained for the creation of the
  network. The `data.frame` does not contain a column called
  "intersection".

- gostObject:

  a `list` created by gprofiler2 that contains the results from an
  enrichment analysis.

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
data(parentalNapaVsDMSOEnrichment)

## Only retained the GO Molecular Function results
results <- parentalNapaVsDMSOEnrichment$result[
        parentalNapaVsDMSOEnrichment$result$source == "GO:MF", ]

## Extract node and edge information
information <-
    enrichViewNet:::extractNodesAndEdgesInformation(
        gostResults=results, gostObject=parentalNapaVsDMSOEnrichment)
```
