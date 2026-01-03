# Extract node and edge information to be used to create network when interaction column is present

Create a list containing all node and edge information needed to create
the network

## Usage

``` r
extractInformationWhenIntersection(gostResults)
```

## Arguments

- gostResults:

  a `data.frame` containing the terms retained for the creation of the
  network. The `data.frame` must contain a column called "intersection".

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
    enrichViewNet:::extractInformationWhenIntersection(
        gostResults=results)
```
