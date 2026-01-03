# Transform the node and edge information into a format easy to process to create a CX JSON text file

Extract information about nodes and edges that is necessary to create
the CX JSON text representing the network

## Usage

``` r
formatInformationForCXJSON(results)
```

## Arguments

- results:

  a `list` containing the information about the nodes and edges present
  in the network. TODOS

## Value

a `list` containing 4 entries:

- `"nodes"`: a `data.frame` containing the information about the nodes
  present in the network.

- `"edges"`: a `data.frame` containing the information about the edges
  present in the network.

- `"nodeAttributes"`: a `data.frame` containing the attributes
  associated to the nodes present in the network.

- `"edgesAttributes"`: a `data.frame` containing the attributes
  associated to the edges present in the network

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
info <- enrichViewNet:::extractNodesAndEdgesInformation(gostResults=results, 
            gostObject=parentalNapaVsDMSOEnrichment)
            
## Format node and edge information
information <- enrichViewNet:::formatInformationForCXJSON(result=info)
```
