# Filter the results to retain only the selected terms

Filter the enrichment results to retain only the selected terms and
remove root term if requested.

## Usage

``` r
filterResults(gostResults, source, termIDs, removeRoot)
```

## Arguments

- gostResults:

  a `data.frame` containing the enriched terms that should be filtered.

- source:

  a `character` string representing the selected source that will be
  used to generate the network. To hand-pick the terms to be used,
  "TERM_ID" should be used and the list of selected term IDs should be
  passed through the `termIDs` parameter. The possible sources are
  "GO:BP" for Gene Ontology Biological Process, "GO:CC" for Gene
  Ontology Cellular Component, "GO:MF" for Gene Ontology Molecular
  Function, "KEGG" for Kegg, "REAC" for Reactome, "TF" for TRANSFAC,
  "MIRNA" for miRTarBase, "CORUM" for CORUM database, "HP" for Human
  phenotype ontology and "WP" for WikiPathways.

- termIDs:

  a `vector` of `character` strings that contains the term IDS retained
  for the creation of the network or `NULL`.

- removeRoot:

  a `logical` that specified if the root terms of the selected source
  should be removed (when present).

## Value

a `data.frame` of filtered terms with or without the root term.

## Author

Astrid Deschênes

## Examples

``` r

## Loading dataset containing result from an enrichment analysis done with
## gprofiler2
data(demoGOST)

## Only retained the GO - Molecular Function results
results <- demoGOST$result

## Remove WIKIPATHWAYS root term
selected <- enrichViewNet:::filterResults(gostResults=results, source="WP", 
    termIDs=NULL, removeRoot=TRUE)

```
