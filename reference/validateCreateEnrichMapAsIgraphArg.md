# Validate arguments passed to createEnrichMapAsIgraph() function

Validate the arguments passed to createEnrichMapAsIgraph()
function.First, the object containing the enrichment results must
correspond to a object created by `gprofiler2` software. Second, the
selected source must at least have one enriched term in the results.
Then, if the source is 'TERM_ID', the listed terms must be present in
the enrichment results.

## Usage

``` r
validateCreateEnrichMapAsIgraphArg(
  gostObject,
  query,
  source,
  termIDs,
  removeRoot,
  showCategory,
  similarityCutOff
)
```

## Arguments

- gostObject:

  a `list` created by `gprofiler2` that contains the results from an
  enrichment analysis.

- query:

  a `character` string representing the name of the query that is going
  to be used to generate the graph. The query must exist in the
  `gostObject` object.

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

  a `vector` of `character` strings that contains the term IDs retained
  for the creation of the network. This parameter is only used when
  `source` is set to "TERM_ID".

- removeRoot:

  a `logical` that specified if the root terms of the selected source
  should be removed (when present).

- showCategory:

  a positive `integer` or a `vector` of `characters` representing terms.
  If a `integer`, the first `n` terms will be displayed. If `vector` of
  terms, the selected terms will be displayed.

- similarityCutOff:

  a positive `numeric` between 0 and 1 indicating the minimum level of
  similarity between two terms to have an edge linking the terms.

## Value

`TRUE` when all arguments are valid

## Author

Astrid Deschênes

## Examples

``` r

## Load the result of an enrichment analysis done with gprofiler2
data(demoGOST)

## Check that all arguments are valid
enrichViewNet:::validateCreateEnrichMapAsIgraphArg(
    gostObject=demoGOST, query="query_1", source="GO:BP", termIDs=NULL, 
    removeRoot=FALSE, showCategory=20, similarityCutOff=0.5)
#> [1] TRUE
```
