# Validate arguments passed to createEnrichMapMultiBasic() function

Validate the arguments passed to createEnrichMapMultiBasic() function.
First, the object containing the enrichment results must correspond to a
object created by `gprofiler2` software. Second, the selected source
must at least have one enriched term in the results. Then, if the source
is 'TERM_ID', the listed terms must be present in the enrichment
results.

## Usage

``` r
validateCreateEnrichMapMultiBasicGOSTArgs(
  gostObjectList,
  queryList,
  source,
  termIDs,
  removeRoot
)
```

## Arguments

- gostObjectList:

  a `list` of `gprofiler2` objects that contain the results from an
  enrichment analysis. The list must contain at least 2 entries. The
  number of entries must correspond to the number of entries for the
  `queryList` parameter.

- queryList:

  a `list` of `character` strings representing the names of the queries
  that are going to be used to generate the graph. The query names must
  exist in the associated `gostObjectList` objects and follow the same
  order. The number of entries must correspond to the number of entries
  for the `gostObjectList` parameter.

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

## Value

`TRUE` when all arguments are valid

## Author

Astrid Deschênes

## Examples

``` r

## Load the result of an enrichment analysis done with gprofiler2
data(parentalNapaVsDMSOEnrichment)
data(rosaNapaVsDMSOEnrichment)

## Check that all arguments are valid
enrichViewNet:::validateCreateEnrichMapMultiBasicGOSTArgs(
    gostObjectList=list(parentalNapaVsDMSOEnrichment, 
                            rosaNapaVsDMSOEnrichment),
    queryList=list("parental_napa_vs_DMSO", "rosa_napa_vs_DMSO"), 
    source="GO:BP", termIDs=NULL, removeRoot=TRUE)
#> [1] TRUE
```
