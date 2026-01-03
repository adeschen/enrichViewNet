# Validate arguments passed to creatNetwork() function

Validate the arguments passed to creatNetwork() function. First, the
object containing the enrichment results must correspond to a object
created by `gprofiler2` software. Second, the selected source must at
least have one enriched term in the results. Then, if the source is
'TERM_ID', the listed terms must be present in the enrichment results.

## Usage

``` r
validateCreateNetworkArguments(
  gostObject,
  source,
  termIDs,
  removeRoot,
  query,
  title,
  collection,
  fileName
)
```

## Arguments

- gostObject:

  a `list` created by `gprofiler2` that contains the results from an
  enrichment analysis.

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

- query:

  a `character` string that specified the retained query to generate the
  network or `NULL`.

- title:

  a `character` string representing the name assigned to the network.

- collection:

  a `character` string representing the collection name assigned to the
  network.

- fileName:

  a `character` string representing the name of the CX JSON file that is
  created when Cytoscape is not running. The name must have a '.cx'
  extension.

## Value

`TRUE` when all arguments are valid

## Author

Astrid Deschênes

## Examples

``` r

## Load the result of an enrichment analysis done with gprofiler2
data(demoGOST)

## Check that all arguments are valid
enrichViewNet:::validateCreateNetworkArguments(gostObject=demoGOST,
    source="GO:BP", termIDs=NULL, removeRoot=FALSE, query=NULL, 
    title="Network graph Test",
    collection="test collection", fileName="test.cx")
#> [1] TRUE
```
