# Validate arguments related to GOST passed to validateCreateEnrichMapMultiComplexArg() function

Validate the arguments passed to
validateCreateEnrichMapMultiComplexArg() function. First, the object
containing the enrichment results must correspond to a object created by
`gprofiler2` software. Second, the selected source must at least have
one enriched term in the results. Then, if the source is 'TERM_ID', the
listed terms must be present in the enrichment results.

## Usage

``` r
validateCreateEnrichMapMultiComplexGostPartTwo(gostObjectList, queryInfo)
```

## Arguments

- gostObjectList:

  a `list` of `gprofiler2` objects that contain the results from an
  enrichment analysis. The list must contain at least 2 entries. The
  number of entries must correspond to the number of entries for the
  `queryList` parameter.

- queryInfo:

  a `data.frame` contains one row per group being displayed. The number
  of rows must correspond to the number of entries for the
  `gostObjectList` parameter. The mandatory columns are:

  - `queryName`: a `character` string representing the name of the query
    retained for this group). The query names must exist in the
    associated `gostObjectList` objects and follow the same order.

  - `source`: a `character` string representing the selected source that
    will be used to generate the network. To hand-pick the terms to be
    used, "TERM_ID" should be used and the list of selected term IDs
    should be passed through the `termIDs` parameter. The possible
    sources are "GO:BP" for Gene Ontology Biological Process, "GO:CC"
    for Gene Ontology Cellular Component, "GO:MF" for Gene Ontology
    Molecular Function, "KEGG" for Kegg, "REAC" for Reactome, "TF" for
    TRANSFAC, "MIRNA" for miRTarBase, "CORUM" for CORUM database, "HP"
    for Human phenotype ontology and "WP" for WikiPathways. Default:
    "TERM_ID".

  - `removeRoot`: a `logical` that specified if the root terms of the
    selected source should be removed (when present).

  - `termIDs`: a `character` strings that contains the term IDS retained
    for the creation of the network separated by a comma ',' when the
    "TERM_ID" source is selected. Otherwise, it should be a empty string
    ("").

  - `groupName`: a `character` strings that contains the name of the
    group to be shown in the legend. Each group has to have a unique
    name.

## Value

`TRUE` when all arguments are valid

## Author

Astrid Deschênes

## Examples

``` r

## Load the result of an enrichment analysis done with gprofiler2
data(parentalNapaVsDMSOEnrichment)
data(rosaNapaVsDMSOEnrichment)

queryDataFrame <- data.frame(queryName=c("parental_napa_vs_DMSO", 
    "rosa_napa_vs_DMSO"), source=c("KEGG", "WP"), removeRoot=c(TRUE, TRUE),
    termIDs=c("", ""), groupName=c("parental - KEGG", "rosa - WP"), 
    stringsAsFactors=FALSE)

## Check that all arguments are valid
enrichViewNet:::validateCreateEnrichMapMultiComplexGostPartTwo(
    gostObjectList=list(parentalNapaVsDMSOEnrichment, 
                            rosaNapaVsDMSOEnrichment),
    queryInfo=queryDataFrame)
#> [1] TRUE
```
