# Using functional enrichment results in gprofiler2 format to create an enrichment map with multiple groups from same or different enrichment analyses

User selected enrichment terms are used to create an enrichment map. The
selection of the term can by specifying by the source of the terms
(GO:MF, REAC, TF, etc.) or by listing the selected term IDs. The map is
only generated when there is at least on significant term to graph.

## Usage

``` r
createEnrichMapMultiComplex(
  gostObjectList,
  queryInfo,
  showCategory = 30L,
  categoryLabel = 1,
  categoryNode = 1,
  line = 1,
  ...
)
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

- showCategory:

  a positive `integer` or a `vector` of `characters` representing terms.
  If a `integer`, the first `n` terms will be displayed. If `vector` of
  terms, the selected terms will be displayed. Default: `30L`.

- categoryLabel:

  a positive `numeric` representing the amount by which plotting
  category nodes label size should be scaled relative to the default
  (1). Default: `1`.

- categoryNode:

  a positive `numeric` representing the amount by which plotting
  category nodes should be scaled relative to the default (1). Default:
  `1`.

- line:

  a non-negative `numeric` representing the scale of line width.
  Default: `1`.

- ...:

  additional arguments that will be pass to the
  [`emapplot`](https://rdrr.io/pkg/enrichplot/man/emapplot.html)
  function.

## Value

a `ggplot` object which is the enrichment map for enrichment results.

## Author

Astrid Deschênes

## Examples

``` r

## Loading dataset containing results from 2 enrichment analyses done with
## gprofiler2
data(parentalNapaVsDMSOEnrichment)
data(rosaNapaVsDMSOEnrichment)

## Create list of enrichment objects required to generate the enrichment map
gostObjectList=list(parentalNapaVsDMSOEnrichment, 
    parentalNapaVsDMSOEnrichment, rosaNapaVsDMSOEnrichment, 
    rosaNapaVsDMSOEnrichment)
    
## Create data frame containing required information enabling the 
## selection of the retained enriched terms for each enrichment analysis.
## One line per enrichment analyses present in the gostObjectList parameter
## With this data frame, the enrichment results will be split in 4 groups:
## 1) KEGG significant terms from parental napa vs DMSO (no root term)
## 2) REACTOME significant terms from parental napa vs DMSO (no root term)
## 3) KEGG significant terms from rosa napa vs DMSO (no root term)
## 4) REACTOME significant terms from rosa napa vs DMSO (no root term)
queryDataFrame <- data.frame(queryName=c("parental_napa_vs_DMSO", 
        "parental_napa_vs_DMSO", "rosa_napa_vs_DMSO", "rosa_napa_vs_DMSO"), 
    source=c("KEGG", "REAC", "KEGG", "REAC"), 
    removeRoot=c(TRUE, TRUE, TRUE, TRUE), termIDs=c("", "", "", ""), 
    groupName=c("parental - KEGG", "parental - Reactome", 
        "rosa - KEGG", "rosa - Reactome"), stringsAsFactors=FALSE)
    
## Create graph for KEGG and REACTOME significant results from 
## 2 enrichment analyses
createEnrichMapMultiComplex(gostObjectList=gostObjectList, 
    queryInfo=queryDataFrame, line=1.5)
#> Coordinate system already present.
#> ℹ Adding new coordinate system, which will replace the existing one.
#> Warning: ggrepel: 10 unlabeled data points (too many overlaps). Consider increasing max.overlaps

```
