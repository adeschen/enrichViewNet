# Using functional enrichment results from gprofiler2 to create a Cytoscape network

User selected enrichment terms are used to create a Cytoscape network
where the selected terms and the genes that where part of the enrichment
analysis are all represented as nodes. Edges are linking the genes to
their terms. The selection of the term can by specifying the source of
the terms (GO:MF, REAC, TF, etc.) or by listing the selected term IDs.
The network is only generated when there is at least on significant term
to graph. When the enrichment analysis contains more than one query,
only one query can be selected to generate the network.

## Usage

``` r
createNetwork(
  gostObject,
  source = c("TERM_ID", "GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "TF", "MIRNA", "HPA",
    "CORUM", "HP", "WP"),
  termIDs = NULL,
  removeRoot = TRUE,
  query = NULL,
  title = "gprofiler network",
  collection = "enrichment results",
  fileName = "gprofilerNetwork.cx"
)
```

## Arguments

- gostObject:

  a `list` created by gprofiler2 that contains the results from an
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
  phenotype ontology and "WP" for WikiPathways. Default: "TERM_ID".

- termIDs:

  a `vector` of `character` strings that contains the term IDS retained
  for the creation of the network. Default: `NULL`.

- removeRoot:

  a `logical` that specified if the root terms of the selected source
  should be removed (when present). Default: `TRUE`.

- query:

  a `character` string that specified the retained query to generate the
  network. When `NULL`, the query present in the result is retained;
  `NULL` cannot be used when more than one query is present. Default:
  `NULL`.

- title:

  a `character` string representing the name assigned to the network.
  Default: "gprofiler network".

- collection:

  a `character` string representing the collection name assigned to the
  network. Default: "enrichment results".

- fileName:

  a `character` string representing the name of the CX JSON file that is
  created when Cytoscape is not running. The name must have a '.cx'
  extension. Default: "gprofilerNetwork_01.cx".

## Value

`TRUE`

## Author

Astrid Deschênes

## Examples

``` r

## Loading dataset containing result from an enrichment analysis done with
## gprofiler2
data(parentalNapaVsDMSOEnrichment)

## Some of the enrichment results present in the dataset
head(parentalNapaVsDMSOEnrichment$result)
#>                   query significant      p_value term_size query_size
#> 1 parental_napa_vs_DMSO        TRUE 2.198177e-11      3187        157
#> 2 parental_napa_vs_DMSO        TRUE 2.139995e-10      1238        157
#> 3 parental_napa_vs_DMSO        TRUE 2.848209e-10      3583        157
#> 4 parental_napa_vs_DMSO        TRUE 8.790372e-10        48        157
#> 5 parental_napa_vs_DMSO        TRUE 3.046614e-09      4216        157
#> 6 parental_napa_vs_DMSO        TRUE 3.550218e-09      4119        157
#>   intersection_size  precision     recall    term_id source
#> 1                59 0.37579618 0.01851271 GO:0048523  GO:BP
#> 2                35 0.22292994 0.02827141 GO:0008219  GO:BP
#> 3                61 0.38853503 0.01702484 GO:0048519  GO:BP
#> 4                10 0.06369427 0.20833333 GO:0070059  GO:BP
#> 5                65 0.41401274 0.01541746 GO:0080090  GO:BP
#> 6                64 0.40764331 0.01553775 GO:0051171  GO:BP
#>                                                                           term_name
#> 1                                           negative regulation of cellular process
#> 2                                                                        cell death
#> 3                                         negative regulation of biological process
#> 4 intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress
#> 5                                           regulation of primary metabolic process
#> 6                                 regulation of nitrogen compound metabolic process
#>   effective_domain_size source_order      parents
#> 1                 24022        13544 GO:00099....
#> 2                 24022         3218   GO:0009987
#> 3                 24022        13540 GO:00081....
#> 4                 24022        16817 GO:00349....
#> 5                 24022        18790 GO:00192....
#> 6                 24022        14316 GO:00068....
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          evidence_codes
#> 1                                                                                                                                                                            IMP,IDA IGI IBA,IDA,IDA TAS,IDA IMP TAS,IDA,IDA IMP IBA NAS,ISS,ISS,IDA IMP ISS IBA,ISS,IDA IMP ISS,IDA IMP IGI IEP IBA TAS,IDA,TAS,IDA IMP,IDA IMP ISS,IGI ISS,ISS,IMP,IMP ISS IBA,IDA IBA NAS,IDA IGI TAS,IBA,IDA IBA,IDA TAS,IDA IMP IBA TAS,ISS,IDA ISS,IDA,IDA,IBA,IDA IMP IBA TAS,IDA,IDA,IDA IMP ISS IBA IC,IMP ISS IBA,ISS,IDA IMP IBA,ISS,ISS,IDA,IMP,TAS,IDA IMP IGI ISS ISO IBA,IBA,IMP,IDA IMP IBA,IDA ISS IBA,IBA,IMP TAS,NAS,IDA,ISS,IDA ISS NAS,IDA IMP TAS,IDA IMP TAS,IDA ISS ISO TAS,IBA
#> 2                                                                                                                                                                                                                                                                                                                                                                                                                            IDA IGI IBA,TAS,IDA,IMP,IDA IMP IGI ISS IBA TAS,IDA NAS,ISS,ISS,IMP,TAS,IDA,TAS,IMP,IGI ISS TAS,IMP,IMP,IDA,IMP,TAS,IDA IMP TAS,IMP,IMP ISS TAS,NAS,IDA IMP IGI IBA TAS,IMP,IDA ISS,ISS,TAS,TAS,IDA IMP ISS IBA TAS IC,IMP,TAS,IDA IMP TAS,IDA IMP TAS,IBA
#> 3                                                                                                                                                    IMP,IDA IGI IBA,IDA,IDA TAS,IDA IMP TAS,IDA IMP,ISS,IDA IMP IBA NAS,ISS,ISS,IDA IMP ISS IBA,ISS,IDA IMP ISS,IDA IMP IGI IEP ISS IBA TAS,IDA,TAS,IDA IMP,IDA IMP ISS TAS,IGI ISS,ISS,IMP,IMP ISS IBA,IDA IBA NAS,IDA IGI TAS,IBA,IDA IBA,IDA TAS,IDA IMP IBA TAS,ISS,IDA ISS,IDA,IDA,IBA,IDA IMP IBA TAS,IDA,IDA,IMP,IDA IMP ISS IBA IC,IMP ISS IBA,ISS,IDA IMP IBA,IGI ISS,ISS,IDA,IMP,TAS,IDA IMP IGI ISS ISO IBA,IBA,IMP,IDA IMP IBA,IDA ISS IBA,IBA,IMP TAS,NAS,IDA,ISS,IDA ISS NAS,IDA IMP TAS,IDA IMP TAS,IDA ISS ISO TAS,IBA
#> 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                IDA IBA,TAS,IDA ISS IBA TAS,ISS,IMP,TAS,TAS,TAS,IDA IMP IBA TAS IC,IDA
#> 5 IDA IMP IBA,IDA IBA,IDA IBA TAS IC,ISS,IDA IBA,IDA,IDA ISS IBA,IBA,IPI,IDA IMP IBA NAS,ISS,IDA ISS IBA,IBA,IDA,IDA IMP IGI ISS IBA,IDA ISS IBA,IDA IBA,IDA IMP ISS IBA,IDA IMP IEP IBA TAS,IBA TAS,IDA IMP IBA,IDA IMP ISS,IDA IMP ISS ISO IBA NAS IC,IMP,IDA IMP,NAS,IDA IBA,IDA IMP IGI ISS IBA NAS,IDA,IDA IMP IBA,IDA,ISS,ISS IBA,TAS,IBA,IDA,IDA IBA,IDA,IDA ISS IBA NAS,IDA,IMP,IDA IMP ISS IBA,IBA,ISS IBA,IDA IGI ISS,IDA IMP IGI ISS IBA TAS NAS,IBA,IDA IBA TAS NAS,EXP IDA IMP ISS IBA NAS,IDA IMP ISS IBA,IDA IMP IGI ISS ISO IBA TAS NAS,IBA,IBA,IDA IBA,IDA ISS,IMP,IDA IBA NAS,TAS,IBA,IMP ISS IBA NAS,IDA IMP IBA,IDA IMP IBA,IDA IMP IGI ISS ISO IBA TAS NAS,IBA,IBA
#> 6     IDA IMP IBA,IDA IBA,IDA IBA TAS IC,ISS,IDA IBA,IDA,IDA ISS IBA,IBA,IPI,IDA IMP IBA NAS,ISS,IDA ISS IBA,IBA,IDA,IDA IMP IGI ISS IBA,IDA ISS IBA,IDA IBA,IDA IMP ISS IBA,IDA IMP IEP IBA TAS,IBA TAS,IDA IMP IBA,IDA IMP ISS,IDA IMP ISS ISO IBA NAS IC,IMP,IDA IMP,NAS,IDA IBA,IDA IMP IGI ISS IBA NAS,IDA,IDA IMP IBA,IDA,ISS,ISS IBA,IBA,IDA,IDA IBA,IDA,IDA ISS IBA NAS,IDA,IMP,IDA IMP ISS IBA,IBA,ISS IBA,IDA IGI ISS,IDA IMP IGI ISS IBA TAS NAS,IBA,IDA IBA TAS NAS,EXP IDA IMP ISS IBA NAS,IDA IMP ISS IBA,IDA IMP IGI ISS ISO IBA TAS NAS,IBA,IBA,IDA IBA,IDA ISS,IMP,IDA IBA NAS,TAS,IBA,IMP ISS IBA NAS,IDA IMP IBA,IDA IMP IBA,IDA IMP IGI ISS ISO IBA TAS NAS,IBA,IBA
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      intersection
#> 1                                                                                                 ENSG00000007944,ENSG00000051108,ENSG00000059728,ENSG00000087074,ENSG00000100292,ENSG00000105327,ENSG00000113916,ENSG00000116741,ENSG00000119508,ENSG00000120129,ENSG00000123358,ENSG00000124216,ENSG00000124762,ENSG00000125266,ENSG00000125740,ENSG00000125812,ENSG00000128016,ENSG00000128272,ENSG00000128590,ENSG00000128965,ENSG00000130766,ENSG00000133639,ENSG00000136826,ENSG00000138166,ENSG00000140961,ENSG00000141232,ENSG00000141582,ENSG00000142178,ENSG00000143878,ENSG00000148926,ENSG00000157557,ENSG00000158050,ENSG00000159388,ENSG00000162772,ENSG00000163273,ENSG00000163874,ENSG00000164056,ENSG00000164920,ENSG00000168140,ENSG00000169242,ENSG00000172216,ENSG00000172602,ENSG00000173334,ENSG00000173530,ENSG00000175197,ENSG00000177873,ENSG00000179388,ENSG00000183508,ENSG00000183691,ENSG00000184545,ENSG00000184557,ENSG00000185022,ENSG00000197019,ENSG00000198576,ENSG00000204103,ENSG00000204388,ENSG00000204389,ENSG00000245848,ENSG00000256223
#> 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000051108,ENSG00000087074,ENSG00000099860,ENSG00000100292,ENSG00000105327,ENSG00000113916,ENSG00000119508,ENSG00000120738,ENSG00000124216,ENSG00000124762,ENSG00000125266,ENSG00000125657,ENSG00000128016,ENSG00000128272,ENSG00000128965,ENSG00000133639,ENSG00000136826,ENSG00000140961,ENSG00000141582,ENSG00000141682,ENSG00000142178,ENSG00000143878,ENSG00000144655,ENSG00000153094,ENSG00000162772,ENSG00000163874,ENSG00000170345,ENSG00000172216,ENSG00000173530,ENSG00000175197,ENSG00000179388,ENSG00000184557,ENSG00000204388,ENSG00000204389,ENSG00000241839
#> 3                                                                 ENSG00000007944,ENSG00000051108,ENSG00000059728,ENSG00000087074,ENSG00000100292,ENSG00000105327,ENSG00000113369,ENSG00000113916,ENSG00000116741,ENSG00000119508,ENSG00000120129,ENSG00000123358,ENSG00000124216,ENSG00000124762,ENSG00000125266,ENSG00000125740,ENSG00000125812,ENSG00000128016,ENSG00000128272,ENSG00000128590,ENSG00000128965,ENSG00000130766,ENSG00000133639,ENSG00000136826,ENSG00000138166,ENSG00000140961,ENSG00000141232,ENSG00000141582,ENSG00000142178,ENSG00000143878,ENSG00000148926,ENSG00000157557,ENSG00000158050,ENSG00000159388,ENSG00000162772,ENSG00000163273,ENSG00000163659,ENSG00000163874,ENSG00000164056,ENSG00000164920,ENSG00000168140,ENSG00000169242,ENSG00000172216,ENSG00000172602,ENSG00000173334,ENSG00000173530,ENSG00000175197,ENSG00000177873,ENSG00000179388,ENSG00000183508,ENSG00000183691,ENSG00000184545,ENSG00000184557,ENSG00000185022,ENSG00000197019,ENSG00000198576,ENSG00000204103,ENSG00000204388,ENSG00000204389,ENSG00000245848,ENSG00000256223
#> 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000051108,ENSG00000087074,ENSG00000105327,ENSG00000128272,ENSG00000128965,ENSG00000141682,ENSG00000153094,ENSG00000172216,ENSG00000175197,ENSG00000204389
#> 5 ENSG00000051108,ENSG00000059728,ENSG00000087074,ENSG00000100292,ENSG00000102554,ENSG00000103044,ENSG00000105327,ENSG00000113070,ENSG00000113369,ENSG00000113916,ENSG00000116741,ENSG00000119508,ENSG00000119630,ENSG00000120129,ENSG00000120738,ENSG00000122877,ENSG00000123358,ENSG00000124216,ENSG00000124762,ENSG00000125740,ENSG00000125812,ENSG00000128016,ENSG00000128272,ENSG00000128965,ENSG00000130766,ENSG00000133639,ENSG00000135625,ENSG00000136826,ENSG00000141232,ENSG00000141582,ENSG00000141682,ENSG00000142178,ENSG00000144655,ENSG00000148926,ENSG00000152433,ENSG00000153094,ENSG00000157557,ENSG00000159388,ENSG00000162772,ENSG00000162783,ENSG00000163659,ENSG00000163874,ENSG00000164056,ENSG00000164920,ENSG00000169242,ENSG00000170345,ENSG00000170684,ENSG00000171223,ENSG00000172216,ENSG00000173334,ENSG00000175197,ENSG00000177873,ENSG00000179388,ENSG00000183508,ENSG00000183691,ENSG00000184557,ENSG00000185022,ENSG00000197019,ENSG00000197566,ENSG00000204103,ENSG00000204388,ENSG00000204389,ENSG00000245848,ENSG00000256223,ENSG00000257446
#> 6                 ENSG00000051108,ENSG00000059728,ENSG00000087074,ENSG00000100292,ENSG00000102554,ENSG00000103044,ENSG00000105327,ENSG00000113070,ENSG00000113369,ENSG00000113916,ENSG00000116741,ENSG00000119508,ENSG00000119630,ENSG00000120129,ENSG00000120738,ENSG00000122877,ENSG00000123358,ENSG00000124216,ENSG00000124762,ENSG00000125740,ENSG00000125812,ENSG00000128016,ENSG00000128272,ENSG00000128965,ENSG00000130766,ENSG00000133639,ENSG00000135625,ENSG00000136826,ENSG00000141232,ENSG00000141582,ENSG00000141682,ENSG00000142178,ENSG00000144655,ENSG00000152433,ENSG00000153094,ENSG00000157557,ENSG00000159388,ENSG00000162772,ENSG00000162783,ENSG00000163659,ENSG00000163874,ENSG00000164056,ENSG00000164920,ENSG00000169242,ENSG00000170345,ENSG00000170684,ENSG00000171223,ENSG00000172216,ENSG00000173334,ENSG00000175197,ENSG00000177873,ENSG00000179388,ENSG00000183508,ENSG00000183691,ENSG00000184557,ENSG00000185022,ENSG00000197019,ENSG00000197566,ENSG00000204103,ENSG00000204388,ENSG00000204389,ENSG00000245848,ENSG00000256223,ENSG00000257446

## Create network for Gene Ontology - Molecular Function related results
## in Cytoscape (when the application is opened)
## Otherwise, create a CX file in the temporary directory
## The file can be opened in Cytoscape
createNetwork(gostObject=parentalNapaVsDMSOEnrichment, source="KEGG", 
    removeRoot=FALSE, title="KEGG Graph", 
    fileName=file.path(tempdir(), "KEGG_demo.cx"))
#> Unable to connect to Cytoscape. 
#> CX JSON file will be created.
#> Preparing information for generating CX JSON file.
#> CX JSON file "/tmp/Rtmp8gkWVv/KEGG_demo.cx" has been created.
#> [1] TRUE

```
