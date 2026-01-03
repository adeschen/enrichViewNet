# The result of a differential expression analysis done between napabucasin treated and DMSO control parental MiaPaCa2 cells. The cells were treated for 2 hours with 0.5 uM napabucasin. The protocol to generate the RNA-seq is described in Froeling F.E.M. et al 2019.

The object is a `data.frame` with 24184 rows and 4 columns. Each row
correspond to a tested gene.

## Usage

``` r
data(parentalNapaVsDMSODEG)
```

## Format

a `data.frame` containing the results of a differential expression
analysis between napabucasin treated and DMSO control parental MiaPaCa2
cells for all 24184 genes tested. The 4 columns are:

- `"EnsemblID"`: a `character` string representing the unique Ensembl
  identifier for the tested gene

- `"log2FoldChange"`: a `numeric` representing the expression difference
  (in log2FoldChange) between the napabucasin treatment and the DMSO
  control for the tested gene

- `"padj"`: a `numeric` representing the adjusted p-value associated to
  the difference in expression for the tested gene

- `"GeneName"`: a `character` string representing the name of the tested
  gene

## Source

The original RNA-sequencing data is available at the Gene Expression
Omnibus (GEO) under the accession number GSE135352.

## Value

a `data.frame` containing the results of a differential expression
analysis between napabucasin treated and DMSO control parental MiaPaCa2
cells for all 24184 genes tested. The 4 columns are:

- `"EnsemblID"`: a `character` string representing the unique Ensembl
  identifier for the tested gene

- `"log2FoldChange"`: a `numeric` representing the expression difference
  (in log2FoldChange) between the napabucasin treatment and the DMSO
  control for the tested gene

- `"padj"`: a `numeric` representing the adjusted p-value associated to
  the difference in expression for the tested gene

- `"GeneName"`: a `character` string representing the name of the tested
  gene

## Details

The differentially expressed genes between napabucasin-treated cells
(0.5 uM) and DMSO as vehicle control are reprinted from Clinical Cancer
Research, 2019, 25 (23), 7162–7174, Fieke E.M. Froeling, Manojit Mosur
Swamynathan, Astrid Deschênes, Iok In Christine Chio, Erin Brosnan,
Melissa A. Yao, Priya Alagesan, Matthew Lucito, Juying Li, An-Yun Chang,
Lloyd C. Trotman, Pascal Belleau, Youngkyu Park, Harry A. Rogoff, James
D. Watson, David A. Tuveson, Bioactivation of napabucasin triggers
reactive oxygen species–mediated cancer cell death, with permission from
AACR.

## See also

- [createNetwork](createNetwork.md) for transforming functional
  enrichment results from gprofiler2 into a Cytoscape network

- [createEnrichMap](createEnrichMap.md) for transforming functional
  enrichment results from gprofiler2 into an enrichment map

## Examples

``` r

## Required library
library(gprofiler2)

## Loading data set containing the results of a differentially expressed 
## analysis between 2-hour treatment with 0.5 uM napabucasin and 
## DMSO vehicle control parental MiaPaCa2 cells
data(parentalNapaVsDMSODEG)

allGenes <- unique(parentalNapaVsDMSODEG$EnsemblID)

## Select the significantly differentially expressed genes
selection <- which(abs(parentalNapaVsDMSODEG$log2FoldChange) > 1 & 
                            parentalNapaVsDMSODEG$padj < 0.05)
                            
selectedGenes <- unique(parentalNapaVsDMSODEG$EnsemblID[selection])

## Run an enrichment analysis using WikiPathways dataset
gostres <- gost(query = list(parental_napa_vs_DMSO=selectedGenes),
    organism="hsapiens",
    correction_method = "g_SCS",
    sources=c("WP"), significant=TRUE, evcodes=TRUE,
    custom_bg=allGenes, exclude_iea=TRUE)
#> Detected custom background input, domain scope is set to 'custom'.


```
