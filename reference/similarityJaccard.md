# Calculate the Jaccard similarity coefficient for the terms present in a data frame

The function calculates the Jaccard similarity coefficient for the
enriched terms present i a data frame.

## Usage

``` r
similarityJaccard(resultDF)
```

## Arguments

- resultDF:

  a `data.frame` containing the enrichment terms that are going to be
  graphed as an enrichment map. The `data.frame` must contain the
  `geneID` and `Description` columns.

## Value

a `matrix` with the lower triangle section containing the Jaccard
similarity coefficients as `numeric`. The upper triangle section is
filled with `NA` as well as the diagonal.

## Author

Astrid Deschênes

## Examples

``` r

## Data frame with the enriched terms
resultData <- data.frame( 
    ID=c("REAC:R-HSA-9031628", "REAC:R-HSA-9648895", "REAC:R-HSA-9614085",
            "REAC:R-HSA-166520"),
    Description=c("NGF-stimulated transcription", 
            "Response of EIF2AK1 (HRI) to heme deficiency",
            "FOXO-mediated transcription", "Signaling by NTRKs"),
    geneID=c(paste0("ENSG00000120738/ENSG00000122877/ENSG00000125740/",
        "ENSG00000135625/ENSG00000170345/ENSG00000171223/ENSG00000173334/", 
        "ENSG00000179388/ENSG0000019857"),
        paste0("ENSG00000087074/ENSG00000128272/ENSG00000128965/", 
        "ENSG00000162772/ENSG00000172216/ENSG00000175197"),
        paste0("ENSG00000105327/ENSG00000113916/ENSG00000124762/", 
        "ENSG00000133639/ENSG00000136826/ENSG00000153094/ENSG00000175197"),
        paste0("ENSG00000105327/ENSG00000113916/ENSG00000124762/", 
        "ENSG00000133639/ENSG00000136826/ENSG00000153094/ENSG0000017519/",
        "ENSG00000120738/ENSG00000122877/ENSG00000125740/ENSG00000135625/",
        "ENSG00000170345/ENSG00000171223/ENSG00000173334/ENSG00000179388/", 
        "ENSG00000198576")))

## Calculate the Jaccard similarity coefficient for all pairwise 
## term combination
enrichViewNet:::similarityJaccard(resultDF=resultData)
#>                                              NGF-stimulated transcription
#> NGF-stimulated transcription                                           NA
#> Response of EIF2AK1 (HRI) to heme deficiency                    0.0000000
#> FOXO-mediated transcription                                     0.0000000
#> Signaling by NTRKs                                              0.4705882
#>                                              Response of EIF2AK1 (HRI) to heme deficiency
#> NGF-stimulated transcription                                                           NA
#> Response of EIF2AK1 (HRI) to heme deficiency                                           NA
#> FOXO-mediated transcription                                                    0.08333333
#> Signaling by NTRKs                                                             0.00000000
#>                                              FOXO-mediated transcription
#> NGF-stimulated transcription                                          NA
#> Response of EIF2AK1 (HRI) to heme deficiency                          NA
#> FOXO-mediated transcription                                           NA
#> Signaling by NTRKs                                             0.3529412
#>                                              Signaling by NTRKs
#> NGF-stimulated transcription                                 NA
#> Response of EIF2AK1 (HRI) to heme deficiency                 NA
#> FOXO-mediated transcription                                  NA
#> Signaling by NTRKs                                           NA
    
```
