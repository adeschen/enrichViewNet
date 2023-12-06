<!-- badges: start -->
[![R-CMD-check-bioc](https://github.com/adeschen/enrichViewNet/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/adeschen/enrichViewNet/actions/workflows/check-bioc.yml)
[![codecov](https://codecov.io/gh/adeschen/enrichViewNet/graph/badge.svg?token=N3RA2934V5)](https://codecov.io/gh/adeschen/enrichViewNet)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic%202.0-0298c3.svg)](https://opensource.org/licenses/Artistic-2.0)
<!-- badges: end -->


The **enrichViewNet** package enables the transformation of 
functional enrichment results, formatted as the results obtained  by [gprofiler2](https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html), 
into networks. 

First, the **enrichViewNet** package enables the 
visualization of enrichment results a [Cytoscape](https://cytoscape.org/) 
network where significant terms 
(Gene Ontology/Reactome/Transcription Factor/etc..) and 
genes are nodes. In this network, edges connect genes that are part of the
enrichment request to their terms. 

In addition, the **enrichViewNet** package also provides the option to 
create enrichment maps from functional enrichment results. 
Enrichment maps, as introduced in the Bioconductor 
[enrichplot](https://bioconductor.org/packages/release/bioc/html/enrichplot.html) package, 
enable the visualization of enriched terms 
into a network with edges connecting overlapping genes. Thus, enriched terms 
with overlapping genes cluster together.

## Bioconductor Package ##

[![Bioconductor Time](https://bioconductor.org/shields/years-in-bioc/enrichViewNet.svg)](https://bioconductor.org/packages/enrichViewNet)

The **enrichViewNet** package is now an official package 
of [Bioconductor](http://bioconductor.org/). 

The current Bioconductor release can be directly downloaded from their website:
[Current release](https://bioconductor.org/packages/enrichViewNet/)


## Authors ##

[Astrid Desch&ecirc;nes](http://ca.linkedin.com/in/astriddeschenes "Astrid Desch&ecirc;nes"), 
[Pascal Belleau](http://ca.linkedin.com/in/pascalbelleau "Pascal Belleau"), 
[Robert L Faure](https://www.crchudequebec.ulaval.ca/en/research/researchers/robert-l-faure/), 
[Maria J Fernandes](https://www.crchudequebec.ulaval.ca/en/research/researchers/maria-fernandes/),
[Alexander Krasnitz](https://www.cshl.edu/research/faculty-staff/alexander-krasnitz/) and 
[David A Tuveson](https://tuvesonlab.labsites.cshl.edu/)

## License ##

This package and the underlying **enrichViewNet** code are distributed under 
the Artistic license 2.0. You are free to use and redistribute this software. 

For more information on Artistic 2.0 License see
[http://opensource.org/licenses/Artistic-2.0](http://opensource.org/licenses/Artistic-2.0)


## Documentation ##

[enrichViewNet Website](https://adeschen.github.io/enrichViewNet/)

[enrichViewNet Get Started](https://adeschen.github.io/enrichViewNet/articles/enrichViewNet.html)


## Installation ##

To install the latest version accessible, the 
[devtools](https://cran.r-project.org/web/packages/devtools/index.html) 
package is required.

     ## Load required package
     library(devtools)

     ## Install the latest version of enrichViewNet
     devtools::install_github('adeschen/enrichViewNet')


To install this package 
from [Bioconductor](https://bioconductor.org), start R 
(version 4.3 or later) and enter: 

     if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

     BiocManager::install("enrichViewNet")


## Bugs/Feature requests ##

If you have any bugs or feature requests, 
[let us know](https://github.com/adeschen/enrichViewNet/issues). 

Thanks!
