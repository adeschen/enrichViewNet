<!-- badges: start -->
[![R-CMD-check-bioc](https://github.com/adeschen/enrichViewNet/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/adeschen/enrichViewNet/actions/workflows/check-bioc.yml)
[![Codecov test coverage](https://codecov.io/gh/adeschen/enrichViewNet/branch/main/graph/badge.svg)](https://codecov.io/gh/adeschen/enrichViewNet?branch=main)
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
**_enrichplot_** package, enable the visualization of enriched terms 
into a network with edges connecting overlapping genes. Thus, enriched terms 
with overlapping genes cluster together.


## Authors ##

[Astrid Desch&ecirc;nes](http://ca.linkedin.com/in/astriddeschenes "Astrid Desch&ecirc;nes"), 
[Pascal Belleau](http://ca.linkedin.com/in/pascalbelleau "Pascal Belleau"), 
[Robert L Faure](https://www.crchudequebec.ulaval.ca/en/research/researchers/robert-l-faure/), 
[Maria J Fernandes](https://www.crchudequebec.ulaval.ca/en/research/researchers/maria-fernandes/) and 
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


## Bugs/Feature requests ##

If you have any bugs or feature requests, 
[let us know](https://github.com/adeschen/enrichViewNet/issues). 

Thanks!
