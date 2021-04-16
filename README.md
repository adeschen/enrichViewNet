<!-- badges: start -->
[![R-CMD-check-bioc](https://github.com/adeschen/gprofiler2cytoscape/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/adeschen/gprofiler2cytoscape/actions/workflows/check-bioc.yml)
[![Codecov test coverage](https://codecov.io/gh/adeschen/gprofiler2cytoscape/branch/main/graph/badge.svg)](https://codecov.io/gh/adeschen/gprofiler2cytoscape?branch=main)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic%202.0-0298c3.svg)](https://opensource.org/licenses/Artistic-2.0)
<!-- badges: end -->


The **_gprofiler2cytoscape_** package enables the transformation of 
functional enrichment results obtained by [gprofiler2](https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html) 
into a 
[Cytoscape](https://cytoscape.org/) network where significant terms (Gene Ontology/Reactome/Kegg/Transcription Factor/etc..) and 
genes are nodes. In this network, edges connect genes that are part of the
enrichment request to their terms. The **_gprofiler2cytoscape_** package 
offers the option to generate a network for only a portion of the 
significant terms by selecting the source or by providing a 
specific list of terms.


## Authors ##

[Astrid Desch&ecirc;nes](http://ca.linkedin.com/in/astriddeschenes "Astrid Desch&ecirc;nes") and
[Pascal Belleau](http://ca.linkedin.com/in/pascalbelleau "Pascal Belleau")


## License ##

This package and the underlying **_gprofiler2cytoscape_** code are distributed under 
the Artistic license 2.0. You are free to use and redistribute this software. 

For more information on Artistic 2.0 License see
[http://opensource.org/licenses/Artistic-2.0](http://opensource.org/licenses/Artistic-2.0)


## Documentation ##

[gprofiler2cytoscape Website](https://adeschen.github.io/gprofiler2cytoscape/)

[gprofiler2cytoscape Get Started](https://adeschen.github.io/gprofiler2cytoscape/articles/gprofiler2cytoscape.html)


## Installation ##

To install the latest version accessible, the  [devtools](https://cran.r-project.org/web/packages/devtools/index.html) 
package is required.

     ## Load required package
     library(devtools)

     ## Install the latest version of gprofiler2cytoscape
     devtools::install_github('adeschen/gprofiler2cytoscape')


## Bugs/Feature requests ##

If you have any bugs or feature requests, 
[let us know](https://github.com/adeschen/gprofiler2cytoscape/issues). 

Thanks!
