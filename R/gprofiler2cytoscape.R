#' From functional enrichment results to biological networks
#'
#' The \code{gprofiler2cytoscape} package enables the visualization 
#' of enrichment results obtained by \code{gprofiler2} 
#' (\url{https://cran.r-project.org/web/packages/gprofiler2/index.html})  
#' under the form of \code{Cytoscape} network (\url{https://cytoscape.org/}). 
#' 
#' In those networks, both gene datasets (GO terms/pathways/protein complexes)
#' and genes are represented as nodes. A edge connect a gene to its datasets.
#' In the current version, only genes present in at least one gene dataset 
#' are retained.
#'
#' @docType package
#'
#' @name gprofiler2cytoscape-package
#'
#' @aliases gprofiler2cytoscape-package gprofiler2cytoscape
#'
#' @author Astrid Deschênes and
#' Pascal Belleau
#'
#' Maintainer:
#' Astrid Deschênes <adeschen@hotmail.com>
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{createNetwork}} {for transforming functional 
#'     enrichment results from gprofiler2 into a Cytoscape network}
#' }
#' 
#' @encoding UTF-8
#' @keywords package
NULL


#' The result of a functional enrichment analysis done with \code{gprofiler2} 
#' (\url{https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html}).
#' 
#' The object is a \code{list} with 2 entries. It contains the results of the
#' enrichment analysis as well as the metadata related to the analysis.
#'
#' This dataset can be
#' used to test the \code{\link{createNetwork}} function.
#'
#' @name demoGOST
#'
#' @docType data
#'
#' @aliases demoGOST
#'
#' @format The \code{list} contains two entries. The \code{result} entry 
#' contains a \code{data.frame} with the significant results obtained by
#' an enrichment analysis done with \code{gprofiler2}. The \code{meta} entry 
#' contains a named list with all the metadata for the query.
#'
#' @return  A \code{list} containing two entries. The \code{result} entry 
#' contains a \code{data.frame} with the significant results obtained by
#' an enrichment analysis done with \code{gprofiler2}. 
#' The \code{meta} entry contains a named list with all the 
#' metadata for the query.
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{createNetwork}} {for transforming functional 
#'     enrichment results from gprofiler2 into a Cytoscape network}
#' }
#'
#' @usage data(demoGOST)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading dataset containing result from an enrichment analysis done with
#' ## gprofiler2
#' data(demoGOST)
#'
#' \dontrun{
#' 
#' ## Create network for WikiPathways results
#' createNetwork(gostObject = demoGOST, source="WP", title="Wikipathways")
#' 
#' }
#'
#'
NULL


#' The result of a differential expression analysis done between 
#' napabucasin treated and DMSO control parental 
#' MiaPaCa2 cells. The cells were treated for 2 hour with 0.5 uM napabucasin.
#' The protocol to generate the RNA-seq is described 
#' in Froeling F.E.M. et al 2019.
#' 
#' The object is a \code{data.frame} with 24184 rows and 4 columns. 
#' Each row correspond to a tested gene.
#'
#' @name parentalNapaVsDMSODEG
#'
#' @docType data
#'
#' @aliases parentalNapaVsDMSODEG
#'
#' @format a \code{data.frame} containing the results of a differential 
#' expression analysis between napabucasin treated and DMSO control parental 
#' MiaPaCa2 cells for all 24184 genes tested. The 4 columns are:
#' \itemize{
#' \item{EnsemblID} {a \code{character} string representing the unique Ensembl 
#' identifier for the tested gene}
#' \item{EnsemblID} {a \code{numeric} representing the expression difference 
#' (in log2FoldChange) between the napabucasin treatment and the DMSO control 
#' for the tested gene}
#' \item{padj} {a \code{numeric} representing the adjusted p-value associated  
#' to the difference in expression for the tested gene; \code{NA} when the 
#' adjusted p-value as not been calculated (equivalent to not significant)}
#' \item{GeneName} {a \code{character} string representing the name of 
#' the tested gene}
#' }
#'
#' @return  a \code{data.frame} containing the results of a differential 
#' expression analysis between napabucasin treated and DMSO control parental 
#' MiaPaCa2 cells for all 24184 genes tested. The 4 columns are:
#' \itemize{
#' \item{EnsemblID} {a \code{character} string representing the unique Ensembl 
#' identifier for the tested gene}
#' \item{EnsemblID} {a \code{numeric} representing the expression difference 
#' (in log2FoldChange) between the napabucasin treatment and the DMSO control 
#' for the tested gene}
#' \item{padj} {a \code{numeric} representing the adjusted p-value associated  
#' to the difference in expression for the tested gene; \code{NA} when the 
#' adjusted p-value as not been calculated (equivalent to not significant)}
#' \item{GeneName} {a \code{character} string representing the name of 
#' the tested gene}
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{createNetwork}} {for transforming functional 
#'     enrichment results from gprofiler2 into a Cytoscape network}
#' }
#'
#' @usage data(parentalNapaVsDMSODEG)
#'
#' @keywords datasets
#'
#' @details
#' 
#' The differentially expressed genes between napabucasin-treated 
#' cells (0.5 uM) and DMSO as vehicle control are reprinted from Clinical 
#' Cancer Research, 2019, 25 (23), 7162–7174, Fieke E.M. Froeling, Manojit 
#' Mosur Swamynathan, Astrid Deschênes, Iok In Christine Chio, Erin Brosnan, 
#' Melissa A. Yao, Priya Alagesan, Matthew Lucito, Juying Li, An-Yun Chang, 
#' Lloyd C. Trotman, Pascal Belleau, Youngkyu Park, Harry A. Rogoff, 
#' James D. Watson, David A. Tuveson, Bioactivation of napabucasin triggers 
#' reactive oxygen species–mediated cancer cell death, with permission 
#' from AACR.
#' 
#' @source 
#' 
#' The original RNA-sequencing data is available at the Gene Expression 
#' Omnibus (GEO) under the accession number GSE135352.
#' 
#' @examples
#'
#' ## Required library
#' library(gprofiler2)
#' 
#' ## Loading data set containing the results of a differentially expressed 
#' ## analysis between 2-hour treatment with 0.5 uM napabucasin and 
#' ## DMSO vehicle control parental MiaPaCa2 cells
#' data(parentalNapaVsDMSODEG)
#' 
#' allGenes <- unique(parentalNapaVsDMSODEG$EnsemblID)
#' 
#' ## Select the significantly differentially expressed genes
#' selection <- which(abs(parentalNapaVsDMSODEG$log2FoldChange) > 1 & 
#'                             parentalNapaVsDMSODEG$padj < 0.05)
#'                             
#' selectedGenes <- unique(parentalNapaVsDMSODEG$EnsemblID[selection])
#' 
#' ## Run an enrichment analysis using WikiPathways dataset
#' gostres <- gost(query = list(parental_napa_vs_DMSO=selectedGenes),
#'     organism="hsapiens",
#'     correction_method = "g_SCS",
#'     sources=c("WP"), significant=TRUE, evcodes=TRUE,
#'     custom_bg=allGenes, exclude_iea=TRUE)
#' 
#' 
#'
NULL


#' The result of a differential expression analysis done between 
#' napabucasin treated and DMSO control MiaPaCa2 cells stably expressing 
#' the Rosa26 control vector. The cells were treated for 2 hour 
#' with 0.5 uM napabucasin.
#' The protocol to generate the RNA-seq is described 
#' in Froeling F.E.M. et al 2019.
#' 
#' The object is a \code{data.frame} with 23542 rows and 4 columns. 
#' Each row correspond to a tested gene.
#'
#' @name rosaNapaVsDMSODEG
#'
#' @docType data
#'
#' @aliases rosaNapaVsDMSODEG
#'
#' @format a \code{data.frame} containing the results of a differential 
#' expression analysis between napabucasin treated and DMSO control MiaPaCa2 
#' cells stably expressing the Rosa26 control vector for all 23542 genes 
#' tested. The 4 columns are:
#' \itemize{
#' \item{EnsemblID} {a \code{character} string representing the unique Ensembl 
#' identifier for the tested gene}
#' \item{EnsemblID} {a \code{numeric} representing the expression difference 
#' (in log2FoldChange) between the napabucasin treatment and the DMSO control 
#' for the tested gene}
#' \item{padj} {a \code{numeric} representing the adjusted p-value associated  
#' to the difference in expression for the tested gene}
#' \item{GeneName} {a \code{character} string representing the name of 
#' the tested gene}
#' }
#'
#' @return  a \code{data.frame} containing the results of a differential 
#' expression analysis between napabucasin treated and DMSO control MiaPaCa2 
#' cells stably expressing the Rosa26 control vector for all 23542 genes 
#' tested. The 4 columns are:
#' \itemize{
#' \item{EnsemblID} {a \code{character} string representing the unique Ensembl 
#' identifier for the tested gene}
#' \item{log2FoldChange} {a \code{numeric} representing the expression 
#' difference (in log2FoldChange) between the napabucasin treatment and 
#' the DMSO control for the tested gene}
#' \item{padj} {a \code{numeric} representing the adjusted p-value associated  
#' to the difference in expression for the tested gene}
#' \item{GeneName} {a \code{character} string representing the name of 
#' the tested gene}
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{createNetwork}} {for transforming functional 
#'     enrichment results from gprofiler2 into a Cytoscape network}
#' }
#'
#' @usage data(rosaNapaVsDMSODEG)
#'
#' @keywords datasets
#'
#' @details
#' 
#' The differentially expressed genes between napabucasin-treated 
#' cells (0.5 uM) and DMSO as vehicle control are reprinted from Clinical 
#' Cancer Research, 2019, 25 (23), 7162–7174, Fieke E.M. Froeling, Manojit 
#' Mosur Swamynathan, Astrid Deschênes, Iok In Christine Chio, Erin Brosnan, 
#' Melissa A. Yao, Priya Alagesan, Matthew Lucito, Juying Li, An-Yun Chang, 
#' Lloyd C. Trotman, Pascal Belleau, Youngkyu Park, Harry A. Rogoff, 
#' James D. Watson, David A. Tuveson, Bioactivation of napabucasin triggers 
#' reactive oxygen species–mediated cancer cell death, with permission 
#' from AACR.
#' 
#' @source 
#' 
#' The original RNA-sequencing data is available at the Gene Expression 
#' Omnibus (GEO) under the accession number GSE135352.
#' 
#' @examples
#'
#' ## Required library
#' library(gprofiler2)
#' 
#' ## Loading dataset containing the results of a differentially expressed 
#' ## analysis between 2-hour treatment with 0.5 uM napabucasin and 
#' ## DMSO vehicle control MiaPaCa2 cells stably expressing the 
#' ## Rosa26 control vector
#' data(rosaNapaVsDMSODEG)
#' 
#' allGenes <- unique(rosaNapaVsDMSODEG$EnsemblID)
#' 
#' ## Select the significantly differentially expressed genes
#' selection <- which(abs(rosaNapaVsDMSODEG$log2FoldChange) > 1 & 
#'                             rosaNapaVsDMSODEG$padj < 0.05)
#'                             
#' selectedGenes <- unique(rosaNapaVsDMSODEG$EnsemblID[selection])
#' 
#' ## Run an enrichment analysis using Transfac dataset (transcription factor)
#' gostres <- gost(query = list(rosa_napa_vs_DMSO=selectedGenes),
#'     organism="hsapiens",
#'     correction_method = "g_SCS",
#'     sources=c("TF"), significant=TRUE, evcodes=TRUE,
#'     custom_bg=allGenes, exclude_iea=TRUE)
#' 
#' 
#'
NULL


#' The result of an enrichment analysis has been done using the significantly 
#' differentially expressed genes between napabucasin treated and DMSO 
#' control parental MiaPaCa2 cells. 
#' The cells were treated for 2 hour with 0.5 uM napabucasin.  
#' The protocol to generate the RNA-seq is described 
#' in Froeling F.E.M. et al 2019.
#' 
#' The enrichment analysis was done with gprofile2 package 
#' (Kolberg L et al 2020) with database version 'e109_eg56_p17_1d3191d' and 
#' g:SCS multiple testing correction method applying significance 
#' threshold of 0.05 (Raudvere U et al 2019). All tested genes were used 
#' as background.
#' 
#' The object is a named \code{list} with 2 entries. The 'result' entry 
#' contains a \code{data.frame} with the enrichment analysis results and 
#' the 'meta' entry contains metadata information.
#' 
#' @name parentalNapaVsDMSOEnrichment
#'
#' @docType data
#'
#' @aliases parentalNapaVsDMSOEnrichment
#'
#' @format a \code{list} containing 2 entries:
#' \itemize{
#' \item{result} {a \code{data.frame} with the significantly enriched 
#' terms }
#' \item{meta} {a TODO}
#' }
#'
#' @return  a \code{list} containing 2 entries:
#' \itemize{
#' \item{result} {a \code{data.frame} with the significantly enriched 
#' terms }
#' \item{meta} {a TODO}
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{createNetwork}} {for transforming functional 
#'     enrichment results from gprofiler2 into a Cytoscape network}
#' }
#'
#' @usage data(parentalNapaVsDMSOEnrichment)
#'
#' @keywords datasets
#'
#' @details
#' 
#' The dataset used for the enrichment analysis is 
#' associated to this publication:
#' 
#' Froeling F.E.M. et al.Bioactivation of Napabucasin Triggers Reactive Oxygen 
#' Species–Mediated Cancer Cell Death. Clin Cancer Res 
#' 1 December 2019; 25 (23): 7162–7174
#' 
#' The enrichment analysis has been done with gprofile2 package 
#' (Kolberg L et al 2020) with database version 'e109_eg56_p17_1d3191d' and 
#' g:SCS multiple testing correction method applying significance 
#' threshold of 0.05 (Raudvere U et al 2019). All tested genes were used 
#' as background.
#' 
#' @source 
#' 
#' The original RNA-sequencing data is available at the Gene Expression 
#' Omnibus (GEO) under the accession number GSE135352.
#' 
#' @examples
#'
#' ## Required library
#' library(gprofiler2)
#' 
#' ## Loading dataset containing the results of a differentially expressed 
#' ## analysis between 2-hour treatment with 0.5 uM napabucasin and 
#' ## DMSO vehicle control parental MiaPaCa2 cells
#' data(parentalNapaVsDMSOEnrichment)
#' 
#' ## TODO
#' 
#'
NULL