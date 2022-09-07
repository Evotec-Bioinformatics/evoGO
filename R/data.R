#' Differential gene expression: Heart (Left Ventricle) vs. Lung
#'
#' Human tissue gene expression data was obtained from GTEx V8 database. Only 
#' samples with tissue degradation level 0 or 1 and death circumstances 1-3 on 
#' 4-point Hardy Scale were used.
#' In total, 29 heart samples and 14 lung samples were used for the analysis.
#'
#' @format A data frame with 34638 rows and 5 variables:
#' \describe{
#'   \item{GeneID}{Ensembl gene ID}
#'   \item{CPM}{counts per million, average of the normalized expression values}
#'   \item{log2FoldChange}{log2(FoldChange), effect size estimate}
#'   \item{Pvalue}{P value}
#'   \item{Pvalue.adjusted}{P value adjusted using Benjamini-Hochberg procedure}
#' }
#' @source \url{https://gtexportal.org/}
"gtex_heart_lung"

#' Differential gene expression: Heart (Left Ventricle) vs. Brain (Cortex)
#'
#' Human tissue gene expression data was obtained from GTEx V8 database. Only 
#' samples with tissue degradation level 0 or 1 and death circumstances 1-3 on 
#' 4-point Hardy Scale were used.
#' In total, 29 heart samples and 19 brain samples were used for the analysis.
#'
#' @format A data frame with 36697 rows and 5 variables:
#' \describe{
#'   \item{GeneID}{Ensembl gene ID}
#'   \item{CPM}{counts per million, average of the normalized expression values}
#'   \item{log2FoldChange}{log2(FoldChange), effect size estimate}
#'   \item{Pvalue}{P value}
#'   \item{Pvalue.adjusted}{P value adjusted using Benjamini-Hochberg procedure}
#' }
#' @source \url{https://gtexportal.org/}
"gtex_heart_brain"

#' Differential gene expression: Lung vs. Brain (Cortex)
#'
#' Human tissue gene expression data was obtained from GTEx V8 database. Only 
#' samples with tissue degradation level 0 or 1 and death circumstances 1-3 on 
#' 4-point Hardy Scale were used.
#' In total, 14 lung samples and 19 brain samples were used for the analysis.
#'
#' @format A data frame with 36758 rows and 5 variables:
#' \describe{
#'   \item{GeneID}{Ensembl gene ID}
#'   \item{CPM}{counts per million, average of the normalized expression values}
#'   \item{log2FoldChange}{log2(FoldChange), effect size estimate}
#'   \item{Pvalue}{P value}
#'   \item{Pvalue.adjusted}{P value adjusted using Benjamini-Hochberg procedure}
#' }
#' @source \url{https://gtexportal.org/}
"gtex_lung_brain"