# evoGO R package


*evoGO* is an R package developed by Evotec International GmbH that provides an 
advanced functionality for performing a Gene Ontology (GO) enrichment analysis. 
The *evoGO* algorithm deprioritizes GO terms that are redundant or unspecific within 
the current context making interpretation of results less time-consuming and more 
efficient.


## Method

Redundancy is a known issue with results of the classic GO overrepresentation analysis. 
It originates from a significant overlap between gene sets of parent and child GO terms 
and makes it harder to pinpoint the most essential findings. We developed an algorithm 
that deprioritizes redundant terms by taking GO graph topology into account. 

Briefly, the idea behind the *evoGO* algorithm is to reduce the significance of less 
specific terms if their differentially expressed genes (DEGs) contribute to higher 
significance of more specific descendant terms.
To achieve that:
- The algorithm down-weights genes attributed to a term in all its less significant 
ancestors;
- Gene down-weighting is proportional to the ratio of ancestor and descendant term 
p values;
- Sums of DEG weights instead of DEG counts are used to re-calculate p values with 
Fisher’s exact test.

The *evoGO* code is specifically optimized to allow performing batch 
calculations involving of thousands of analyses (e.g., in drug screening studies) within a 
reasonable time frame. Acknowledging the importance of result reproducibility, the *evoGO* provides
functionality for extracting and using specific versions of Gene Ontology and Ensembl gene 
annotation. Please refer to the documentation for details.


## Dependencies

- assertthat
- biomaRt
- curl
- dplyr
- ontologyIndex
- parallel


## Installation
devtools::install_github("Evotec-Bioinformatics/evoGO")


## Quick start

Before starting your work with *evoGO* package it is recommended to acquire an up-to-date 
GO annotation for species of interest. *evoGO* package has built-in means for downloading 
latest version of Gene Ontology from [GO Consortium website](http://geneontology.org/) and 
GO term annotation from [Ensembl database](https://www.ensembl.org/), which together are used
to construct the GO annotation. To download the annotation use:

```r
goAnnotation <- getGOAnnotation("hsapiens")
```

Please note that downloading data from Ensembl server may take a while. We recommend 
to have the latest version of `biomaRt` package installed for optimal performance. 
For a quick test of package functions you can also download an example H. sapiens annotation 
file using `exampleGOAnnotation()`. This annotation is provided only for testing purposes 
and not updated regularly.

By default, the annotation files are stored in the `extdata` directory located inside 
*evoGO* installation directory, which can be found using `path.package("evoGO")` after 
the package was loaded. A list of the stored annotations can be viewed:

```r
listGOAnnotations()
```

You can quickly load latest previously stored annotation:

```r
goAnnotation <- loadGOAnnotation("hsapiens")
```

In the following example of GO enrichment analysis, we use one of the sample datasets 
included with the package. The data was obtained by performing differential expression 
analysis for heart and brain human tissue samples available from the 
[GTEx V8 database](https://gtexportal.org/) (see the package documentation for details).

```r
# Get IDs of differentially expressed genes
degs <- gtex_heart_brain$GeneID[gtex_heart_brain$Pvalue.adjusted < 0.05]
# Get IDs of all available genes
universe <- gtex_heart_brain$GeneID
# Perform GO enrichment analysis
result <- calcGOenrichment(goAnnotation, deGenes = degs, domain = "BP", universe = universe)
head(result, 10) 

#         id                                   name                                    def fisher.pvalue evogo.pvalue annotated significant
# GO:0048731                     system development "The process whose specific outcome...  4.908115e-35 4.908115e-35      3936        3295
# GO:0007399             nervous system development "The process whose specific outcome...  3.320425e-32 3.320425e-32      2398        2054
# GO:0048699                  generation of neurons "The process in which nerve cells a...  1.022475e-25 1.022475e-25      1434        1250
# GO:0065008       regulation of biological quality "Any process that modulates a quali...  3.645834e-25 3.645834e-25      3188        2659
# GO:0051179                           localization "Any process in which a cell, a sub...  9.227587e-25 9.227587e-25      3799        3141
# GO:0023051                regulation of signaling "Any process that modulates the fre...  1.096125e-23 1.096125e-23      2990        2495
# GO:0030182                 neuron differentiation "The process in which a relatively ...  2.999252e-23 2.999252e-23      1366        1187
# GO:0010646       regulation of cell communication "Any process that modulates the fre...  1.136421e-22 1.136421e-22      2948        2457
# GO:0051128 regulation of cellular component or... "Any process that modulates the fre...  2.166147e-22 2.166147e-22      2107        1784
# GO:0030030           cell projection organization "A process that is carried out at t...  1.136513e-20 1.136513e-20      1504        1291
```

As you can see, top 10 enriched terms contain *"nervous system development"*, 
*"generation of neurons"*, *"neuron differentiation"*, which reflects the 
fundamental differences between brain and heart tissue.

You can also find that more than 300 of redundant terms no longer appear as significantly 
enriched when *evoGO* algorithm is applied:

```r
sum(result$fisher.pvalue < 0.05)
# 1520
sum(result$evogo.pvalue < 0.05)
# 1209
```


## License

GPL-2


## References

Ashburner et al. Gene ontology: tool for the unification of biology. Nat Genet. May 2000;25(1):25-9. [[abstract](https://www.ncbi.nlm.nih.gov/pubmed/10802651) | [full text](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3037419/)]

The Gene Ontology resource: enriching a GOld mine. Nucleic Acids Res. Jan 2021;49(D1):D325-D334. [[abstract](https://pubmed.ncbi.nlm.nih.gov/33290552/) | [full text](https://academic.oup.com/nar/article-pdf/49/D1/D325/35364517/gkaa1113.pdf)]

