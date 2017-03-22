This pipeline is to relate tumor samples and cell lines using gene expression data ([Chen B. BMC medical genomics, 2015](https://bmcmedgenomics.biomedcentral.com/articles/10.1186/1755-8794-8-S2-S5)). Tumors samples that are not correlated to cell lines should be removed in creating disease gene expression signatures. Cell lines that are correlated to tumor samples are suggested to used in validation.

1. Correlate tumor samples from TCGA and cell lines from CCLE using their gene expression profiles
2. Compute HCC gene expression signatures

run workflow.R to reproduce results. We have provided the HCC signatures in the data file, you don't have to run this component.


