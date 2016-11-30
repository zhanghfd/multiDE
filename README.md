# multiDE

multiDE: A dimension reduced model based statistical method for differential expression analysis using RNA-sequencing data from multiple conditions

Based on a dimension reduced ANOVA model, the R package ‘multiDE’ conducts differential expression (DE) analysis using high throughput next-generation RNA-seq read count data generated from samples of multiple conditions. The samples across all conditions can be either matched or unmatched. The package provides log2-fold change estimates and their standard errors, Wald test p-value of DE analysis for each gene.

The manual file is "multiDE-manual.pdf". 

Installation of multiDE in R:

> library(‘devtools’);

> install_github(‘zhanghfd/multiDE’);
