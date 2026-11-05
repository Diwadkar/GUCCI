# GUCCI
Genetic determinants Underlying Cell-Cell Interaction (GUCCI) model to identify cell-cell interaction QTLs (cci-QTLs)

This repository provides an implementation of the GUCCI framework for identifying cell-cell interaction quantitative trait loci (cci-QTLs) by integrating single-cell RNA sequencing (scRNA-seq) and bulk eQTL summary statistics. GUCCI leverages information across cell types to refine genetic effect size estimates and identify genetic variants whose regulatory effects on receptor gene expression depend on the expression levels of their corresponding ligands.

The repository includes functions to estimate cell-type weights using non-negative least squares (NNLS) and to implement the GUCCI model for jointly refining genetic effect sizes (β) and their standard errors (SEs), enabling the identification of cci-QTLs.

## Run GUCCI

### Step 1: Cell type weights

To estimate cell type proportion weights, we use a non-negative least squares (NNLS) approach. Users need to provide a vector of beta effect sizes (beta) and their corresponding standard error estimates (se) as input. The bulk data estimates should be the first entry in each input vector. Provide estimates across all SNPs tested genome-wide.

**Parameters:**
- `beta`: A numeric vector of length k+1 representing the beta effect size estimates of bulk tissue and k cell types.
- `se`: A numeric vector of length k+1 representing the standard error estimates of bulk tissue and k cell types.

```r
weights <- get_weights(beta, se)
```

### Step 2: Refine effect size estimates and identify cci-QTLs

The GUCCI model requires genotype–expression association summary statistics obtained from linear regression, with receptor gene expression as the outcome and genotype as the predictor, adjusting for relevant covariates.

Association analyses should be performed separately after stratifying individuals based on the expression levels of the corresponding ligand gene. Users can specify a single ligand expression threshold or evaluate multiple thresholds, as implemented in the published manuscript using percentiles ranging from the 10th to the 90th percentile.

GUCCI integrates bulk and cell-type-specific eQTL summary statistics to refine genetic effect size estimates and identify cell–cell interaction QTLs (cci-QTLs).

**Parameters:**
- `beta`: A numeric vector of length `k + 1` containing the estimated genetic effect sizes from bulk tissue and `k` cell types for the ligand–receptor pair of interest.
- `se`: A numeric vector of length `k + 1` containing the corresponding standard errors of the genetic effect size estimates.
- `weight`: A numeric vector of length `k` representing the relative contributions of the `k` cell types. The weights must sum to 1.

```r
gucci_cci <- refine.cci(beta, se, weight)
```

## Data Usage
- eQTLGen: https://www.eqtlgen.org
- OneK1K: https://onek1k.org/
- Ligand-Receptor pairs: https://github.com/LewisLabUCSD/Ligand-Receptor-Pairs
- Lymph Node spatial dataset: https://www.10xgenomics.com/datasets/human-lymph-node-1-standard-1-1-0
- Immune single-cell RNASeq datasets: https://cellxgene.cziscience.com/datasets


## Software Usage
- COLOC: https://chr1swallace.github.io/coloc/index.html
- Cell2location: https://github.com/BayraktarLab/cell2location
- Seurat: https://satijalab.org/seurat/
- ScDesign2: https://github.com/JSB-UCLA/scDesign2

## Contact

For questions or issues, please contact: ard5982@psu.edu
