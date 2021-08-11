# scGate

scGate is a method to filter specific cell types from single-cell dataset.

Several methods exist for multiclass classification of single cells (e.g. singleR); but they do not perform well.
Aiming to achieve high resolution for detecting one cell type (or a subset of cell types), we propose to simplify 
multiclass classifiction to a two-class classification, i.e. the cell type(s) of interest vs. everything else.
Details of the algorithm (TO DO).

### Installation
```
remotes::install_github("carmonalab/UCell")
remotes::install_github("carmonalab/scGate")
```

### Functions
scGate is composed of two main functions:
* `train_scGate`: given a set of signatures and a reference set, calculate expected mean and SD for all cell types 
in the reference set. Deviations from this expected distribution (Z-score) can then be used to gate for specific cell types 
in any query dataset
* `scGate`: apply a scGate gating model to filter specific cell types in a query dataset

See the R help for all the parameters.

### Gating models
Pre-trained scGate models are available with the scGate package in: `scGate::scGate_DB`

For example, to filter T cells from a human dataset run:
```
query.data <- scGate(query.data, gating.model=scGate_DB$human$CD8.Tcell)
```

Pure/Impure classifications are stored in the object metadata (column `is.pure`), together with the most likely annotation for gated-out cells (column `scGate.annotation`).


