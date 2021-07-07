## scGate - filter specific cell types from single-cell datasets

scGate is a method to filter specific cell types from single-cell dataset.

To install scGate run:
```
remotes::install_github("carmonalab/UCell")
remotes::install_git("https://gitlab.unil.ch/carmona/scGate.git")
```

Pre-trained scGate models are available with the scGate package in: `scGate::scGate_DB`

For example, to filter T cells from a human dataset run:
```
query.data <- scGate(query.data, gating.model=scGate_DB$human$T.cell)
```




