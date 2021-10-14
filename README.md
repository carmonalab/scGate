# scGate

scGate is a method to filter specific cell types from single-cell dataset.

Several methods exist for multiclass classification of single cells based on profiles learned from a training set (e.g. singleR). However, the performance of these methods highly depend on the quality, annotation, and completeness of the training set. They can only detect the cell types present in the training set and are difficult to customize (i.e. retraining is necessary). With scGate, we aim at providing a training-free method for detecting specific cell types, using signature-based, hierachical models that can be easily customized by the user. 


### Installation
```
remotes::install_github("carmonalab/UCell")
remotes::install_github("carmonalab/scGate")
install.packages("ggparty")
```

### Testing the package

Obtain a single-cell dataset for testing (**UPDATE LINK**) and run scGate to filter a specific cell type
```
test.set <- readRDS(ADD_LINK)

#Load scGate and upload a database of models
library(scGate)
models.DB <- get_scGateDB()

#For example, filter B cells from this dataset
model.Bcell <- models.DB$human$generic$PanBcell 
test.set <- scGate(test.set, model = model.Bcell)
DimPlot(test.set)
```

### Gating models
A database of gating models for scGate is available on [GitHub](https://github.com/carmonalab/scGate_models), and can loaded directly into your R workspace with the following function:
```
library(scGate)
models.DB <- get_scGateDB()
```
The first time you run this command the scGate database will be downloaded from the repository. On successive calls it will load your local version of the DB.

You may freely edit the available models or create new models for your cell type of interest. You can then load your custom model into R using:
```
my.model <- load_scGate_model("path_to_my.model")
```

You can also apply the `plot_tree` function to visualize the hierarchical structure of one of the models.
```
scGate::plot.tree(scGate::plot_tree(models.DB$human$generic$Tcell))
scGate::plot.tree(scGate::plot_tree(my.model)
```

### Demos and tutorials

Add this section

### References

Add references



