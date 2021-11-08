# scGate: marker-based purification of cell types from heterogeneous single-cell RNA-seq datasets

scGate is an R package that automatizes the typical manual marker-based approach to cell type annotation, to enable accurate and intuitive purification of a cell population of interest from a complex scRNA-seq dataset, without requiring reference gene expression profiles or training data. 

scGate builds upon our recent method [UCell](https://github.com/carmonalab/UCell) for robust single-cell signature scoring, and [Seurat](https://github.com/satijalab/seurat/), a comprehensive and powerful framework for single-cell omics analysis.

Briefly, scGate takes as input: i) a gene expression matrix stored in a Seurat object and ii) a “gating model” (GM), consisting of a set of marker genes that define the cell population of interest. The GM can be as simple as a single marker gene, or a combination of positive and negative markers. More complex GMs can be constructed in a hierarchical fashion, akin to gating strategies employed in flow cytometry. 
scGate evaluates the strength of the markers' gene expression in each cell using the rank-based method UCell, and then performs k-nearest neighbor (kNN) smoothing by calculating the mean UCell score across neighboring cells. kNN-smoothing aims at compensating for the large degree of sparsity in scRNA-seq data. Finally, a universal threshold over kNN-smoothed signature scores is applied in binary decision trees generated from the user-provided gating model, to annotate cells as either “pure” or “impure”, with respect to the cell population of interest. 

### Installation

```
install.packages("remotes")
remotes::install_github("carmonalab/UCell")
remotes::install_github("carmonalab/scGate")
```

### Testing the package

Use scGate to purify a cell population of interest using manually defined marker genes

```
library(scGate)

#Get a test scRNA-seq dataset (as a list of Seurat objects)
sample.data.seurat.list <- scGate::get_testing_data()

seurat_object <- sample.data.seurat.list$Satija

#Manually define a simple scGate's gating model to purify eg. Natural Killer (NK) cells, using a positive marker KLRD1 and negative marker CD3D
my_scGate_model <- gating_model(name = "NK", signature = c("KLRD1","CD3D-"))  

#scGate it!
seurat_object <- scGate(data = seurat_object, model = my_scGate_model)

#Use Seurat to visualize "Pure" and "Impure" cells
DimPlot(seurat_object, group.by = "is.pure")

#Use Seurat to subset pure cells
seurat_object_purified <- subset(seurat_object, subset = `is.pure` == "Pure" )
```

### Pre-defined Gating models

A database of gating models for scGate is available on [scGate_models](https://github.com/carmonalab/scGate_models) and can be loaded using ``get_scGateDB()`
```
#Get scGate's database of pre-defined gating models
scGate_models_DB <- get_scGateDB()

#For example, filter abT cells using one of scGate's pre-defined gating models
seurat_object <- scGate(seurat_object, model = models.DB$human$generic$Tcell.alphabeta)

DimPlot(seurat_object)
```
The first time you run `get_scGateDB()`  the database will be downloaded from the repository. On successive calls it will load your local version of the DB.

You may manually edit the available models (eg in Excel) or create new models for your cell type of interest. You can then load your custom model into R using:
```
my_scGate_model <- load_scGate_model("path_to_my.model")
```

You can use the `plot_tree` function to visualize the hierarchical structure of one of the models (requires [ggparty](https://cran.r-project.org/package=ggparty)).

```
install.packages("ggparty")
scGate::plot_tree(models.DB$human$generic$Tcell.alphabeta)
```

### Demos and tutorials

Check out the scGate's DEMO for a reproducible analysis, construction of hierarchical gating models, tools for performance evaluation and other advanced features

### References

scGate: marker-based purification of cell types from heterogeneous single-cell RNA-seq datasets. Massimo Andreatta, Ariel Berenstein, Santiago J. Carmona (in preparation)
