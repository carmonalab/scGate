\name{NEWS}
\title{News for package \pkg{scGate}}

\section{New in version 1.7.0}{
\itemize{
    \item Update maintainer emails
    \item Bug fixes
}}

\section{New in version 1.6.2}{
\itemize{
    \item Enable parallelization over a list of scGate models (in multi-classifier mode)
    \item Fix bug on parallelization when specifying BPPARAM
    \item Added cell ontology links for controlled cell type dictionary
    \item Update accessor functions for compatibility with Seurat v5
}}

\section{New in version 1.4.0}{
\itemize{
    \item Use function \code{SmoothKNN()} from the \code{UCell} package for nearest neighbor smoothing of signature scores
    \item Clean up code for CRAN submission 
}}

\section{New in version 1.2.0}{
\itemize{
    \item Simultaneously evaluate multiple gating models, given as a list to model parameter. This returns individual gatings as well as a combined cell type annotation from all the models (scGate_multi metadata column), and allows running scGate as a multi-class classifier and annotate cell types in a dataset.
    \item New \code{reduction} parameter: this i) Allows nearest-neighbor smoothing on an existing dimensionality reduction (e.g. "pca" or "harmony") and ii) Can significantly speed up computation time and annotation of integrated objects
}}
