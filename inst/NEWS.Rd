\name{NEWS}
\title{News for package \pkg{scGate}}


\section{New in version 1.8.0}{
\itemize{
    \item Remove dependency from ggparty. Use self-contained ggplot2 function to visualize model trees.
    \item Update figshare link (changed url format)
    \item Fix kNN-smoothing being applied cumulatively across levels to signatures used at multiple levels; each level now smooths the original UCell scores. This may slightly change results at levels > 1
    \item Model levels are now evaluated in numerical order (level1, level2, ...) regardless of row order in the model table. Level names must be level1..levelN, without gaps, and are case-insensitive (e.g. 'Level4' is the same as 'level4')
    \item Stored UCell rankings (keep.ranks=TRUE) are only re-used if computed with the same assay, slot and maxRank; otherwise rankings are recomputed
    \item get_scGateDB: fix caching and force_update for tagged versions; remove system call to 'rm' (Windows compatibility)
    \item Fix performance.metrics with return_contingency=TRUE
    \item scGate with a single model no longer overwrites the scGate_multi metadata column
}}

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
