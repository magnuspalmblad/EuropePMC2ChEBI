# EuropePMC2ChEBI
The main parts of this project are the ChEBIplot R package and the KNIME workflow for extracting ChEBI annotations from a Europe PMC literature search and visualizing these as functions of mass and hydrophobicity. ChEBIplot provides functions for visualizing a single set of ChEBIs, comparisons of two sets in two colors, and tripartite comparisons as [RGB plots](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0102903).

The KNIME workflow should be configured with the search query, number of results pages to download and the output file name.

The original version that was used to generate the figures in the *Analytical Chemistry* [paper](https://pubs.acs.org/doi/10.1021/acs.analchem.8b05818) emplyed for-loops that, though clear, are inefficient in R. These constructs should probably be fully vectorized or implemented in C++ using [Rcpp](http://www.rcpp.org/cpp "Rcpp's homepage").

In the Examples folder, you can find data and the R code to visualize the applicability of atmospheric pressure chemical ionization (APCI), electrospray ionization (ESI) and electron impact/ionization (EI). There is also a PowerPoint template for presenting the visualizations produced by ChEBIplot.

If you are interested in using ChEBIplot, the KNIME workflow or contributing to this project, don't be shy!

-Magnus
