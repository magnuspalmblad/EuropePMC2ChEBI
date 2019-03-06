# EuropePMC2ChEBI
The main part of this project is the ChEBIplot R package and the KNIME workflow for extracting ChEBI annotations from a Europe PMC literature search. From these, one can for example look up the calculated masses and predicted octanol-water partition coefficients and visualize these. ChEBIplot provides functions for visualizing a single set of ChEBIs, comparisons of two sets in two colors, and tripartite comparisons as [RGB plots](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0102903).

The KNIME workflow should be configured with the search query, number of results pages to download and the output file name.

The original version that was used to generate the figures in the *Analytical Chemistry* emplyed for-loops that, though clear, are inefficient in R. These constructs should be fully vectorized or implemented in C++ using [Rcpp](http://www.rcpp.org/cpp "Rcpp's homepage").
