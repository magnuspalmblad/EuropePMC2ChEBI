There are a number of relatively easy things that could and should be done next:

1. Extract document frequencies from all publication years (from ca. 1800 until 2020) to provide all data for TFIDF and allow time-matched TFIDF normalization.
2. Move more hard-coded constants to function header, defining sensible defaults.
3. Accelerate the plotting by vectorizing loop constructs or refactoring part of the code in C++ using Rcpp.
4. Match KNIME workflow output to ChEBI plot input (additional filters).
5. Interactive HTML5 visualization with tooltip revealing and linking to most abundant ChEBI entry in (log P/mass) bin.
