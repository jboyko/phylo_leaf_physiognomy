# DiLP preprocessing

`code/dilp_preprocessing.R` contains code from the `dilp` R package by
Matthew J. Butrim, Alexander J. Lowe, Andrew G. Flynn, Aly Baumgartner,
Daniel J. Peppe, Dana L. Royer, and James D. Boyko. Matthew J. Butrim is
identified as the copyright holder in the package metadata.

- Column names and specimen trait calculations are from `R/utils.R` and
  `R/dilp_functions.R` at the
  [fork commit used by the analysis](https://github.com/jboyko/dilp/tree/88b568ab910aa7a6494e3eb407f424c341099d71).
- The morphotype and site summaries use those same rules, arranged as in the
  verified development fork. The range warning checks for actual out-of-range
  values. The wrapper returns the three processed trait tables used in the
  analysis.
- The preprocessing source is available under GPL version 3 or later.
  The [license](DiLP-GPL-3.md) accompanies this copy.
