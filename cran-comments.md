The code in rtn1.cpp was taken from Jonathan Olmsted’s GitHub repository at https://github.com/olmjo/RcppTN/blob/master/src/rtn1.cpp, hence his cph designation.  All other code was hand-written by one of the authors or uses code that is part of other R / Rcpp packages.

This is a maintenance update to address CRAN check WARNINGS introduced by the recent upgrade of RcppArmadillo (Armadillo 15).
* Replaced deprecated `arma::is_finite(val)` calls with `std::isfinite(val)` in C++ source files (e.g., `three_utility_probit_helper_functions.cpp`), as required by the upstream RcppArmadillo change.
