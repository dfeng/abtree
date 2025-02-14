# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rcpp_BuildTree <- function(y, x, trt, ordering, ncat, ntrt, min_bucket, min_split, max_depth, mtry, split_cond) {
    .Call('_abtree_rcpp_BuildTree', PACKAGE = 'abtree', y, x, trt, ordering, ncat, ntrt, min_bucket, min_split, max_depth, mtry, split_cond)
}

rcpp_Prune <- function(xptr, valid_y, valid_x, valid_trt, ncat, cp_table) {
    .Call('_abtree_rcpp_Prune', PACKAGE = 'abtree', xptr, valid_y, valid_x, valid_trt, ncat, cp_table)
}

rcpp_Predict <- function(xptr, test_y, test_x, test_trt, ncat, ntrt, pred_max_depth) {
    .Call('_abtree_rcpp_Predict', PACKAGE = 'abtree', xptr, test_y, test_x, test_trt, ncat, ntrt, pred_max_depth)
}

