#include <Rcpp.h>
#include <float.h> // DBL_MAX

#include "vector.h" // shorthands for std::vector<type>
#include "node.h" // custom data types
#include "abtree.h" // master header file

using namespace Rcpp;

/*
 *
 * Idioms:
 *  - use Rcpp data structures whenever the object will be in contact with R
 *  - use native cpp data structures o/w
 * 
 */

// [[Rcpp::export]]
List rcpp_BuildTree(NumericVector y, NumericMatrix x,
                    IntegerVector trt, IntegerMatrix ordering,
                    IntegerVector ncat, int ntrt,
                    int min_bucket, int min_split, int max_depth,
                    int mtry) {

  int ncol = x.ncol();
  int nrow = x.nrow();

  // Initialize the root
  NumericVector root_y(ntrt);
  IntegerVector root_n(ntrt);
  for (int i = 0; i < nrow; i++) {
    root_y[trt[i]] += y[i];
    root_n[trt[i]] ++;
  }
  Block b(root_y, root_n);

  Node *root = new Node(ntrt);
  root->blok = b;

  // Rcout << "Partition" << std::endl;
  // TODO: do we want to do the ordering matrix in cpp?
  // Right now this isn't the bottleneck
  Partition(root, y, x, trt, ordering,
            ntrt, ncat,
            ncol, 0, nrow,
            min_bucket, min_split, max_depth,
            mtry,
            0);

  // Rcout << "CP" << std::endl;
  // DoubleMat cp_table;
  // // if no good splits were found even for the root
  // // then skip cp_table
  // if (root->split_col != -1) {
  //   cp_table = TreeComplexity(root);
  // }

  // Rcout << "Export" << std::endl;
  DoubleMat tree_df;
  ExportTree(root, tree_df);

  // Pointer to our cpp tree struct
  XPtr<Node, PreserveStorage, DeleteTree> xptr(root, true);

  // Rcout << "All Done" << std::endl;

  List ret;
  ret["cpp.tree"] = wrap(tree_df);
  ret["cpp.ptr"]  = xptr;
  // ret["cp.table"] = wrap(cp_table);
  return ret;
}

// [[Rcpp::export]]
List rcpp_Prune(SEXP xptr,
                NumericVector valid_y, NumericMatrix valid_x,
                IntegerVector valid_trt, IntegerVector ncat,
                NumericMatrix cp_table) {

  XPtr<Node, PreserveStorage, DeleteTree> root(xptr); // recover root

  PredictPrune(root, valid_y, valid_x, valid_trt, ncat, cp_table);

  // find maximum cp_value
  double argmax_cp;
  double max_profit = -DBL_MAX;
  for (int i = 0; i < cp_table.nrow(); i++) {
    if (cp_table(i,1) >= max_profit) {
      max_profit = cp_table(i,1);
      argmax_cp = cp_table(i,0);
    }
  }
  // Rcout << "argmax_cp" << argmax_cp << std::endl;
  // Rcout << "max_profit " << max_profit << std::endl;
  PruneTree(root, argmax_cp);

  DoubleMat new_tree_df;
  ExportTree(root, new_tree_df);

  List z;
  z["cpp.prune.tree"] = wrap(new_tree_df);
  z["cp.table"] = cp_table;
  return z;
}

// [[Rcpp::export]]
List rcpp_Predict(SEXP xptr,
                  NumericVector test_y, NumericMatrix test_x,
                  IntegerVector test_trt, IntegerVector ncat, int ntrt) {
  XPtr<Node, PreserveStorage, DeleteTree> root(xptr); // recover root
  // IntegerVector pred_trt = Predict(root, test_y, test_x, test_trt, ncat, ntrt);
  int nrow = test_x.nrow();
  IntegerVector pred_trt(nrow);
  NumericMatrix pred_prob(nrow, ntrt);
  // for each data point
  for (int i = 0; i < nrow; i++) {
    Node *leaf = PredictNode(root, test_x(i,_), ncat);
    pred_prob(i,_) = leaf->blok.p;
    // Rcout << "p: " << leaf->blok.p << std::endl;
    // Rcout << "n: " << leaf->blok.total_n << std::endl;
    pred_trt[i] = leaf->blok.opt_trt;
    // leaf->predict_y[test_trt[i]] += test_y[i];
    // leaf->predict_n[test_trt[i]]++;
  }
  // for each leaf
  // TODO: remind myself what this is for? don't think it does anything useful right now...
  // FillTest(root);

  // Pretty sure we don't need to export this, at least for now...?
  // DoubleMat test_tree_df;
  // ExportTree(root, test_tree_df);

  List z;
  // z["test"] = wrap(test_tree_df);
  z["predict.trt"] = wrap(pred_trt);
  z["predict.prob"] = wrap(pred_prob);
  return z;
}