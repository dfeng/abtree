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
 *  - use natire cpp data structures o/w
 * 
 */

// [[Rcpp::export]]
List rcpp_BuildTree(NumericVector y, NumericMatrix x,
                    IntegerVector trt, IntegerMatrix ordering,
                    IntegerVector ncat, int ntrt,
                    int min_bucket, int min_split, int max_depth) {
  // Rprintf("Start of something new\n");

  int ncol = x.ncol();
  int nrow = x.nrow();

  // Rprintf("Init the root\n");
  // Initialize the root
  NumericVector root_y(ntrt);
  IntegerVector root_n(ntrt);
  for (int i = 0; i < nrow; i++) {
    root_y[trt[i]] += y[i];
    root_n[trt[i]] ++;
  }
  Block b(root_y, root_n);

  Node *root = new Node();
  root->blok = b;

  // Rprintf("Time to partition\n");
  // TODO: do we want to do the ordering matrix in cpp?
  Partition(root, y, x, trt, ordering,
            ntrt, ncat,
            ncol, 0, nrow,
            min_bucket, min_split, max_depth,
            0);
  // Rprintf("Done partition\n");

  DoubleMat cp_table;
  // if no good splits were found even for the root
  // then skip cp_table
  if (root->split_col != -1) {
    cp_table = TreeComplexity(root);
  }

  DoubleMat tree_df;
  ExportTree(root, tree_df);

  // // Free the Ents!
  // DeleteTree(root);

  // Pointer to our cpp tree struct
  XPtr<Node> ptr(root);

  List ret;
  ret["cpp.tree"] = wrap(tree_df);
  ret["cpp.tree.ptr"] = ptr;
  ret["cp.table"] = wrap(cp_table);
  return ret;
}

// // [[Rcpp::export]]
// List rcpp_Prune(NumericMatrix tree_df,
//                 NumericVector valid_y, NumericMatrix valid_x,
//                 IntegerVector valid_trt, IntegerVector ncat,
//                 int ntrt, NumericMatrix cp_tbl) {
//   Node *root = ImportTree(tree_df, ntrt);
//   DoubleMat cp_table = NumToDoubleMat(cp_tbl);

//   PredictPrune(root, valid_y, valid_x, valid_trt, ncat, cp_table);

//   // find maximum cp_value
//   double argmax_cp;
//   double max_profit = -DBL_MAX;
//   for (int i = 0; i < cp_table.size(); i++) {
//     if (cp_table[i][1] >= max_profit) {
//       max_profit = cp_table[i][1];
//       argmax_cp = cp_table[i][0];
//     }
//   }
//   // Rprintf("argmax_cp %0.5f\n", argmax_cp);
//   // Rprintf("max_profit %0.5f\n", max_profit);
//   PruneTree(root, argmax_cp);

//   DoubleMat new_tree_df;
//   ExportTree(root, new_tree_df);
//   // return wrap(new_tree_df);

//   // Free the Ents!
//   DeleteTree(root);

//   List z;
//   z["tree"] = wrap(new_tree_df);
//   z["cp.table"] = wrap(cp_table);
//   return z;
// }

// // [[Rcpp::export]]
// List rcpp_Predict(NumericMatrix tree, std::vector<double> y,
//                   NumericMatrix x, std::vector<int> trt,
//                   std::vector<int> ncat) {
//   Node * root = ImportTree(tree);

//   IntegerVector pred_trt = Predict(root, y, x, trt, ncat);

//   DoubleMat testtree_df;
//   ExportTree(root, testtree_df, TRUE);

//   List z;
//   z["test"] = wrap(testtree_df);
//   z["trt"] = wrap(pred_trt);
//   return z;
// }