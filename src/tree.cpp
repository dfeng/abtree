#include <Rcpp.h>

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
                    IntegerVector ncat,
                    int min_bucket, int min_split, int max_depth) {
  Rprintf("Start of something new\n");

  int ncol = x.ncol();
  int nrow = x.nrow();
  int ntrt = max(trt) + 1;

  Rprintf("Init the root\n");
  // Initialize the root
  NumericVector root_y(ntrt);
  IntegerVector root_n(ntrt);
  for (int i = 0; i < nrow; i++) {
    root_y[trt[i]] += y[i];
    root_n[trt[i]] ++;
  }
  Block b(root_y, root_n);

  // TODO: make sure that this is the right way of instantiating an object
  // I think so - basically dynamic memory is when you don't know the size of the object at compile time, which is not our case.
  // Node *root = new Node();
  // root->blok = b;
  Node root = Node();
  root.blok = b;

  Rprintf("Time to partition\n");
  // TODO: do we want to do the ordering matrix in cpp?
  Partition(&root, y, x, trt, ordering,
            ntrt, ncat,
            ncol, 0, nrow,
            min_bucket, min_split, max_depth,
            0);
  Rprintf("Done partition\n");
  // if no good splits were found even for the root
  // NumericMatrix cp_table;
  // if (root->split_col != -1) {
  //   cp_table = TreeComplexity(root);
  // }

  DoubleMat tree_df;
  ExportTree(&root, tree_df);

  List ret;
  ret["tree"] = wrap(tree_df);
  // ret["cp.table"] = wrap(cp_table);
  return ret;
}

// // [[Rcpp::export]]
// List rcpp_Prune(NumericMatrix tree_df, std::vector<double> y,
//                 NumericMatrix x_temp, std::vector<int> trt,
//                 std::vector<int> ncat,
//                 NumericMatrix cp_tbl) {
//   Node * root = ImportTree(tree_df);
//   DoubleMat cp_table = NumToDoubleMat(cp_tbl);

//   DoubleMat valid;
//   for (int j = 0; j < x_temp.ncol(); j++) {
//     DoubleVec xx;
//     for (int i = 0; i < x_temp.nrow(); i++) {
//       xx.push_back(x_temp(i,j));
//     }
//     valid.push_back(xx);
//   }

//   PredictPrune(root, y, valid, trt, ncat, cp_table);

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

// // NumericMatrix to DoubleMat
// DoubleMat NumToDoubleMat(NumericMatrix m) {
//   DoubleMat out;
//   for (int i = 0; i < m.nrow(); i++) {
//     DoubleVec o_row;
//     for (int j = 0; j < m.ncol(); j++) {
//       o_row.push_back(m(i,j));
//     }
//     out.push_back(o_row);
//   }
//   return out;
// }
