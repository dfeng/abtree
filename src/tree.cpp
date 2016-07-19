#include <Rcpp.h>
#include <vector>

#include "node.h"
#include "abtree.h"

using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_BuildTree(std::vector<double> y, NumericMatrix x_temp,
                    std::vector<int> trt, IntegerMatrix ordering_temp,
                    std::vector<int> ncat, int min_bucket,
                    int min_split, int max_depth) {
  int ncol = ncat.size();
  int nrow = y.size();

  DoubleMat x;
  IntMat ordering;
  // TODO: better way of handling *Matrix
  for (int j = 0; j < ncol; j++) {
    IntVec o;
    DoubleVec xx;
    for (int i = 0; i < nrow; i++) {
      o.push_back(ordering_temp(i,j));
      xx.push_back(x_temp(i,j));
    }
    ordering.push_back(o);
    x.push_back(xx);
  }
  // GetRNGstate();

  // Initialize the root
  double y0 = 0.0, y1 = 0.0;
  int n0 = 0, n1 = 0;
  for (int i = 0; i < nrow; i++) {
    if (trt[i] == 0) {
      y0 += y[i];
      n0++;
    } else {
      y1 += y[i];
      n1++;
    }
  }
  Block b(y0, y1, n0, n1);
  Node *root = new Node(b);
  Partition(root, y, x, trt, ordering, ncat,
            ncol, 0, nrow,
            min_bucket, min_split, max_depth, 0);

  // if no good splits were found even for the root
  DoubleMat cp_table;
  if (root->split_col != -1) {
    cp_table = TreeComplexity(root);
  }

  // Prune
  // Rprintf("%d\n", root->split_col);

  // Rprintf("cp_size: %d", cp_table.size());
  // PredictPrune(root, x, cp_table, ncol, nrow);

  // // find maximum cp_value
  // double argmax_cp;
  // double max_profit = -DBL_MAX;
  // for (int i = 0; i < cp_table.size(); i++) {
  //   if (cp_table[i][2] > max_profit) {
  //     max_profit = cp_table[i][2];
  //     argmax_cp = cp_table[i][1];
  //   }
  // }
  // // Rprintf("cp arg max: %0.2f", argmax_cp);
  // PruneTree(root, argmax_cp);
  // Rprintf("did we get here?\n");

  // PutRNGstate();

  DoubleMat tree_df;
  ExportTree(root, tree_df);

  List z;
  z["tree"] = wrap(tree_df);
  z["cp.table"] = wrap(cp_table);
  return z;
}

// [[Rcpp::export]]
List rcpp_Prune(NumericMatrix tree_df, std::vector<double> y,
                NumericMatrix x_temp, std::vector<int> trt,
                std::vector<int> ncat,
                NumericMatrix cp_tbl) {
  Node * root = ImportTree(tree_df);
  DoubleMat cp_table = NumToDoubleMat(cp_tbl);

  DoubleMat valid;
  for (int j = 0; j < x_temp.ncol(); j++) {
    DoubleVec xx;
    for (int i = 0; i < x_temp.nrow(); i++) {
      xx.push_back(x_temp(i,j));
    }
    valid.push_back(xx);
  }

  PredictPrune(root, y, valid, trt, ncat, cp_table);

  // find maximum cp_value
  double argmax_cp;
  double max_profit = -DBL_MAX;
  for (int i = 0; i < cp_table.size(); i++) {
    if (cp_table[i][1] >= max_profit) {
      max_profit = cp_table[i][1];
      argmax_cp = cp_table[i][0];
    }
  }
  // Rprintf("argmax_cp %0.5f\n", argmax_cp);
  // Rprintf("max_profit %0.5f\n", max_profit);
  PruneTree(root, argmax_cp);

  DoubleMat new_tree_df;
  ExportTree(root, new_tree_df);
  // return wrap(new_tree_df);
  List z;
  z["tree"] = wrap(new_tree_df);
  z["cp.table"] = wrap(cp_table);
  return z;
}

// [[Rcpp::export]]
List rcpp_Predict(NumericMatrix tree, std::vector<double> y,
                  NumericMatrix x, std::vector<int> trt,
                  std::vector<int> ncat) {
  Node * root = ImportTree(tree);

  IntegerVector pred_trt = Predict(root, y, x, trt, ncat);

  DoubleMat testtree_df;
  ExportTree(root, testtree_df, TRUE);

  List z;
  z["test"] = wrap(testtree_df);
  z["trt"] = wrap(pred_trt);
  return z;
}

// NumericMatrix to DoubleMat
DoubleMat NumToDoubleMat(NumericMatrix m) {
  DoubleMat out;
  for (int i = 0; i < m.nrow(); i++) {
    DoubleVec o_row;
    for (int j = 0; j < m.ncol(); j++) {
      o_row.push_back(m(i,j));
    }
    out.push_back(o_row);
  }
  return out;
}
